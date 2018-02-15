#!/usr/bin/env python
from collections import defaultdict, OrderedDict
import dateutil.parser
import datetime
import operator
import sys,traceback
import csv
import numpy

import diagnoses
#from __builtin__ import None


class Dispensing:
    """Details of a dispensing"""
    
    def __init__(self, date, days, dose = None):
        """(self, string, int, float) -> None"""
        
        self.date = dateutil.parser.parse(date,dayfirst=True)
        self.days = days
        self.dose = dose
        self.last_date = self.date + datetime.timedelta(days=self.days)
        
    def __lt__(self, other):
        """Date -> Boolean"""
        return self.date < other.date
    
    def __repr__(self):
        return "{}:{}".format(self.date,self.days)

class Dispensings:
    """Methods to add and summarise dispensings that a particular individual has had"""
    
    def __init__(self, nhi, age=None, sex = None, birthdate = None,
                 continuityFile=None, classificationFile=None,
                 diagnosis="Empty",local_diagnosis="Empty",moh_diagnosis="Empty",
                 medical_registrar=None):
        
        self.nhi = nhi
        self.age = age
        self.sex = sex
        self.birthdate = birthdate
        self.date_of_death = None
        self.ethnicity = defaultdict(int)
        self.dhb = defaultdict(int)
        self.nyears = 0
        self.dispensings = defaultdict(list)
        self.durations = defaultdict(list)
        self.fOutContinuity = continuityFile
        self.fOutClassification = classificationFile
        self.diagnosis = diagnosis
        self.local_diagnosis = local_diagnosis
        self.moh_diagnosis = moh_diagnosis
        self.medical_registar = medical_registrar
        
        self.first_year = None # Calculated and set by drugs_received_by_year (called from classify)
        self.first_month = 12
        self.final_classification = None
        
        self.l_dopa = ['Sinemet','Madopar','Sindopa','Kinson']
        self.da_agonist = ['Lisuride', 'Pergolide', 'Ropinirole', 'Bromocriptine', 
                              'Apomorphine', 'Pramipexole']
    
    def age_at_year(self,year):
        
        at_date = datetime.datetime(year,12,1)
        dob = dateutil.parser.parse(self.birthdate,dayfirst=True)
        return (at_date-dob).days/365.0
    
    def primary_dhb(self):
        count = 0
        dhb = None
        for item in self.dhb:
            if self.dhb[item] > count:
                count = self.dhb[item]
                dhb = item
        return dhb

    def primary_ethnicity(self):
        count = 0
        ethnicity = 'Unknown'
        for item in self.ethnicity:
            if (item != 'Unknown') and (self.ethnicity[item] > count):
                count = self.ethnicity[item]
                ethnicity = item
        return ethnicity    
    
    def add_dispensing(self,drug,date,days,dose=None):
        self.dispensings[drug].append(Dispensing(date,days,dose))
    
    def years_receieved_drugs(self):
        years = set()

        for drug in self.dispensings:
            for dispensing in self.dispensings[drug]:
                years.add(dispensing.date.year)
        
        return years      
    
    def days_unmedicated_before_death(self):
        
        if self.date_of_death in (None,''):
            return 'NA'
        else:
            longer_than_human_lifespan = 2000000
            days_unmedicated = longer_than_human_lifespan
            
            for drug in self.dispensings:
                if (drug in self.l_dopa) or (drug in self.da_agonist):
                    for dispensing in self.dispensings[drug]:
                        diff = (dateutil.parser.parse(self.date_of_death,dayfirst=True) - dispensing.last_date).days
                        days_unmedicated = min(days_unmedicated,diff)
                        
            if days_unmedicated == longer_than_human_lifespan:
                return 'NA'
            else:
                return days_unmedicated
    
    def drugs_received_by_year(self):
        
        drugs_received_by_year=defaultdict(lambda: defaultdict(int))
        
        for drug in self.dispensings:
            for dispensing in self.dispensings[drug]:
                drugs_received_by_year[dispensing.date.year][drug]+=dispensing.days
        
        
        years = sorted(drugs_received_by_year.keys())
        self.first_year = years[0]
        
        # Determine first month of received drug in year
        for drug in self.dispensings:
            for dispensing in self.dispensings[drug]:
                if dispensing.date.year == self.first_year and dispensing.date.month < self.first_month:
                    self.first_month = dispensing.date.month
        
        nyears = len(years)
        
        for i in xrange(nyears-1,0,-1):
            for j in xrange(0,i):
                #print i, j
                for drug in drugs_received_by_year[years[j]]:
                    drugs_received_by_year[years[i]][drug] += drugs_received_by_year[years[j]][drug]
            
            
        
    
    def max_dose(self,drug,year=None):
        
        max_dose = None
        for dispensing in self.dispensings[drug]:
            if year:
                if dispensing.date.year == year and dispensing.dose > max_dose:
                    max_dose = dispensing.dose
            else:
                if dispensing.dose > max_dose:
                    max_dose = dispensing.dose
                
        return max_dose       
        #if not max_dose:
        #    return 0
        #else:
        #    return max_dose
    
    def drugs_received(self):
        """ Return dictionary of drugs received as keys with days as the value """
        
        drugs_received=defaultdict(int)
        for drug in self.dispensings:
            for dispensing in self.dispensings[drug]:
                drugs_received[drug]+=dispensing.days
        return drugs_received
    
    def write_presciription_block_line(self, drug, block, Ndispensings, first_dispensing, last_dispensing):

        duration = (last_dispensing.date - first_dispensing.date).days + last_dispensing.days

        data={'nhi':self.nhi,
              'ethnicity':self.primary_ethnicity(),
              'drug':drug,
              'start_date':str(first_dispensing.date),
              'duration':duration,
              'block':block,
              'dispensings':Ndispensings}
        
        self.fOutContinuity.writerow(data)
        
        
        self.durations[drug].append(duration)
        
    def process_dispensings(self):
        ## Maximum time between prescriptions before we consider it a break
        maximum_time_between_prescriptions = datetime.timedelta(weeks=20)
        
        for drug in self.dispensings:
            block = 0
            Ndispensings = 0
            
            ## Sort dispensings
            self.dispensings[drug].sort()
            dispensings = self.dispensings[drug]

            first_dispensing = dispensings[0]
            previous_dispensing = dispensings[0]
            for dispensing in dispensings:
                if dispensing.date > (previous_dispensing.date + maximum_time_between_prescriptions):
                    # Have had a break and now back on
                    self.write_presciription_block_line(drug,
                                                        block,
                                                        Ndispensings,
                                                        first_dispensing,
                                                        previous_dispensing)
                    block +=1
                    Ndispensings = 1
                    first_dispensing=dispensing
                else:
                    Ndispensings +=1
                previous_dispensing = dispensing
                
            # Write out last/only prescription of drug
            self.write_presciription_block_line(drug, block, Ndispensings, first_dispensing, previous_dispensing)
    
    def classify(self,by_year=False):


        if by_year:
            
            sorted_years = sorted(self.years_receieved_drugs())
            nyears = len(sorted_years)
            first_year = sorted_years[0]
            last_year = sorted_years[-1]
            
            drugs_received = self.drugs_received()
            classification, subclassification, dose = self.classify_worker(drugs_received)
            year_in_data = 0
            
            self.final_classification = classification
            
            for year in sorted_years:
                
                year_in_data += 1
                
                # If this is the last year of data explain why
                if year == last_year:
                    if self.date_of_death:
                        if dateutil.parser.parse(self.date_of_death,dayfirst=True).year <= year+1:
                            future_status = 'LAST_YEAR_DECEASED'
                        else:
                            future_status = 'LAST_YEAR_UNKNOWN'
                    else:
                        if year == 2014:
                            future_status = 'NO_MORE_FOLLOWUP'
                        else:
                            future_status = 'LAST_YEAR_UNKNOWN'
                else:
                    if year+1 in sorted_years:
                        future_status = 'PRESENT_NEXT_YEAR'
                    else:
                        future_status = 'MISSING_BUT_RETURN'
                
                data={'nhi': self.nhi,
                     'age': "{:.1f}".format(self.age_at_year(year)),
                     'sex': self.sex,
                     'ethnicity': self.primary_ethnicity(),
                     'dhb': self.primary_dhb(),
                     'year': year,
                     'classification': classification,
                     'subclassification': "{}-{}".format(classification,subclassification),
                     'years_of_data': nyears,
                     'year_in_data': year_in_data,
                     'year_first_seen':first_year,
                     'age_first_seen':"{:.1f}".format(self.age_at_year(first_year)),
                     'year_last_seen':last_year,
                     'diagnosis':self.diagnosis,
                     'local_diagnosis':self.local_diagnosis,
                     'moh_diagnosis':self.moh_diagnosis,
                     'dose':dose,
                     'ldopa_on_days':self.ldopa_days,
                     'ldopa_first_to_last_days':self.ldopa_period,
                     'days_unmedicated_before_death':self.days_unmedicated_before_death(),
                     'future_status':future_status
                }

                
                self.fOutClassification.writerow(data)
                
    def classify_worker(self,drugs_received,year=None):
        pd_only_drugs = ['Apomorphine','Pergolide','Tolcapone','Entacapone']
        da_agonist = ['Bromocriptine','Lisuride','Pramipexole']
        anticholinergic=['Benztropine','Procyclidine','Orphenadrine']
        
        ropinirole_dose = self.max_dose('Ropinirole',year)
        pramipexole_dose = self.max_dose('Pramipexole',year)
        
        
        on_anticholinergic = False
        anticholinergic_duration = 0
        for drug in anticholinergic:
            if drug in drugs_received:
                on_anticholinergic = True
                anticholinergic_duration += drugs_received[drug]
        
        on_ldopa = False
        ldopa_duration = 0
        for drug in self.l_dopa:
            if drug in drugs_received:
                on_ldopa = True
                ldopa_duration = max(ldopa_duration, max(self.durations[drug]))
        
        if on_ldopa:

            # Determine first and last date on ldopa
            first_date = datetime.datetime(2100,1,1)
            last_date = datetime.datetime(1900,1,1)
                        
            for drug in self.l_dopa:
                if drug in drugs_received:
                    
                    for dispensing in self.dispensings[drug]:
                        
                        if dispensing.date < first_date:
                            first_date = dispensing.date
                            
                        final_date = dispensing.date + datetime.timedelta(days=dispensing.days)
                        if final_date > last_date:
                            last_date = final_date
            
            self.ldopa_period = (last_date-first_date).days
            
            # Determine days on drug (can overlap so not a simple case of adding days together)
            if self.ldopa_period > 0:
                
                # Create an array of zeros over the period on ldopa
                days_over_period = numpy.zeros(self.ldopa_period)
                
                # Fill in days where on ldopa
                for drug in self.l_dopa:
                    
                    if drug in drugs_received:
                    
                        for dispensing in self.dispensings[drug]:
                            start_day = (dispensing.date-first_date).days
                            finish_day = start_day + dispensing.days
                            days_over_period[start_day:finish_day]=1
                
                # add together all days where on ldopa
                self.ldopa_days = sum(days_over_period)

            if self.ldopa_period == 0 or self.ldopa_days == 0:
                # Missing days data
                self.ldopa_period = 'NA'
                self.ldopa_days = 'NA'
            
        else:
            # Not on ldopa
            self.ldopa_period = 'NA'
            self.ldopa_days = 'NA'


                
        
        
        on_agonist = False
        agonist_duration = 0
        for drug in da_agonist:
            if drug in drugs_received:
                on_agonist = True
                agonist_duration = max(agonist_duration, max(self.durations[drug]))

        ### Consider Definite Cases

        for drug in pd_only_drugs:
            if drug in drugs_received:
                return "Very probable",drug, "NA"

        if on_ldopa:
            if 'Amantadine' in drugs_received:
                return "Very probable","L-Dopa & Amantadine", "NA"
            if 'Selegiline' in drugs_received:
                return "Very probable","L-Dopa & Selegiline", "NA"
            if on_agonist:
                return "Very probable","L-Dopa & other DA agonist", "NA"
            if on_anticholinergic:
                return "Very probable","L-Dopa & Anticholinergic", "NA"
        
        
        if 'Pramipexole' in drugs_received and pramipexole_dose >= 0.75:
            return "Very probable","Pramipexole >= 0.75 mg/day", pramipexole_dose
            
        if on_agonist or ('Ropinirole' in drugs_received):
            if 'Amantadine' in drugs_received:
                return "Very probable","Agonist & Amantadine", "NA"
            if on_anticholinergic:
                return "Very probable","Agonist & Anticholinergic", "NA"
            
        
        ### Consider Probable Cases
        
        if 'Selegiline' in drugs_received:
            return "Probable","Selegiline", "NA"
        
        if on_ldopa and ('Ropinirole' in drugs_received):
            if ropinirole_dose > 0.6:
                return "Probable","L-Dopa & Ropinirole > 0.6 mg/day", "NA"
            else:
                return "Possible","L-Dopa & Ropinirole <= 0.6 mg/day or unknown", "NA"
                
        if on_ldopa:
            if ldopa_duration > 180:
                return "Probable","L-DopaOnly more than 180 days", "NA"
            else:
                return "Possible","L-DopaOnly less than 180 days", "NA"
        
        if 'Lisuride' in drugs_received:
            return "Possible","Lisuride", "NA"
        
        if 'Pramipexole' in drugs_received and pramipexole_dose < 0.75:
            return "Possible", "Pramipexole < 0.75 mg/day or unknown dose", pramipexole_dose
        
        
        ### Consider Possible Cases
        if 'Ropinirole' in drugs_received:            
            if max(self.durations['Ropinirole']) > 180:
                return "Possible", "RopiniroleOnly more than 180 days", "NA"
            else:
                return "Unlikely", "RopiniroleOnly less than 180 days", "NA"
        
        if 'Amantadine' in drugs_received:
            return "Unlikely", "Amantadine Only", "NA"
            
        if ('Bromocriptine' in drugs_received):
            if ('Ropinirole' in drugs_received):
                return "Unlikely", "Bromocriptine and Ropinirole", "NA"
            else:
                return "Unlikely", "Bromocriptine Only", "NA"
        
        if on_anticholinergic:
            return "Unlikely","Anticholineric Only", "NA"
        
        
        print "Not classified: ", drugs_received
        return "Not classified", "None", "NA"
        
def process_prescriptions_csv(inFile,outContinuity,outClassification,inDiagnoses,inMohDiagnoses):
    
    #Continuity of drugs
    fOutContinuity = open(outContinuity, "w")
    fields = [('nhi',1), 
              ('ethnicity',2), 
              ('drug',3), 
              ('start_date',4),('duration',5),('block',6),('dispensings',7)]
    dwcont = csv.DictWriter(fOutContinuity,delimiter=',',restval='NA',fieldnames=OrderedDict(fields))
    dwcont.writeheader()
    
    fOutClassification = open(outClassification, "w")
    fields = [('nhi',1), ('age',2), ('sex',3), ('ethnicity',4), ('dhb',5), 
              ('year',6), ('classification',7), ('subclassification',8),
              ('years_of_data',9), ('year_in_data',10), ('year_first_seen',11),
              ('age_first_seen',12), ('year_last_seen',13), ('diagnosis',14), 
              ('local_diagnosis',14),('moh_diagnosis',14),('dose',15),
              ('ldopa_on_days',16), ('ldopa_first_to_last_days',17), 
              ('days_unmedicated_before_death',18),('future_status',19)]
    dwclass = csv.DictWriter(fOutClassification,delimiter=',',restval='NA',fieldnames=OrderedDict(fields))
    dwclass.writeheader()
    
    # Summary of providers
    fOutProviders = open(outProviders, "w")
    fields =[('nhi',None),('nproviders',None),('dhb',None)]
    dwp = csv.DictWriter(fOutProviders, delimiter=',',restval='NA',fieldnames=OrderedDict(fields))
    dwp.writeheader()
    
    # Incidence
    fOutIncidence = open("output/incidence.csv", "w")
    fields =[('nhi',None),
             ('age',None),
             ('year',None),
             ('month',None),
             ('classification',None)]
    ordered_fieldnames = OrderedDict(fields)
    dwi = csv.DictWriter(fOutIncidence, delimiter=',',restval='NA',fieldnames=ordered_fieldnames)
    dwi.writeheader()
        
    all_diagnoses = diagnoses.Diagnoses(inDiagnoses,inMohDiagnoses)
    providers = Providers(inMedicalCouncil)
    
    
    with open(inFile, "r") as f:        
        records = csv.DictReader(f)
        
        previous_nhi=None
        for record in records:
                
            nhi = record['nhi']
            
            ## If this nhi is new we need to process dispsensing for previous nhi
            ## Then create an empty dispensings object for new nhi
            ## Corner case in that for first line the nhi will be different but 
            ## their will be no previous dispensings to process
            if previous_nhi != nhi:
                
                if previous_nhi != None:
                    ## Process data collected
                    dispensings.process_dispensings()
                    dispensings.classify(by_year=True)
                    #everyone.append(dispensings)
                    #print "DHB {}, Providers {}, Ethnicity {}".format(dispensings.dhb,dispensings.provider,dispensings.ethnicity)
                    dwp.writerow({'nhi':nhi,
                                  'dhb':dispensings.primary_dhb(),
                                  #'nproviders':dispensings.total_number_providers()
                                  })
                    
                    dwi.writerow({'nhi':nhi,
                                  'age':dispensings.age,
                                  'year':dispensings.first_year,
                                  'month':dispensings.first_month,
                                  'classification':dispensings.final_classification})
                    #dispensings.check_for_unknown_providers()
                    
                ## Setup for new records
                age = float(record['age'])
                sex = record['sex']
                
                # Use CDHB/Clinic diagnoses as default, if don't have use MoH diagnoses
                diagnosis = all_diagnoses.getDiagnosis(nhi)
                local_diagnosis = all_diagnoses.getLocalDiagnosis(nhi)
                moh_diagnosis = all_diagnoses.getMohDiagnosis(nhi)
                
                dispensings = Dispensings(nhi,
                                          age,
                                          sex,
                                          record['birthdate'],
                                          dwcont,
                                          dwclass,
                                          diagnosis,
                                          local_diagnosis,
                                          moh_diagnosis,
                                          providers)
            
            previous_nhi=nhi
            
            dispensings.ethnicity[record['ethnicity']] += 1
            dispensings.dhb[record['dhb']] += 1
            
            if record['date_of_death'] != 'NA':
                dispensings.date_of_death = record['date_of_death']
            
            if 'NA' in record['dose_mg']:
                dose = None
            else:
                dose = float(record['dose_mg'])
                
            if record['days_supply'] == 'NA':
                days = 0
            else:
                days = int(record['days_supply'])
                
            ## Add dispensing
            dispensings.add_dispensing(drug = record['drug'].replace('"',''),
                                       date = record['date'],
                                       days = days,
                                       dose = dose)
        
        fOutContinuity.close()
        fOutClassification.close()
        
        print "Unknown IDs: {}".format(providers.number_unknown())



if __name__ == '__main__':
    inFile = "output/included_records.csv"
    inDiagnoses = "input/diagnoses_all_sources.csv" 
    inMohDiagnoses = "output/moh_diagnoses.csv" 
    
    outContinuity = "output/continuity.csv"
    outClassification = "output/classification.csv"
    outProviders = "output/providers.csv"

    
    process_prescriptions_csv(inFile,outContinuity,outClassification,
                                inDiagnoses,inMohDiagnoses)
