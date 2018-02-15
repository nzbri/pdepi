from collections import defaultdict, OrderedDict
import dateutil.parser
import datetime
import operator
import sys,traceback
import csv

import pharmacdata

class MOHData:
    """ Process the raw MOH data and output into a csv file useful for diagnoses """
    
    def __init__(self):
        
        nhi_all =set()
        nhi_pd = set()
        nhi_pharmac=set()
        
        self.pharms = pharmacdata.PharmacData()
        
        pharmac_missing_mortality=defaultdict(int)
        pharmac_missing_admission=defaultdict(int)

        fields =OrderedDict([('nhi',1),
                             ('age',2),
                             ('year',3),
                             ('sex',4),
                             ('ethnicity',5),
                             ('dhb',6),
                             ('source',7)
                             ])
        
        f_out_m = open('output/moh_missing_in_pharms.csv',"w")
        dwm = csv.DictWriter(f_out_m, delimiter=',',restval='NA',fieldnames=fields)
        dwm.writeheader()
        
        # Read in NHIs from pharmac data
        fname = 'output/classification.csv'
        with open(fname, "r") as f:
            records = csv.DictReader(f)
            
            for record in records:
                #print record
                if record['year_in_data']=='1':
                    nhi_pharmac.add(record['nhi'])
        
        
        ## Process mortality data
        
        deceased_count = 0
        pd_deceased_count = 0
        no_pd_mortality = set()
        
        files = ('raw/mos3358all/mos3358.csv',
                 'raw/mos3464/mos3464.csv')
        for fname in files:
            with open(fname, "r") as f:
                records = csv.DictReader(f)
                
                for record in records:
                    deceased_count +=1
                    pd = False
                    try:
                        nhi = record['MAST_NHI']
                    except KeyError:
                        nhi = record['PRIM_HCU']
                        
                    for field in ('icda', #Underlying cause of death
                                  'icdd', #Underlying cause of death
                                  'icdf1','icdf2','icdf3','icdf4', #Other relevant diseases present
                                  'icdg1','icdg2', #Other contributing causes
                                  'icdc1','icdc2','icdj1','icdj2' #Cancer as non-contributing cause
                                  ):
                        nhi_all.add(nhi)
                        
                        try:
                            value = record[field]
                        except: 
                            continue
                        
                        if 'G20' in value:
                            pd = True
                            nhi_pd.add(nhi)
                            pd_deceased_count+=1
                            if nhi not in nhi_pharmac:
                                pharmac_missing_mortality[record['REGYR']]+=1
                            #print "Parkinson's found: {}".format(record['MAST_NHI'])
                            else:
                                print "Parkinson's in mortality and Pharmac: {} {}".format(nhi,record['DOD'])
                            if nhi not in nhi_pharmac:
                                dhb = self.pharms.map_item(record['DHBDOM'],self.pharms.dhb_mapping)
                                date = dateutil.parser.parse(record['DOD'],dayfirst=True)
                                dwm.writerow({'age':record['AGE_AT_DEATH_YRS'],
                                              'year':date.strftime("%Y"),
                                              'nhi':nhi,
                                              'sex':record['SEX'],
                                              'source':'Mortality',
                                              'dhb':dhb})
                            
                        
                    if not pd:
                        no_pd_mortality.add(nhi)
            
        print "Number deceased with PD: {} from a total of {} records".format(pd_deceased_count,deceased_count)


        ## Process admission data

        diagnoses=defaultdict(list)
        admission_count = 0
        pd_not_noted_on_death_count = set()
        
        diagfields = ['diag{:02.0f}'.format(i) for i in xrange(1,31)]
        
        fname = 'raw/pus9058all/pus9058.csv'
        with open(fname, "r") as f:
            records = csv.DictReader(f)
            
            for record in records:
                admission_count +=1
                nhi_all.add(record['MAST_NHI'])
                for field in diagfields:
                    diagnosis = record[field]
                    if diagnosis:
                        diagnoses[diagnosis].append(record['MAST_NHI'])
                    if diagnosis == 'G20':
                        if record['MAST_NHI'] not in nhi_pharmac:
                            admission_date = dateutil.parser.parse(record['EVSTDATE'])
                            pharmac_missing_admission[admission_date.year]+=1
                            date = dateutil.parser.parse(record['EVSTDATE'],dayfirst=True)
                            dwm.writerow({'age':record['AGE_DSCH'],
                                          'year':date.strftime("%Y"),
                                          'nhi':record['MAST_NHI'],
                                          'sex':record['GENDER'],
                                          'dhb':self.pharms.map_item(record['DHBDOM'],self.pharms.dhb_mapping),
                                          'source':'Admissions'})
                        if record['MAST_NHI'] in no_pd_mortality:
                            pd_not_noted_on_death_count.add(record['MAST_NHI'])

                        
        print "Number admissions with PD: {} total from {} unique individuals (total of {} admissions)".format(len(diagnoses['G20']),
                                                                                                               len(set(diagnoses['G20'])),
                                                                                                               admission_count)
        with open("output/admission_diagnoses.csv", "w") as f:
            
            for code in sorted(diagnoses.keys()):
                output= "{},{},{}\n".format(code,
                                              len(diagnoses[code]),
                                              len(set(diagnoses[code]))
                                              )
                f.write(output)
        
        nhi_pd = nhi_pd | set(diagnoses['G20'])
        
        print "Number of unique PD identified from mortality/admissions: {}".format(len(nhi_pd))
        
        
        
        print "Total NHI in pharmac data: {}".format(len(nhi_pharmac))
        print "Number of PD (identified from mortality/admissions) in pharmac: {}".format(len(nhi_pharmac & nhi_pd))
        print "Number of other diagnoses (identified from mortality/admissions) in pharmac: {}".format(len(nhi_pharmac & (nhi_all - nhi_pd)))
        print "Total PD (identified from mortality/admissions) not in pharmac: {}".format(len(nhi_pd-nhi_pharmac))
        
        print "Admission data shows PD, has died, but PD not shown in mortality data: {}".format(len(pd_not_noted_on_death_count))
        
        
        print "Missing in pharmac but in mortality by year"
        print pharmac_missing_mortality
        print "Missing in pharmac but in admission by year"
        print pharmac_missing_admission
        
        ## Write out diagnoses to file
        
        fields =OrderedDict([('nhi',1),
                             ('diagnosis',2),
                             ('ethnicity',3),
                             ])
        
        f_out = open('output/moh_diagnoses.csv',"w")
        dwd = csv.DictWriter(f_out, delimiter=',',restval='NA',fieldnames=fields)
        dwd.writeheader()
        
        for nhi in nhi_all:
            
            data = {'nhi':nhi}
            
            if nhi in nhi_pd:
                data['diagnosis'] = 'PD'
            else:
                data['diagnosis'] = 'Other'
            
            dwd.writerow(data)
        
if __name__ == '__main__':

    mortality = MOHData()
 