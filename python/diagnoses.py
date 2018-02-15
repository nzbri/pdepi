import sys,traceback
import csv
from collections import defaultdict

class Diagnoses:
    def __init__(self, 
                 local_diagnoses_filename = "input/diagnoses_all_sources.csv", 
                 #local_diagnoses_filename = "", 
                 moh_diagnoses_filename = "output/moh_diagnoses.csv"):
        
        self.all_diagnoses = defaultdict(list)
        
        ### Import local diagnoses
        
        self.local_diagnoses = dict()
                    
        
        ### Import moh diagnoses
        
        self.moh_diagnoses = dict()
        
        try:
            inputFile = open(moh_diagnoses_filename)
            reader = csv.DictReader(inputFile)
        except:
            traceback.print_exc(file=sys.stdout)
            print("Unable to open {0}, all MOH diagnoses will be empty".format(moh_diagnoses_filename))
            return
            
        for row in reader:
            self.moh_diagnoses[row['nhi']]=row['diagnosis']
            self.all_diagnoses[row['nhi']].append(row['diagnosis'])
        
        
        
        def process_diagnosis_file(name="", filename="", nhi_field_name="",
                                   diagnosis_field_name="",diagnosis_detail=""):
            
            print("Processing updated {} file \n".format(name))
            
            new = 0
            diags = dict()
            
            for row in csv.DictReader(open(filename)):
                
                nhi = row[nhi_field_name].replace(" ","")
                diag = row[diagnosis_field_name].replace(" ","")
                
                diags[nhi]=diag
                
                if diag == 'other':
                    diag = "Other"
                elif diag == 'MH':
                    diag = "Other"
                elif diag == 'unknown':
                    continue
                
                # ignore Neurology database others with no additional detail, unreliable
                if name == 'Neurology database' and diag == 'Other' and row[diagnosis_detail]=='':
                    continue
                
                try:
                    current_diagnosis = self.local_diagnoses[nhi]
                    
                    self.all_diagnoses[nhi].append(diag)
                    
                    if current_diagnosis == 'NA':
                        current_diagnosis = diag
                    
                    if name in ('Research database'):
                        current_diagnosis = diag
                    
                    if current_diagnosis != diag:
                        print "{}: Current/{} diagnosis different: {}/{}".format(nhi,name,current_diagnosis,diag)
                except KeyError:
                    print "{} New known diagnosis from {}: {}".format(name,nhi,diag)
                    self.all_diagnoses[nhi].append(diag)
                    self.local_diagnoses[nhi] = diag
                    new +=1
            
            print("\n{} new diagnoses found from {}\n".format(new,name))
            
            return diags

        rd_diags = process_diagnosis_file(name="Research database",
                               filename='input/diagnoses_alice_2016.csv',
                               nhi_field_name='NHI',
                               diagnosis_field_name = 'DiseaseGroup')     

        process_diagnosis_file(name="Tim PP",
                               filename='input/diagnoses_tim_pp_2015.csv',
                               nhi_field_name='nhi',
                               diagnosis_field_name = 'Tim_diag2')   

        process_diagnosis_file(name="MSPD Society",
                               filename='input/diagnoses_mspd_2015.csv',
                               nhi_field_name='nhi',
                               diagnosis_field_name = 'mspd_diag2')      
        
        process_diagnosis_file(name="Clinics",
                               filename='input/diagnoses_clinic_2015.csv',
                               nhi_field_name='nhi',
                               diagnosis_field_name = 'diag2')                

        process_diagnosis_file(name="CDHB",
                               filename='input/diagnoses_cdhb_2014.csv',
                               nhi_field_name='nhi',
                               diagnosis_field_name = 'dhb_diag')   

        process_diagnosis_file(name="Neurology database",
                               filename='input/diagnoses_neurology_2015.csv',
                               nhi_field_name='nhi',
                               diagnosis_field_name = 'diag1',
                               diagnosis_detail = 'diag2')          
        

        try:
            inputFile = open(local_diagnoses_filename)
            reader = csv.DictReader(inputFile, delimiter=",")
            
            for row in reader:
                try:
                    current_diagnosis = self.local_diagnoses[row['nhi']]
                except:
                    pass
                    #print row
                
            
#                 for field in ('research_diag','Tim_pp_diag','mspd_diag',
#                               'neurology_diag','clinic_diag','Cant_diag'):
#                     
#                     diag = row[field].replace(" ","")
#                     if diag not in ('NA','unknown'):
#                         if diag == "PD":
#                             self.local_diagnoses[row['nhi']]='PD'
#                             self.all_diagnoses[row['nhi']].append('PD')
#                         else:
#                             self.local_diagnoses[row['nhi']]='Other'
#                             self.all_diagnoses[row['nhi']].append('Other')
#                         break
        except:
            traceback.print_exc(file=sys.stdout)
            print("Unable to open {0}, all diagnoses will be empty".format(local_diagnoses_filename))
        
        
        
        
        in_moh = 0
        diff_diag = 0
        for nhi in rd_diags:
            try:
                moh_diag = self.moh_diagnoses[nhi]
                in_moh +=1
                if moh_diag != 'PD':
                    diff_diag +=1
            except KeyError:
                pass
            
        print "{} diags from research database, {} in MOH, {} different".format(len(rd_diags.keys()),
                                                                                in_moh,diff_diag)       
        
        multiple_diagnoses = 0
        diff_diags = 0
        for nhi in self.all_diagnoses:
            diags = self.all_diagnoses[nhi]
            if len(diags) > 1:
                multiple_diagnoses +=1 
                for i in xrange(0,len(diags)-1):
                    if diags[i] != diags[i+1]:
                        diff_diags += 1
                        break
                 
        print "Number of diagnoses {} multiple {} different {}".format(len(self.all_diagnoses.keys()), multiple_diagnoses, diff_diags)
    
    def getLocalDiagnosis(self,nhi):
        
        # Get local diagnosis if exists
        try:
            return self.local_diagnoses[nhi]
        except KeyError:
            return 'NA'  
    
    def getMohDiagnosis(self,nhi):
        
        # Get MOH diagnosis if exists
        try:
            return self.moh_diagnoses[nhi]
        except KeyError:
            return 'NA'
        
    def getDiagnosis(self, nhi):
        
        local_diagnosis = self.getLocalDiagnosis(nhi)
        
        moh_diagnosis =  self.getMohDiagnosis(nhi)
        
        # Print if different
        if local_diagnosis != 'NA' and moh_diagnosis != 'NA':
            if local_diagnosis != moh_diagnosis:
                print "{}: Local/MOH diagnosis: {}/{}".format(nhi, local_diagnosis, moh_diagnosis)
        
        #print local_diagnosis, moh_diagnosis
        
        if local_diagnosis != 'NA':
            return local_diagnosis
        else:
            return moh_diagnosis
        
if __name__ == '__main__':
    diagnoses = Diagnoses()
            