from collections import defaultdict, OrderedDict
import dateutil.parser
import datetime
import operator
import sys,traceback
import csv
import sqlite3

def dict_from_row(row):
    return dict(zip(row.keys(), row))

class PharmacData:
    
    def __init__(self, datasets = None, outfname = None, exclude_under_20 = True):
        

        
        self.datasets = datasets
        self.outfname = outfname
        self.exclude_under_20 = exclude_under_20
        
        
        self.dbconn = sqlite3.connect('output/pharmac.db')
        self.dbconn.row_factory = sqlite3.Row
        
        self.db = self.dbconn.cursor()

        # Create table
        self.db.execute('DROP TABLE IF EXISTS dispensings')
        
        self.db.execute('''CREATE TABLE dispensings
                     (nhi text, birthdate text, date_of_death text, age real, sex text,
                      ethnicity text, dhb text, date text, drug text, drug_group text, dose_mg text, 
                      days_supply text)''')
        self.db.execute('''CREATE INDEX Idx1 ON dispensings(nhi)''')
        
        # Map from ethnic ID to ethnicity
        self.ethnic_mapping = {
                'European':('10', '11', '12', '54', '61'), # include other in European, primarily New Zealander
                'Maori':('21',),
                'Pacific':('30', '31', '32', '33', '34', '35', '36', '37'),
                'Asian':('40', '41', '42', '43', '44'), 
                'Other':('51', '52', '53'),
                'Unknown':('94', '95', '97', '99', 'un')
            }
        
        #Map from DHB ID to DHB Name
        self.dhb_mapping = {
                'Northland':('011', '11'),
                'Waitemata':('021', '21'),
                'Auckland':('022', '22'),
                'Counties':('023', '23'),
                'Waikato':('031', '31'),
                'Lakes':('042', '42'),
                'BoP':('047', '47'), 
                'Tairawhiti':('051', '51'),
                'HawkesBay':('061', '61'), 
                'Taranaki':('071', '71'),
                'MidCentral':('081', '81'),
                'Whanganui':('082', '82'),
                'CapitalCoast':('091', '91'),
                'Hutt':('092', '92'),
                'Wairarapa':('093', '93'),
                'NelsonMarlb':('101',),
                'WCoast':('111',),
                'Canterbury':('121',),
                'SCanterbury':('123',),
                'Southern':('131', '160', '141'),
                'Unknown':('UNK',),
            }
        
        self.excluded_drugs = ('Clozapine','Donepezil hydrochloride','Quetiapine')
        
        # IDs to drugs
        self.drugid_mapping = {
                'Biperiden': ('57133',),
                'Kinson':('81042',),
                'Orphenadrine': ('57132', '61698', '78800'),
                'Procyclidine': ('57131', '62207'),
                'Sindopa': ('60314', '60315', '60316'),
                'Entacapone': ('57107', '73250', '76220', '79398','79548'), 
                'Amantadine': ('57128', '72011', '72012', '72013', '72014', '78396'),
                'Tolcapone': ('57126', '67200', '67201', '67474', '67475', '78386'), 
                'Lisuride': ('57111', '69703', '69704', '69705', '69706', '69707', '69708','79808'),
                #'Rivastigmine':('57263', '57264', '81297', '81298', '81325', '81326'),
                #'Donepezil':('77775','77776','81399','81400'),
                'Apomorphine': ('57101', '63560', '75198', '75228', '76160', '76166', '76179', 
                                '77117', '78664'), 
                'Pergolide': ('57129', '57130', '64541', '64542', '64543', '64544', '75473', 
                              '75474'), 
                'Benztropine': ('57134', '57135', '58610', '58611', '58612', '58613', '66401', 
                                '66402', '73154', '73311', '76349','78757'), 
                'Sinemet': ('59005', '59006', '57112', '57113', '57114', '62226', '62227', 
                            '62228', '60314', '60315', '60316', '69668', '69669', '76723', 
                            '76794', '79731', '79732', '79733'), 
                'Pramipexole':('57122', '81236', '81237', '81238', '81239', '81240', '78817', 
                               '78818', '78819', '78848', '78849', '78850', '79877', '79901', 
                               '80501', '80502'),
                'Selegiline': ('57109', '57110', '60371', '60372', '60373', '60374', '60375', 
                               '60376', '60377', '60378', '66419', '66420', '66421', '66422', 
                               '66423', '66424', '69790', '69791', '69792', '69793', '69795', 
                               '69796', '73856', '69794', '71562', '71563', '71564', '77785',
                               '78500', '80038', ''),
                'Madopar': ('58563', '58564', '58565', '58566', '58567', '58568', '58569', 
                            '58570', '58571', '58572', '58573', '58574', '58575', '58576', 
                            '57115', '57116', '57117', '57118', '57119', '57120', '57121', 
                            '69566', '69567', '69568', '69569', '69570', '69571', '68947', 
                            '68948', '68949', '68950', '68951', '58562'), 
                'Bromocriptine': ('58553', '58554', '58555', '58556', '58557', '58558', '58559', 
                                  '58560', '57123', '57124', '57125', '61120', '61121', '61122', 
                                  '61124', '61125', '61126', '61127', '61128', '61129', '61130', 
                                  '61131', '65590', '65591', '65592', '65593', '65594', '65595', 
                                  '70398', '70399', '70401', '70402', '70403', '70404', '70405', 
                                  '70406', '58561', '72163', '72164', '72165', '72162', '70400', 
                                  '61123', '76633', '76668', '76682', '76894', '76905'), 
                'Ropinirole': ('73276', '73277', '73278', '73279', '76162', '76163', '76164', 
                               '76165', '76279', '76280', '76281', '76282', '77344', '77345', 
                               '77346', '77347', '74417', '76287', '74418', '76288', '80481',
                               '80482', '80483', '80484', '80642', '80644', '80646', '80648'),
            }
        
        # Map from Drug to Group
        self.drug_mapping = {
                'L-dopa':('Sinemet', 'Sindopa', 'Madopar','Kinson'),
                'COMT':('Entacapone', 'Tolcapone'),
                'DA agonist':('Lisuride', 'Pergolide', 'Ropinirole', 'Bromocriptine', 
                              'Apomorphine', 'Pramipexole'),
                'Anticholinergic':('Orphenadrine', 'Benztropine', 'Procyclidine'),
                'MAOI':('Selegiline',),
                'Amantadine':('Amantadine',),
                'Dementia':('Rivastigmine','Donepezil'),
                'Gout':('Allopurinol','Colchicine'),
                'CCB':('Nifedipine','Felodipine','Isradipine','Amlodipine'),
                'Metformin':('Metformin hydrochloride'),
            }
        
        self.dose_mapping = {
                
                "0.125":('78849',),
                "0.20":('57111', '69703', '69704', '69705', '69706', '69707', '69708'),
                "0.25":('57130', '64543', '64544', '75474', '57122', '57106', '73276', '76143', '76162',
                        '76279', '77344', '78848', '80501', '81236'),
                "0.50":('78850',),
                "1.0":('57129', '64541', '64542', '75473', '57105', '73277', '76144', '76163', '76280', '77345', 
                       '73154', '79901','80502'),
                "2.0":('58612', '58613', '57135', '66402', '73311', '66401', '76349', '57104', '73278', '76145', 
                      '76164', '76281', '77346', '58610', '58611', '57134', '57109', '71562', '71563', '71564'),
                "2.5":('58557', '58558', '58559', '58560', '57125', '61120', '61121', '61122', '61124', '61125',
                      '65590', '65591', '65592', '65593', '65594', '65595', '70398', '70399', '70401', '70402',
                      '58561', '72163', '72164', '72165', '72162', '70400', '61123', '76894', '76905'),
                "4.6":('81326',),
                "5.0":('76633', '76668', '57131', '62207', '57110', '60371', '60372', '60373', '60374', '60375', 
                       '60376', '60377', '60378', '66419', '66420', '66421', '66422', '66423', '66424', '69790',
                       '69791', '69792', '69793', '69795', '69796', '73856', '69794', '77785', '57108', '73279',
                       '57103', '76146', '76165', '76282', '77347', '57131'),
                "5.5":('57133',),
                "9.5":('81325',),
                "10.0":('58553', '58554', '58555', '58556', '57123', '57124', '61126', '61127', '61128', '61129', '61130', 
                        '61131', '70403', '70404', '70405', '70406', '76682', '57101', '63560', '78664', '63560', '78757'),
                "20.0":('75228', '76166', '77117'),
                "50.0":('58563', '58564', '58565', '58566', '57117', '58562', '57132', '61698', '57116', '68947', '68948',
                        '68949', '68950', '68951','78800'),
                "100.0":('57128', '72011', '72012', '72013', '72014', '78396', '58567', '58568', '58569', '58570', '58571', 
                         '58572', '57119', '57118', '69566', '69567', '69568', '69569', '69570', '69571', '57126', '67200', 
                         '67474', '67475', '67201', '78386', '57114', '62226', '62227', '62228', '60314', '60315', '60316', '78396'),
                "200.0":('58573', '58574', '58575', '58576', '57120', '57107', '73250', '76220', '57121', '57112', '69668',
                         '69669', '76794'),
                "250.0":('57115', '59005', '59006', '57113', '76723'),
                 ##these last lot are the ropinirole packs  
                 "0.0":('74348', '74417', '76157', '76287', '74349', '74418', '76158', '76288')
                 }                                  

        
    def map_item(self,item,mapping):
    
        for group in mapping:
            if str(item) in mapping[group]:
                return group
        
        return 'ATTN'
        #if 'ATTN' in item:
        #    return '{}'.format(item)
        #else:
        #    return 'ATTN-{}'.format(item)
        
    def process_raw(self):
        
        doderrors_file = 'output/disepensing_after_dod.csv'
        singledisp_file = 'output/single_dispensing.csv'
        
        # Processed records file
        fields =[('nhi',1),
                 ('birthdate',2),
                 ('date_of_death',2),
                 ('age',2),
                 ('sex',2),
                 ('ethnicity',3),
                 ('dhb',3),
                 ('date',4),
                 ('drug',5),
                 ('drug_group',6),
                 ('dose_mg',7),
                 ('days_supply',8)
                 ]
        
        f_out = open(self.outfname,"w")
        dwp = csv.DictWriter(f_out, delimiter=',',restval='NA',fieldnames=OrderedDict(fields))
        dwp.writeheader()
        
        # Processed records file
        fields +=[
                 ('dod_delta',9)
                 ]
        
        fsd_out = open(singledisp_file,"w")
        dwsd = csv.DictWriter(fsd_out, delimiter=',',restval='NA',fieldnames=OrderedDict(fields))
        dwsd.writeheader()
        
        fields =OrderedDict([('nhi',1),
                     ('birthdate',2),
                     ('dod',3),
                     ('date',4),
                     ('days_after_dod',5),
                     ('dispenser_fee',6),
                     ('subsidy_value',7),
                     ('provider_id',8),
                     ('drug',9),
                     ])
        
        dwd = csv.DictWriter(open(doderrors_file,"w"), delimiter=',',restval='NA',fieldnames=fields)
        dwd.writeheader()
     
        
        dispensings = defaultdict(list)
        
        n_excluded_records_nhi = 0
        n_excluded_records_age = 0
        n_excluded_records_dod = 0
        n_excluded_records_drug = 0
        
        n_records = 0
        people = set()
        
        excluded_age = set()
        excluded_age_dob = set()
        excluded_dod = set()
        excluded_drug = set()
        
        excluded_drug_names = set()
        missing_dose = set()
        
        total_records_by_year = defaultdict(int)
        missing_nhi_by_year = defaultdict(int)
        
        missing_by_drug = defaultdict(int)
        total_by_drug = defaultdict(int)
        
        for dataset in self.datasets:
            print "Processing file {}".format(dataset['filename']) 
            with open("raw/"+dataset['filename'], "r") as f:
            
                fk = open("raw/"+dataset['key'], "r")
                
                drug_names = dict()
                keys = csv.DictReader(fk)
                for key in keys: 
                    drug_names[key['DIM_FORM_PACK_SUBSIDY_KEY']]=key['CHEMICAL_NAME']
                
                records = csv.DictReader(f)
                row = 1
                for record in records:
                    
                    
                    ## Testing - only process first 10000
                    #if row > 10000:
                    #    break
                    #else:
                    #    row += 1
    
                    ## Extract data and handle exclusion cases at the records level
                    
                    drug_id = record['DIM_FORM_PACK_SUBSIDY_KEY']
                    #drug = self.map_item(drug_id,self.drugid_mapping)
                    drug = drug_names[drug_id]
                    drug_group = self.map_item(drug,self.drug_mapping)
                    if drug_group == 'ATTN':
                        drug_group = drug
                    nhi = record[dataset['nhi']]
                    date = record['DATE_DISPENSED']
                    date_py = dateutil.parser.parse(date,dayfirst=True)


                    
                    ## Testing between prim_hcu and nhi
                    #nhi2 = record["prim_hcu"]
                    #if nhi!= nhi2:
                    #    nhi_diff[nhi2].add(nhi)
                    
                    n_records += 1
                    total_records_by_year[date_py.year]+=1
                    total_by_drug[drug_group]+=1
                    
                    # Record NHI if known
                    if nhi not in ('','unknown'):
                        people.add(nhi)
                    
                    # Only include if antiparkinson's
                    if drug in self.excluded_drugs:
                        excluded_drug.add(nhi)
                        n_excluded_records_drug += 1
                        try:
                            excluded_drug_names.add("{}-{}".format(drug,drug_id))
                        except:
                            excluded_drug_names.add("{}".format(drug_id))
                        continue
                 
                    
                    # If NHI is empty exclude
                    if nhi in ('','unknown'):
                        n_excluded_records_nhi +=1
                        missing_nhi_by_year[date_py.year]+=1
                        missing_by_drug[drug_group]+=1
                        #print record
                        continue
                    
                    # Only have age if have NHI
                    dob = record['dob']
                    age = (date_py-dateutil.parser.parse(dob,dayfirst=True)).days/365.0                    
                    
                    # if DOD is before dispensing date obviously an error
                    dod = record[dataset['dod']]
                    if dod != '':
                        dod_py = dateutil.parser.parse(dod,dayfirst=True)
                        if dod_py < date_py:
                            data = {'nhi':nhi,
                                    'birthdate':record['dob'],
                                    'date':date,
                                    'dod':dod,
                                    'days_after_dod':(date_py-dod_py).days,
                                    'dispenser_fee':record['DISPENSING_FEE_VALUE'],
                                    'subsidy_value':record['RETAIL_SUBSIDY'],
                                    'provider_id': record['PROVIDER_NUMBER'],
                                    'drug':drug,
                                    }
                            dwd.writerow(data)
                            excluded_dod.add(nhi)
                            n_excluded_records_dod +=1
                            continue
                            
                    # If younger than 20 years exclude
                    if nhi in excluded_age:
                        n_excluded_records_age +=1
                        continue
                    
                    if self.exclude_under_20 and age < 20:
                        excluded_age.add(nhi)
                        excluded_age_dob.add((nhi,dob))
                        n_excluded_records_age +=1
                        continue
                    
                    ethnicity = self.map_item(record['ETHNICGP'],self.ethnic_mapping)
                    dhb = self.map_item(record['DHB_CLAIMANT'],self.dhb_mapping)
                    
                    

                    
                    
                    try:
                        dose_mg = float(self.map_item(record['DIM_FORM_PACK_SUBSIDY_KEY'],self.dose_mapping))
                        days = record['DAILY_DOSE']
                        if days != '':
                            dose_mg *= float(days)
                            dose_mg = "{:0.2f}".format(dose_mg)
                        else:
                            dose_mg = 'NA'
                    except ValueError:
                        missing_dose.add(record['DIM_FORM_PACK_SUBSIDY_KEY'])
                        dose_mg = 'NA-{}'.format(record['DIM_FORM_PACK_SUBSIDY_KEY'])
                    
                    days_supply = record['DAYS_SUPPLY']
                    if days_supply == '0':
                        days_supply = 'NA'
                    
                    summary = {'nhi':nhi,
                               'age':'{:0.1f}'.format(age),
                               'sex':record['GENDER'],
                               'birthdate':record['dob'],
                               'date_of_death':dod,
                               'date':date,
                               'ethnicity':ethnicity,
                               'dhb':dhb,
                               'drug':drug,
                               'drug_group':drug_group,
                               'dose_mg':dose_mg,
                               'days_supply':days_supply
                               }
                    
                    ## OLD: store in a dictionary
                    #dispensings[nhi].append(summary)
                    
                    # New: Put in a DB:
                    self.db.execute('INSERT INTO dispensings ' + 
                                '(nhi, age, sex, birthdate, date_of_death, date, ethnicity, ' +
                                'dhb, drug, drug_group, dose_mg, days_supply) ' +
                                'VALUES (:nhi, :age, :sex, :birthdate, :date_of_death, :date, :ethnicity, ' +
                                ':dhb, :drug, :drug_group, :dose_mg, :days_supply);', summary)
                    
        
        self.dbconn.commit()
                    
        n_final_people = 0
        n_final_records = 0
        
        n_excluded_records_single = 0
        n_excluded_people_single = 0
        
        excluded_single_months = defaultdict(int)
        
        print "All records read in. Now exporting by individual"
        
        #for person in sorted(dispensings.keys()):
        #    sorted_dispensings = sorted(dispensings[person], key=lambda k: k['age'])
        
        persons = self.db.execute("SELECT DISTINCT nhi FROM dispensings ORDER BY nhi")
        for person in persons.fetchall() :
            sorted_dispensings = self.db.execute("SELECT * FROM dispensings WHERE nhi=? ORDER BY age",person).fetchall()
            
            
            
            # Count number of unique dates
            dates = set()
            for dispensing in sorted_dispensings:
                dates.add(dispensing['date'])
            
            # Export data if dispensings on 2 or more dates
            if len(dates) > 1:
                n_final_people += 1
                for dispensing in sorted_dispensings:
                    n_final_records +=1
                    dwp.writerow(dict_from_row(dispensing))
            else:
                for dispensing in sorted_dispensings:
                    dispensing = dict_from_row(dispensing)
                    if dispensing['date_of_death']:
                        date_py = dateutil.parser.parse(dispensing['date'],dayfirst=True)
                        dod_date_py = dateutil.parser.parse(dispensing['date_of_death'],dayfirst=True)
                        dispensing['dod_delta']=(dod_date_py-date_py).days
                    dwsd.writerow(dispensing)
                n_excluded_records_single += len(sorted_dispensings)
                n_excluded_people_single += 1
                date_py = dateutil.parser.parse(dates.pop(),dayfirst=True)
                excluded_single_months[date_py.strftime("%Y-%m")]+=1
                
        n_people = len(people)
        
        
        empty_nhi=set(('','unknown'))
        n_excluded_people_drug = len(excluded_drug - set(dispensings) - empty_nhi)
        n_excluded_people_dod = len(excluded_dod - excluded_drug - set(dispensings) - empty_nhi)
        n_excluded_people_age = len(excluded_age - excluded_drug - excluded_dod - set(dispensings) - empty_nhi)
        
        
        print "{} records from {} people in raw data".format(n_records,n_people)
        
        print "{} records excluded and {} people removed due to only antipsychotic/dementia drug".format(n_excluded_records_drug,
                                                                                                         n_excluded_people_drug)
        
        n_records_remain = n_records-n_excluded_records_drug
        n_people_remain = n_people-n_excluded_people_drug
        
        print "{} records and {} people remain".format(n_records_remain,
                                                       n_people_remain)
        
        print "{} records excluded due to missing NHI".format(n_excluded_records_nhi)
        
        print "{} records excluded and {} people removed due to date of death before dispensing date".format(n_excluded_records_dod,
                                                                                                             n_excluded_people_dod)
        
        print "{} records excluded and {} people removed due to age < 20".format(n_excluded_records_age,
                                                                                 n_excluded_people_age)
        
        print "{} records excluded and {} people removed due to only having a single date of dispensing".format(n_excluded_records_single,
                                                                                                                n_excluded_people_single)
        
        n_records_remain_exlc = n_records_remain - n_excluded_records_nhi - n_excluded_records_dod - n_excluded_records_age - n_excluded_records_single
        n_people_remain_excl = n_people_remain - n_excluded_people_dod - n_excluded_people_age - n_excluded_people_single
        
        print "{} records and {} people in final dataset (based upon exclusion counts)".format(n_records_remain_exlc,
                                                                                         n_people_remain_excl)
        
        print "{} records and {} people in final dataset (based upon actual records exported)".format(n_final_records,n_final_people)
        
        
        
        print "Drugs excluded from final dataset:"
        for drug in sorted(excluded_drug_names):
            print drug
        
        print "Missing doses for these drugs:"
        print missing_dose
        
        print "Years and months of people with only a single prescription:"
        print excluded_single_months
        
        for year in sorted(total_records_by_year.keys()):
            print "{} - missing {:.1f}%".format(year,missing_nhi_by_year[year]*100.0/total_records_by_year[year])
            
        for drug in sorted(total_by_drug.keys()):
            print "{} - missing {:.1f}%".format(drug,missing_by_drug[drug]*100.0/total_by_drug[drug])
        
        #f_age_out = open("age_excluded.txt","w")
        #for item in excluded_age_dob:
        #    f_age_out.write("{}\n".format(item))
        

        
if __name__ == '__main__':

    pd_datasets = (
            {'filename':'phh0256/part1.csv',
             'key':'phh0256/dim_form_pack_subsidy.csv',
             'nhi':'PRIM_HCU',
             'dod':'nhi_dod'
             },
            {'filename':'phh0436/part1.csv',
             'key':'phh0436/dim_form_pack_subsidy.csv',
             'nhi':'prim_hcu',
             'dod':'nhi_dod'
             },
            {'filename':'phh0445/part1.csv',
             'key':'phh0436/dim_form_pack_subsidy.csv',
             'nhi':'prim_hcu',
             'dod':'nhi_dod'
             },
            
            )

    #pharmac = PharmacData(pd_datasets,'output/included_records.csv')
    #pharmac.process_raw()
    
    new_datasets = (
                    {'filename':'phh0563/PHH0563_2005.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2006.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2007.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2008.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2009.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2010.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2011.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2012.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2013.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    {'filename':'phh0563/PHH0563_2014.csv',
                     'key':'phh0563/DIM_FORM_PACK_SUBSIDY.csv',
                     'nhi':'PRIM_HCU',
                     'dod':'DOD'
                     },
                    )
                    
    pharmac = PharmacData(new_datasets,
                          'output/included_records_pd_protection.csv',
                          exclude_under_20 = False
                          )
    pharmac.process_raw()    