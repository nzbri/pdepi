
library(ggplot2)
library(dplyr)
library(tidyr)
library(xlsx)
library(lme4)
library(reshape2)
library(rstan)
library(rstanarm)
library(brms)
library(shinystan)
library(forcats)
library(xlsx)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores ())
seed = 123

## Classification output (generated from Python)
Pall <- read.csv("output/classification_v2.csv")


##
## Age cuts for plotting and standardisation
##
#
# 5 year cuts upto 100:
# Used for nz-wide age plots for 2013 (where know population upto 100 from census)
#
# 5 year cuts upto 95
# Used for future predictions where provided to 95+
#
# 5 year cuts upto 90
# Used for nz-wide age by year plots and ethnicity (where only get upto 90+)
#
# 5 year cuts upto 85
# Used for standardised prevalence and incidence

Pall %<>%
  mutate(ethnicity = relevel(ethnicity,"European")) %>% # Ethnic group with largest population
  mutate(classification = relevel(classification,"Very probable")) %>% # Most likely to have PD
  mutate(agecut5_100 = cut(age,labels=c(20+0:15*5,'100+'), breaks = c(0,25+0:15*5,200))) %>%
  mutate(agecut5_95 = cut(age,labels=c(20+0:14*5,'95+'), breaks = c(0,25+0:14*5,200))) %>%
  mutate(agecut5_90 = cut(age,labels=c(20+0:13*5,'90+'), breaks = c(0,25+0:13*5,200))) %>%
  mutate(agecut5_85 = cut(age,labels=c(20+0:12*5,'85+'), breaks = c(0,25+0:12*5,200))) %>%
  mutate(agecut10_80 = cut(Pall$age,labels=c(20+0:5*10,'80+'), breaks = c(0,30+0:5*10,200)))

##############################
##
## Subsets of data
##

## Dataset of everyone from 2006 to 2103, remove two cases of unknown sex
## Ensures everyone has at least a one year window where data could be collected
P <- subset(Pall,year>2005&year<=2013&sex!="U")
P$sex <- factor(P$sex)

# People in 2013 only
P2013 <- subset(P,year==2013)

# first time in dataset all years (calibration and incidence)
Pfirstall <- subset(Pall,year_in_data==1)

# first time in dataset from 2006 to 2013 (calibration and incidence)
Pfirst <- subset(P,year_in_data==1)

# first time in dataset for 2013 (calibration and incidence)
Pfirst2013 <- subset(Pfirst,year==2013)

####################################
##
## Fixed variables
##

# Correction for percent of population that isn't medicated. 
unmedicated_correction = 1.05

# Factor to convert to per 100,000 rates
to_per_100000 = 100000

# Correction for Level-II multi-count census data
# MoH data is prioritised ethnicity
# Total count over ethnicities divided by total population
# Correct for Census ethnic population data allowing multiple ethnicities
# NOT USED - Corrected by age bands instead
multiple_ethnicities_stated = 0.90

# Correct for unknown ethnicity in MoH data. Assume MAR by ethnicity
unknown_ethnicity = 1.025

# 2013 specific values for missing data corrections

missing_nhi_correction_2013 = 1/0.998 # Prevalence & Incidence - based upon observed NHI missing rate

missing_nhi <- data.frame(year=seq(2006,2013),
                          percent=c(96.5,97.9,98.5,99,99.4,99.6,99.7,99.8)) %>%
  mutate(correction = 100/percent)


# Incidence (some new cases that appear will be due to reduction in missing NHI)
reduced_missing_correction_2013 = 0.1/100  # Based upon observed NHI missing data

reduced_missing_nhi <- data.frame(year=seq(2006,2013),
                                  percent=c(3.7,1.4,0.6,0.5,0.5,0.1,0.1,0.1)) %>%
  mutate(fraction = percent/100)

# Used when standardising by ethnicity for DHB
ethnic_percent = data.frame(ethnicity=c("European","Asian","Pacific","Maori"),
                            percent=c(0.70,0.10,0.06,0.13))

average_over_n_years = 8

##
## Census and standardised data
##

# ERP population by year
popyn <- read.csv("input/pop-by-year-nhi-percent.csv")

# ERP population by year and sex
popys <- read.csv("input/pop-by-year-sex.csv")

# ERP population by year, sex, age
popysa <- read.csv("input/pop-by-year-sex-age5.csv")

# 2013 ERP population upto 100+ (using census proportions for 90, 95, 100+)
# Superseded by popysa with age split extrapolations
# popas2013erp <- read.csv("input/pop-by-sex-age5-2013-erp.csv")

# Census population upto 100+
# Superseded by popysa with age split extrapolations
# popas2013census <- read.csv("input/pop-by-sex-age5-2013-census.csv")

#popa2013census <- popas2013census %>%
#  filter(sex == 'All') %>%
#  select(-sex)

# ERP population by year, age
popya <- popysa %>%
  group_by(year,agecut5_90) %>%
  summarise(pop=sum(pop))

# DHB ERP by age and sex 2013
popdhbas <- read.csv("input/dhb-age-sex-2013.csv")

# DHB ERP
popdhb <- popdhbas %>%
  filter(sex=='All') %>%
  group_by(dhb) %>%
  summarise(pop=sum(pop))

# DHB Census by age, sex, ethnicity
popdhbase_raw <- read.csv("input/dhb-ethnicity-age-sex-2013-processed.csv")
popdhbase <- popdhbase_raw %>%
  filter(!sex=='Total') %>%
  filter(!age=='Total') %>%
  gather(key=ethnicity,value=pop,
         Maori,Pacific,Asian,MELAA,Other,European,
         Total.Stated,Not.Elsewhere.Included,Total) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(pop = as.numeric(sub(",","",pop))) %>%
  mutate(agecut5_85 = cut(age,labels=c(0:16*5,'85+'), breaks = c(0,5+0:16*5,200))) %>%
  group_by(sex,agecut5_85,dhb,ethnicity) %>%
  summarise(pop = sum(pop))


#Population by ethnicity, 2013 (ERP, multiple counts)
pope2013 <- read.csv("input/pop-by-ethnicity-2013.csv")

#Population by ethnicity, 5 year age splits, 2013 (ERP, multiple counts)
popea2013 <- read.csv("input/pop-by-ethnicity-age5-2013.csv")

#Population by ethnicity, 5 year age splits, sex split, 2013 (ERP, multiple counts)
popeas2013 <- read.csv("input/pop-by-ethnicity-age5-sex-2013.csv")

popes2013 <- popeas2013 %>%
  group_by(ethnicity,sex) %>%
  summarise(ethnicpoperp=sum(ethnicpoperp))

#Population by ethnicity, 5 year age splits, sex split, 2013 (ERP, multiple counts)
# Not used
popeas2013cr <- read.csv("input/pop-by-ethnicity-age5-sex-2013-multicount-corrected.csv")

popeas2013c <- popeas2013cr %>%
  gather(ethnicity, poperp, -sex, -agecut) %>%
  filter(!agecut %in% c('0','5','10','15','85',"90+")) %>%
  mutate(ethnicity = factor(ethnicity)) %>%
  filter(!ethnicity %in% c('Total','TotalMC','AsianMC','EuropeanMC',
                           'FactorMC','MaoriMC','OtherMC','PacificMC')) %>%
  mutate(ethnicity = relevel(ethnicity,"European")) %>%
  mutate(ethnicity = factor(ethnicity)) %>%
  mutate(agecut = factor(agecut))

#Population by ethnicity, 5 year age splits, aex splits, 2006 and 2013, ERP, multicount
popeasr <- read.csv("input/pop-by-ethnicity-age5-sex-years.csv")

  popeasc <- popeasr %>%
    ## Generate new 85+ agecut
    spread(agecut,erppop) %>%
    mutate(`85+`=`85`+`90+`) %>%
    tidyr::gather("agecut","erppop",4:23) %>%
    ## Account for multicount
    spread(ethnicity,erppop) %>%
    mutate(totalmc = Asian+European+Maori+Other+Pasifika) %>%
    mutate(mcfrac = All/totalmc) %>%
    mutate(European = European*mcfrac,
           Asian = Asian*mcfrac,
           Maori = Maori*mcfrac,
           Other = Other*mcfrac,
           Pasifika = Pasifika*mcfrac) %>%
    select(-totalmc,-mcfrac) %>%
    tidyr::gather("ethnicity","poperp",4:9) %>%
    filter(!ethnicity == 'All') %>%
    mutate(poperp = as.integer(poperp)) %>%
    ## Generate person-years 
    ## Assume linear population increase from 2006 to 2013
    spread(Year,poperp) %>%
    mutate(personyears=4*`2006`+4*`2013`) %>%
    mutate(`2007`=as.integer(`2006`+1/7*(`2013`-`2006`)),
           `2008`=as.integer(`2006`+2/7*(`2013`-`2006`)),
           `2009`=as.integer(`2006`+3/7*(`2013`-`2006`)),
           `2010`=as.integer(`2006`+4/7*(`2013`-`2006`)),
           `2011`=as.integer(`2006`+5/7*(`2013`-`2006`)),
           `2012`=as.integer(`2006`+6/7*(`2013`-`2006`))
    ) %>%
    tidyr::gather("year","poperp",4:12) %>%
    ## Cleanup
    mutate(ethnicity = factor(ethnicity)) %>%
    mutate(ethnicity = relevel(ethnicity,"European")) %>%
    mutate(agecut = factor(agecut)) %>%
    mutate(year = factor(year))

popec <- popeasc %>%
  filter(!agecut %in% c('85',"90+")) %>%
  group_by(ethnicity,year) %>%
  summarise(poperp=sum(poperp))

#Standard population age (based on 2013 NZ ERP by 5-year age groups)
popstanda <- read.csv("input/standard-pop-age.csv")

#Standard population age-sex (based on 2013 NZ ERP by 5-year age groups)
popstandas <- read.csv("input/standard-pop-age-sex.csv")

## Data used for calibration of model

# Repeated diagnoses for each year
# Not used
calcyse <- P %>%
  filter(!is.na(diagnosis)) %>%
  mutate(pd = as.integer(diagnosis)-1) %>%
  mutate(ethnicity = relevel(ethnicity,"European")) %>%
  mutate(classification = fct_relevel(classification,"Very probable","Probable","Possible","Unlikely")) %>%
  mutate(year = year - 2006)

# Limit to 2013
calcyse2013 <- calcyse %>%
  subset(year == 7)


# Only take each individual once
calcysef <- Pall %>%
  filter(!is.na(diagnosis)) %>%
  filter(ethnicity != "Other") %>%
  droplevels() %>%
  mutate(pd = as.integer(diagnosis)-1) %>%
  mutate(ethnicity = relevel(ethnicity,"European")) %>%
  mutate(classification = fct_relevel(classification,"Very probable","Probable","Possible","Unlikely")) %>%
  subset(year_in_data == 1) %>%
  subset(sex!="U") %>%
  mutate(sex = factor(sex))

# Limit to 2013
calcysc2013 <- calcyse2013 %>%
  count(classification,diagnosis,ethnicity,agecut5_85,sex) %>%
  ungroup %>%
  complete(classification,diagnosis,ethnicity,agecut5_85,sex, fill = list(n=0)) %>%
  spread(key=diagnosis,value=n) %>%
  mutate(total = Other+PD)
