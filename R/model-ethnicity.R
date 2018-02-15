

############################################
##
## Poststratification via Stan -- Ethnicity
##
## Similar to multilevel regression and poststratification
## http://andrewgelman.com/2014/01/15/postdoc-involving-pathbreaking-work-mrp-stan-2014-election/
## http://www.princeton.edu/~jkastell/mrp_primer.html
## https://github.com/malecki/mrp
##
## Examples:
## https://github.com/stan-dev/example-models/wiki/ARM-Models
## https://github.com/stan-dev/example-models/blob/master/ARM/Ch.14/election88_full.stan
##

## Extract data for methods illustration
pd_methods_illustration <- Pfirstall %>%
  select(classification,ethnicity,agecut5_85,sex)
write.csv(pd_methods_illustration,"output/pd_epi_demographics.csv")

ages_85 <- c("20","25","30","35","40","45",
             "50","55","60","65","70","75","80","85+")

# Prevalence Counts by classification - long format
pdcounts2013lp <- P2013 %>%
  filter(ethnicity!="Other") %>%
  droplevels() %>%
  mutate(classification = fct_relevel(classification,"Very probable","Probable","Possible","Unlikely")) %>%
  count(classification,ethnicity,sex,agecut5_85) %>%
  ungroup %>%
  complete(classification,ethnicity,sex,agecut5_85, fill = list(n=0)) %>%
  mutate(type = "Prevalence")

# Incidence Counts over all years by classification - long format
pdcounts2013li <- Pfirst %>%
  filter(ethnicity!="Other") %>%
  droplevels() %>%
  mutate(classification = fct_relevel(classification,"Very probable","Probable","Possible","Unlikely")) %>%
  count(classification,ethnicity,sex,agecut5_85) %>%
  ungroup %>%
  complete(classification,ethnicity,sex,agecut5_85, fill = list(n=0)) %>%
  mutate(type = "Incidence")

# Combine incidence and prevalence counts
pdcounts2013l <- pdcounts2013lp %>%
  rbind(pdcounts2013li)

# Convert to array for stan
pdcounts2013w <- acast(pdcounts2013l, 
                       type ~ classification ~ ethnicity ~ sex ~ agecut5_85, 
                       value.var="n", fill=0)

# NZ Standard pop by age and sexs
pop_nzsl <- popstandas %>%
  filter(!agecut %in% c('0','5','10','15','85',"90","95","90+","100+"))

# Convert to array
pop_nzsw <- acast(pop_nzsl, 
                  sex ~ agecut, 
                  value.var="nz_standard_pop", fill=0)

# NZ Standard pop by age
pop_nzsal <- popstanda %>%
  filter(!agecut %in% c('0','5','10','15','90+','85','90','95','100+')) %>%
  mutate(agecut = fct_relevel(agecut,ages_85))

# Convert to array
pop_nzsaw <- pop_nzsal$nz_standard_pop

# WHO Standard Pop
pop_whol <- popstanda %>%
  filter(!agecut %in% c('0','5','10','15','85',"90","95","90+","100+"))

# population - variable year
# factors personyears (incidence) 2013 (prevalence)
pop_esal <- popeasc %>%
  subset(year %in% c('personyears',2013)) %>%
  filter(ethnicity != "Other") %>%
  droplevels() %>%
  filter(!agecut %in% c('0','5','10','15','85',"90+")) %>%
  mutate(year = relevel(year,"personyears"))

pop_esaw <- acast(pop_esal,
                  year ~ ethnicity ~ sex ~ agecut,
                  value.var="poperp")

# population - variable year
# factors personyears (incidence) 2013 (prevalence)
pop_el <- popec %>%
  filter(ethnicity != "Other") %>%
  droplevels() %>%
  subset(year %in% c('personyears',2013)) %>%
  mutate(year = relevel(year,"personyears"))

pop_ew <- acast(pop_el,
                year ~ ethnicity,
                value.var="poperp")

# Classification model matrix
classmm <- model.matrix( ~ 1 + classification + sex, calcysef)

strat_dat <- list(N = nrow(calcysef), # Number of data points
                  n_beta = 5, # number of predictors
                  n_pdclass = nlevels(calcysef$classification),
                  n_ethnicity = nlevels(calcysef$ethnicity), # Number of ethnicities
                  n_sex = 2, # Male/Female
                  n_types = 2, # Prevalence/Incidence
                  n_age = nlevels(calcysef$agecut5_85), # Number of age groups
                  ethnicity = as.numeric(calcysef$ethnicity), # Ethnicity for each data point
                  age = as.numeric(calcysef$agecut5_85), # Age group for each data point
                  pdclass = as.numeric(calcysef$classification), # Classification for each data point
                  X = classmm, # Calibration model matrix
                  y = calcysef$pd, # Known diagnosis
                  counts_cesa = pdcounts2013w, # Counts by type, classification, etthnicity, sex, age for binomial model
                  pop_e = pop_ew, # Population by ethnicity
                  pop_esa = pop_esaw, # Population by ethnicity, sex, age
                  pop_nzs = pop_nzsw, # Standard population
                  pop_nzsa = pop_nzsaw, # Standard population
                  pop_whos = pop_whol$who/100)


pars <- c('a','e','beta','sigma_a','sigma_e','total_count',
          'ethnic_observed_rates','ethnic_observed_rates_normapprox',
          'ethnic_standardised_rates','ethnic_standardised_rates_who',
          'ethnic_rates_age','ethnic_count','ethnic_corrected_count','ethnic_sex_age_count',
          'mean_standardised_duration','mean_observed_duration',
          'ethnic_standardised_ratios_sex','ethnic_standardised_rates_sex',
          'theta','theta_mean','sigma_theta','theta_as','theta_as_mean','sigma_theta_as',
          'ethnic_rates_sex_age_normapprox','ethnic_rates_sex_age')

# Stan model
# Only linear transformations of parameters are being made hence warnings can be safely ignored 
mod <- stan_model('pd_epi_model_ethnicity_v2.stan')

# Stan VB - approximation (NOT USED, only for initial testing during model development)
# Will NOT work with current iteration of model
#
#fit_vb <- vb(object = mod, data = strat_dat, seed=seed,
#             pars = pars,algorithm="meanfield")
#
#print(fit_vb, digits_summary=3, probs=c(0.025,0.5,0.975),
#      pars = c('ethnic_standardised_rates'))
#
#print(fit_vb, digits_summary=6, probs=c(0.025,0.5,0.975),
#      pars = c('theta'))
#
#print(fit_m2, digits_summary=1, probs=c(0.025,0.5,0.975),
#      pars = c('ethnic_observed_rates','ethnic_observed_rates_normapprox'))
#
#print(fit_m2, digits_summary=1, probs=c(0.025,0.5,0.975),
#      pars = c('theta_mean','sigma_theta'))
#
#print(fit_m2, digits_summary=1, probs=c(0.025,0.5,0.975),
#      pars = c('ethnic_rates_age'))
#
#my_sso <- launch_shinystan(fit_vb)

SEED=123

## Full inference - full run with 4 chains and 500 iterations.
## Takes 2-3 hours on modern machine to run

fit_m2 <- sampling(mod, data = strat_dat, seed = SEED,
                pars = pars, chains = 4, iter = 500,
                control = list(max_treedepth = 15,adapt_delta=0.99))

fit <- fit_m2

# Observed 'rates': direct calcs and proportion via normal approximation to binomial distribution
# Proportions in the prevalence case
print(fit, digits_summary=1, probs=c(0.025,0.5,0.975),
      pars = c('ethnic_observed_rates','ethnic_observed_rates_normapprox'))

# Standardised 'rates' (proportions in the prevalence case)
print(fit, digits_summary=1, probs=c(0.025,0.5,0.975),
      pars=c('ethnic_standardised_rates',
             'ethnic_standardised_rates_who','ethnic_standardised_ratios_sex'))

# Counts
print(fit, digits_summary=1, probs=c(0.025,0.5,0.975),
      pars=c('ethnic_corrected_count'))

# Model parameter
print(fit, digits_summary=4, probs=c(0.025,0.5,0.975),
      pars=c('theta_as'))

## GGplot of prevalence/incidence by age for each ethnicity with intervals

sim_summary <- summary(fit,
                       pars=c("ethnic_rates_age"),
                       probs=c("0.25","0.75"))$summary
estimated_values <- sim_summary[,c("mean", "25%", "75%")]

# Assesmble a data frame to pass to ggplot()
sim_df <- data.frame(type = c(rep("Incidence",14*4),rep("Prevalence",14*4)),
                     agecut5 = rep(ages_85,2*4),
                     Ethnicity = rep(c(rep("European",14),rep("Asian",14),rep("Maori",14),
                                       rep("Pasifika",14)),2))  %>%
  mutate(middle = estimated_values[, "mean"]) %>%
  mutate(lower = estimated_values[, "25%"]) %>%
  mutate(upper = estimated_values[, "75%"]) %>%
  mutate(agecut5 = fct_relevel(agecut5,ages_85)) %>%
  filter(!Ethnicity=='Other') %>%
  mutate(Ethnicity = fct_relevel(Ethnicity,c('European','Asian','Pasifika','Maori'))) %>%
  droplevels()

levels(sim_df$Ethnicity) <- c('European','Asian','Pasifika','Māori')

write.csv(sim_df,"output/pd_ethnicity_age_specific_rates.csv")

### Age-specfific prevalence plot by ethnicity
ggplot(subset(sim_df,type=='Prevalence'),
       aes(x=agecut5,y=middle,group=Ethnicity,colour=Ethnicity))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Ethnicity,colour=NULL),alpha=0.3)+
  geom_point()+geom_line()+ylab("Prevalence (per 100 000)")+
  scale_x_discrete("Age", labels = c("","25","","35",
                                     "","45","","55",
                                     "","65","","75",
                                     "","85+"
  ))+theme(legend.position="none")
ggsave("plots/prevalence-by-age-ethnicity-with-uncertainty-bands.pdf",width=4,height=3,device=cairo_pdf)

### Age-spcific incidence plot by ethnicity
ggplot(subset(sim_df,type=='Incidence'),
       aes(x=agecut5,y=middle,group=Ethnicity,colour=Ethnicity))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Ethnicity,colour=NULL,linetype=NA),show.legend=FALSE,alpha=0.3)+
  geom_point()+geom_line()+ylab("Incidence (per 100 000 person-years)")+
  scale_x_discrete("Age", labels = c("","25","","35",
                                     "","45","","55",
                                     "","65","","75",
                                     "","85+"
  ))+ theme(legend.position=c(.2,.7),legend.key = element_blank(),legend.title = element_blank())
ggsave("plots/incidence-by-age-ethnicity-with-uncertainty-bands.pdf",width=4,height=3,device=cairo_pdf)

#### HRC B&W Plot
#### NOT in paper - B&W plot, single panel, no uncertinaty bands
ggplot(subset(sim_df,type=='Prevalence'),
       aes(x=agecut5,y=middle,group=Ethnicity,shape=Ethnicity))+
  #geom_ribbon(aes(ymin=lower,ymax=upper,fill=Ethnicity,colour=NULL),alpha=0.3)+
  geom_line()+
  geom_point(size=2,fill="white")+ylab("Prevalence (per 100 000)")+
  scale_x_discrete("Age", labels = c("","25","","35",
                                     "","45","","55",
                                     "","65","","75",
                                     "","85+"))+
  scale_shape_manual(values=c(21,22,24,23))+
  theme_bw()+ theme(legend.position=c(.2,.7),legend.key = element_blank(),legend.title = element_blank())
ggsave("plots/hrc-prevalence-by-age-ethnicity-with-uncertainty-bands.pdf",width=4,height=3,device=cairo_pdf)


## Using bars rather than uncertainty bands
ggplot(subset(sim_df,type=='Prevalence'),
       aes(x=agecut5,y=middle,group=Ethnicity,colour=Ethnicity))+
  geom_errorbar(aes(ymax = upper, ymin=lower), width=0)+
  geom_point()+geom_line()+ylab("Prevalence (per 100,000)")+
  xlab("Age")
ggsave("prevalence-by-age-ethnicity-with-uncertainty.pdf",width=5,height=3,device=cairo_pdf)

## Using bars rather than uncertainty bands
ggplot(subset(sim_df,type=='Incidence'),
       aes(x=agecut5,y=middle,group=Ethnicity,colour=Ethnicity))+
  geom_errorbar(aes(ymax = upper, ymin=lower), width=0)+
  geom_point()+geom_line()+ylab("Incidence (per 100,000 per year)")+
  xlab("Age")
ggsave("incidence-by-age-ethnicity-with-uncertainty.pdf",width=5,height=3,device=cairo_pdf)


##########################
## Calibration Curves
## For illustration of methods

sim_summary <- summary(fit,
                       pars=c("calibration"),
                       probs=c("0.025","0.975"))$summary
estimated_values <- sim_summary[,c("mean", "2.5%", "97.5%")]

# Assesmble a data frame to pass to ggplot()
calibration <- data.frame(type = c(rep("Very Probable",14*6*2),rep("Probable",14*6*2),
                              rep("Possible",14*6*2),rep("Unlikely",14*6*2)),
                      Ethnicity = rep(c(rep("European",14*2),rep("Asian",14*2),rep("Maori",14*2),
                                        rep("Other",14*2),rep("Pasifika",14*2),rep("Unknown",14*2)),4),
                      Sex = rep(c(rep("Female",14),rep("Male",14)),6*4),
                      agecut5 = rep(ages_85,2*6*4)
                     )  %>%
  mutate(middle = estimated_values[, "mean"]) %>%
  mutate(lower = estimated_values[, "2.5%"]) %>%
  mutate(upper = estimated_values[, "97.5%"]) %>%
  mutate(agecut5 = fct_relevel(agecut5,ages_85)) %>%
  filter(!Ethnicity=='Other') %>%
  filter(!Ethnicity=='Unknown') %>%
  mutate(Ethnicity = fct_relevel(Ethnicity,c('European','Asian','Pasifika','Maori'))) %>%
  mutate(type = fct_relevel(type,c('Very Probable','Probable','Possible','Unlikely'))) %>%
  droplevels()

#levels(sim_df$Ethnicity) <- c('Asian','European','Māori','Pasifika')
levels(calibration$Ethnicity) <- c('European','Asian','Pasifika','Māori')

write.csv(calibration,"pd_probs_by_vars.csv")

ggplot(subset(calibration,Sex=='Female'),
       aes(x=agecut5,y=middle,group=Ethnicity,colour=Ethnicity))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Ethnicity,colour=NULL),alpha=0.3)+facet_grid(type~Ethnicity)+
  geom_point()+geom_line()+ylab("Probability of Parkinson's given medication classification,\nethnicity, age, and known dignoses")+
  scale_x_discrete("Age", labels = c("","25","","35",
                                     "","45","","55",
                                     "","65","","75",
                                     "","85+"
  ))+theme(legend.position="none")
ggsave("prob-of-pd-with-uncertainty-bands-female.pdf",width=8,height=7,device=cairo_pdf)


#################################
## Overall group probabilities 


class_probs <- Pfirstall 

levels(class_probs$ethnicity) <- c('European','Asian','Māori','Other','Pasifika','Unknown')
levels(class_probs$sex) <- c('Female','Male','Unknown')
levels(class_probs$classification) <- c('Very Probable','Possible','Probable','Unlikely')

class_probs <- class_probs %>%
  select(classification,ethnicity,agecut5_85,sex) %>%
  inner_join(calibration,by=c("classification"="type","ethnicity"="Ethnicity",
                            "agecut5_85"="agecut5","sex"="Sex")) %>%
  mutate(classification = factor(classification),
         ethnicity = factor(ethnicity),
         sex = factor(sex)) %>%
  group_by(classification) %>%
  summarise(mean_prob=mean(middle),lower=mean(lower),upper=mean(upper))




##########################
## Counts by ethnicity,age,sex
## Not in paper, diagnostics.

sim_summary <- summary(fit,
                       pars=c("ethnic_sex_age_count"),
                       probs=c("0.025","0.975"))$summary
estimated_values <- sim_summary[,c("mean", "2.5%", "97.5%")]

sim_df <- data.frame(type = c(rep("Incidence",14*5*2),rep("Prevalence",14*5*2)),
                     Ethnicity = rep(c(rep("European",14*2),rep("Asian",14*2),rep("Maori",14*2),
                                       rep("Pasifika",14*2),rep("Unknown",14*2)),2),
                     Sex = rep(c(rep("Female",14),rep("Male",14)),5*2),
                     agecut5 = rep(ages_85,2*5*2)) %>%
  mutate(middle = estimated_values[, "mean"]) %>%
  mutate(lower = estimated_values[, "2.5%"]) %>%
  mutate(upper = estimated_values[, "97.5%"]) %>%
  mutate(agecut5 = fct_relevel(agecut5,ages_85)) %>%
  droplevels()

#Plot of counts by ethnicity
ggplot(sim_df,
       aes(x=agecut5,y=middle,group=Sex,colour=Sex))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Ethnicity,colour=NULL),alpha=0.3)+facet_grid(type~Ethnicity)+
  geom_point()+geom_line()+ylab("Counts")+scale_y_log10()+
  scale_x_discrete("Age", labels = c("","25","","35",
                                     "","45","","55",
                                     "","65","","75",
                                     "","85+"
  ))+theme(legend.position="none")
ggsave("counts-by-age-sex-type-ethnicity.pdf",width=8,height=7,device=cairo_pdf)


