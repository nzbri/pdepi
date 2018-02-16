// Prevalence and Incidence calculations of Parkinson's by ethnicity in NZ
// 
// Calibration model taking into account medication classification, age, sex, and 
// ethnicity. This uses a Bernoulli model.
//
// Counts derived from calibration model then go into a prevalence incidence model.
// This uses a normal approximation to a binomial model. Can't use a binomial model directly
// as can't cast real parameters back to integers. 
//
// Author: Daniel Myall <daniel.myall@nzbri.org>
// Any queries or feedback on issues with model welcome.
//
// v1.0 : 2016-10-05 - Initial model with diagnositic uncertainty only
// v1.1 : 2017-10-10 - Add full normal approximation to binomial model

data {
  
  
  // Data for calibration model
  
  int<lower=0> N; // Number of known diagnoses
  int<lower=0> n_beta; // number of predictors (intercept, sex, medication class)
  int<lower=0> n_ethnicity; // European, Asian, Maori, Other, Pacific, Unknown
  int<lower=0> n_sex; // Female, male
  int<lower=0> n_age; // age group 20, 25, ..., 80, 85+
  int<lower=0> n_pdclass; // medication class: v prob, prob, poss, unlikely
  int<lower=0> n_types; // incidence, prevalence
  
  // Ethicity, agegroup, and medication class of each datapoint
  int<lower=0,upper=n_ethnicity> ethnicity[N];
  int<lower=0,upper=n_age> age[N];
  int<lower=0,upper=n_pdclass> pdclass[N];
  
  matrix[N, n_beta] X; // Predictor matrix (intercept, sex, medication class)

  // Diagnosis: PD == 1 / Other == 0
  int<lower=0,upper=1> y[N];
  
  // Data for prevalence/incidence model
  
  // Number of people in each stratification:
  // (type,classification, ethnicity, sex, age)
  int counts_cesa[n_types, n_pdclass, n_ethnicity, n_sex, n_age];
  
  // Population by types, ethnicity
  int pop_e[n_types,n_ethnicity-1];
  
  // Population by types, ethnicity, sex, and age
  int pop_esa[n_types, n_ethnicity-1, n_sex, n_age];
  
  // NZ standard population fractions by sex and age
  real pop_nzs[n_sex, n_age];
  real pop_nzsa[n_age];
  
  // WHO standard population fractions by age
  real pop_whos[n_age];
} 

parameters {

  // sex and class effects
  vector[n_beta] beta;

  // Mean for age and ethnicity effect by classification
  vector[n_age-1] a_raw[n_pdclass];
  vector[n_ethnicity] e[n_pdclass];
  
  // Standard deviation for age and ethnicity by classification
  vector<lower=0>[n_pdclass] sigma_a;
  vector<lower=0>[n_pdclass] sigma_e;
  
  // Modelled proportion/rate by ethnicity
  real theta_mean[n_types];
  real theta_re[n_types, n_ethnicity-1];
  real<lower=0> sigma_theta[n_types];

  //Modelled proporation/rate by ethnicity, age, sex
  real theta_as_mean[n_types,n_age];
  real theta_as_re[n_types, n_ethnicity-1, n_sex, n_age];
  real<lower=0> sigma_theta_as[n_types,n_age];
  

  
}

transformed parameters {
  
  vector[N] y_hat;
  vector[N] alpha;
  vector[n_age] a[n_pdclass];
  
  // Modelled frequency overall by ethnicity
  real<lower=0,upper=1> theta[n_types, n_ethnicity-1];

  // Modelled frequency by ethnicity, age and sex
  real<lower=0,upper=1> theta_as[n_types, n_ethnicity-1, n_sex, n_age];
  
  // Calibration by classification, ethnicity, sex, and age
  real calibration[n_pdclass, n_ethnicity, n_sex, n_age];
  
  // Number of PD by ethnicity, sex, and age
  real ethnic_sex_age_count[n_types,n_ethnicity,n_sex,n_age];
  real ethnic_sex_age_corrected_count[n_types,n_ethnicity,n_sex,n_age];
  
  // Number of PD by ethnicity
  real ethnic_count[n_types,n_ethnicity];
  real ethnic_corrected_count[n_types, n_ethnicity-1];
  
  // Corrections
  real correction_missing_nhi[n_types];
  real correction_change_missing_nhi;
  real correction_unmedicated[n_types];
  real correction_unknown_ethnicity[n_types];
  real correction_unknown_ethnicity_age[n_types,n_age];
  
  real tmpt;
  real to_per_100000;
  
  to_per_100000 = 100000.0;
  
  // Place constraint so age effects a[i,:] sum to zero for each i
  for (i in 1:n_pdclass) {
    for (j in 1:(n_age-1)) {
      a[i,j] = a_raw[i,j];
    }
    
    a[i,n_age] = -sum(a_raw[i]);
  }
  
  // Partially Pooled effects - calibration model
  
  for (i in 1:N) {
    alpha[i] = e[pdclass[i]][ethnicity[i]]
              + a[pdclass[i]][age[i]];
  }
              
  y_hat = X * beta + alpha;
  
  // Partially Pooled effects - incidence/prevalence model
  
  for (t in 1:n_types) {
    for (j in 1:(n_ethnicity-1)) {
      theta[t,j] = inv_logit(theta_mean[t]+sigma_theta[t]*theta_re[t,j]);
    }
  }

  for (t in 1:n_types) {
    for (j in 1:(n_ethnicity-1)) {
      for (k in 1:n_sex) {
        for (l in 1:n_age) {
          theta_as[t,j,k,l] = inv_logit(theta_as_mean[t,l]
                      +sigma_theta_as[t,l]*theta_as_re[t,j,k,l]);
        }
      }
    }
  }


  
  // Based on observed missing NHI in 2013 for prevalence (about a 0.2% increase in rates)
  // For incidence the missing percent over 2006-2013 was 1.2%
  correction_missing_nhi[1] = 1.0/0.988; // Incidence
  correction_missing_nhi[2] = 1.0/0.998; // Prevalence
  
  // Some of the new cases observed will be due to now being tracked via NHI
  // rather than actual new cases. In 2006-2013 the "mean" change was 0.8%, which means 
  // 0.8% of the prevalent cases should be subtracted from the incident cases
  // Only applies to incidence
  correction_change_missing_nhi = 0.008; 
  
  // No correction for unmedicated for incidence:
  // Previous years unmedicated will become medicated so with constant incidence 
  // over time (which it largely is at a national level) it will be at a steady state
  // and they will match the missing unmedicated for previous time period
  correction_unmedicated[1] = 1; // No correction for incidence
  correction_unmedicated[2] = 1.05; // 5% unmedicated for prevalence
  
  
  // Calculate fraction with Parkinson's by classification, ethnicity, sex, and age
  
  for (i in 1:n_pdclass) {
    for (j in 1:n_ethnicity) {
      for (k in 1:n_sex) {
        for (l in 1:n_age) {
          
          tmpt = beta[1] + // Intercept
              (k-1)*beta[5];  // Sex

          if (i>1) {
            tmpt = tmpt + beta[i]; // PD Classification group effect
          }
          
          // Age and ethnicity effects
          tmpt = tmpt + a[i][l] + e[i][j];
          
          calibration[i,j,k,l] = inv_logit(tmpt);
          
        }
      }
    }
  }
  
  
  
  // Calculate estimated number of PD by ethnicity (no corrections applied at this stage)
  
  for (t in 1:n_types) {
    for (j in 1:n_ethnicity) {
      ethnic_count[t,j] = 0;
        for (k in 1:n_sex) {
          for (l in 1:n_age) {
            ethnic_sex_age_count[t,j,k,l] = 0;
          }
        }
    }
  }
  
  for (t in 1:n_types) {
      for (j in 1:n_ethnicity) {
        for (k in 1:n_sex) {
          for (l in 1:n_age) {
            for (i in 1:n_pdclass) {
            
            ethnic_count[t,j] = ethnic_count[t,j] + 
                              calibration[i,j,k,l] * counts_cesa[t,i,j,k,l];
            
            ethnic_sex_age_count[t,j,k,l] = ethnic_sex_age_count[t,j,k,l] +
                              calibration[i,j,k,l] * counts_cesa[t,i,j,k,l];
            
          }
        }
      }
    }
  }
  
  // Calculate unknown ethnicity correction - assume uknowns equally distibuted
  // proportion to counts
  
  for (t in 1:n_types) {
    correction_unknown_ethnicity[t] = 1.0 + 
                ethnic_count[t,5]/sum(ethnic_count[t,1:4]);
                
      for (l in 1:n_age) {
        correction_unknown_ethnicity_age[t,l] = 1.0 + 
                sum(ethnic_sex_age_count[t,5,1:2,l]) /
                   ( sum(ethnic_sex_age_count[t,1:4,1,l]) + 
                        sum(ethnic_sex_age_count[t,1:4,2,l]) );
      }
  }

  // Calculate corrected counts

  for (t in 1:n_types) {
    for (j in 1:(n_ethnicity-1)) {
      
      ethnic_corrected_count[t,j] = ethnic_count[t,j]*
                                  correction_unmedicated[t]*
                                  correction_missing_nhi[t]*
                                  correction_unknown_ethnicity[t];
      
      // Correction for changes in missing data on incidence
      if (t == 1) {
        ethnic_corrected_count[1,j] = ethnic_corrected_count[1,j] - 
                  correction_change_missing_nhi*ethnic_count[2,j];
                  
      }
      
      for (l in 1:n_age) {
        for (k in 1:n_sex) {
            ethnic_sex_age_corrected_count[t,j,k,l] = ethnic_sex_age_count[t,j,k,l]*
                                                correction_unmedicated[t]*
                                                correction_missing_nhi[t]*
                                                correction_unknown_ethnicity_age[t,l];

            // Correction for changes in missing data on incidence
            if (t == 1) {
              ethnic_sex_age_corrected_count[1,j,k,l] = ethnic_sex_age_corrected_count[1,j,k,l] -
                        correction_change_missing_nhi*ethnic_sex_age_count[2,j,k,l];
            }
        }
      }
    }
  }
  
} 


model {
  
  // For choice of priors see
  // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  // Student t seems to be the current recommendation
  // Cauchy prior gave similar results but sampling was much slower.
  
  
  
  // Calibration model - 
  // Probability of having Parkinson's given medication classification, sex, age, ethnicity
  
  for (k in 1:n_pdclass) {
    
    // Handle using priors that age_group_{i} will be correlated with 
    // age_group_{i+1} and can't be considered exchangeable. 

    for (j in 1:(n_age-1)) {

      if (j == 1) {
        a_raw[k,1] ~ student_t (3,a_raw[k,2], sigma_a[k]);
      }
      else if (j == (n_age-1)) {
        a_raw[k,n_age-1] ~ student_t (3,(a_raw[k,n_age-2]+a[k,n_age])/2, sigma_a[k]);
      }
      else {
        a_raw[k,j] ~ student_t (3,(a_raw[k,j-1]+a_raw[k,j+1])/2, sigma_a[k]);
      }
    }

    // Standard partial-pooling for ethnicity effects by medication classification
    // Exchangeability assumption valid here

    e[k] ~ student_t (3,0, sigma_e[k]);

    // Priors for variance
    // Note that lower=0 constraint on sigmas gives a half student-t distribution

    sigma_a[k] ~ student_t(4,0,2);
    sigma_e[k] ~ student_t(4,0,2);
    
  }
  
  beta ~ student_t(3,0,4);

  y ~ bernoulli_logit(y_hat);
  
  
  
  // Prevalence / Incidence model
  
  // Set priors on model parameters
  
  for (t in 1:n_types) {
    
    sigma_theta[t] ~ student_t(4,0,2);
    theta_mean[t] ~ normal(-8,2);
    
    for (l in 1:n_age) {
      
      sigma_theta_as[t,l] ~ student_t(4,0,2);
      theta_as_mean[t,l] ~ normal(-8,2);
      
    }
    
  }
  
  
  
  for (t in 1:n_types) {
    for (j in 1:(n_ethnicity-1)) {
      
      // Calculation of observed prevalence/incidence
      
      theta_re[t,j] ~ normal(0,1);
      
      // Normal approximation to bionomial model
      // Add a constant (0.1) to the standard deviation to add extra uncertainty
      // and overcome model convergence issues when low proporation/rates
      // ethnic_corrected_count is a linear transformation of parameters so 
      // there no requirement for an adjsutment of the jacobian
      
      ethnic_corrected_count[t,j] ~ normal(pop_e[t,j]*theta[t,j],
                            sqrt(pop_e[t,j]*theta[t,j]*(1-theta[t,j]))+0.1);

      for (l in 1:n_age) {

        for (k in 1:n_sex) {
          
          // Calculation of age-sex-specific prevalence/incidence
          
          theta_as_re[t,j,k,l] ~ normal(0,1);

          // Normal approximation to bionomial model
          // Add a constant (0.1) to the standard deviation to add extra uncertainty
          // and overcome model convergence issues when low proporation/rates          
          // ethnic_sex_age_corrected_count is a linear transformation of parameters so
          // there no requirement for an adjsutment of the jacobian
          
          ethnic_sex_age_corrected_count[t,j,k,l] ~ normal(pop_esa[t,j,k,l]*theta_as[t,j,k,l],
              sqrt(pop_esa[t,j,k,l]*theta_as[t,j,k,l]*(1-theta_as[t,j,k,l]))+0.1);
        }

      }
    }
  }
  

}

generated quantities {
  
  // For temporary calcs
  real tmp;
  real tmp2;
  real tmp3;
  
  real total_count[n_types];

  // Observed rates by ethnicity
  // These only incorporate diagnostic uncertainty and thus can't be generalised beyond population
  // These are only for diagnostic purposes are are not used in any calculations
  real ethnic_observed_rates[n_types, n_ethnicity-1];
  
  // Observered rates per 100 000 modelled using a normal approximation to binomial model
  real ethnic_observed_rates_normapprox[n_types,n_ethnicity-1];
  
  // Ethnic rates by sex and age
  // These only incorporate diagnostic uncertainty and thus can't be generalised beyond population
  // These are only for diagnostic purposes are are not used in any calculations
  real ethnic_rates_sex_age[n_types, n_ethnicity-1, n_sex, n_age];
  
  // Age-sex-specific rates per 100 000 modelled using a normal approximation to binomial model
  real ethnic_rates_sex_age_normapprox[n_types,n_ethnicity-1, n_sex, n_age];
  
  // Ethnic sex-standardised rates by age
  real ethnic_rates_age[n_types, n_ethnicity-1, n_age];
  
  // Age-sex-standardised rates by ethnicity
  real ethnic_standardised_rates[n_types, n_ethnicity-1];
  real ethnic_standardised_rates_who[n_types, n_ethnicity-1];
  
  real ethnic_standardised_rates_sex[n_types, n_ethnicity-1, n_sex];
  real ethnic_standardised_ratios_sex[n_types, n_ethnicity-1];
  
  real mean_standardised_duration[n_ethnicity-1];
  real mean_observed_duration[n_ethnicity-1];  

  // Observed rates - by ethnicity
  
  for (t in 1:n_types) {
    
    total_count[t] = 0;
    
    for (j in 1:(n_ethnicity-1)) {
     
      total_count[t] = total_count[t] + ethnic_corrected_count[t,j];
                                  
      ethnic_observed_rates[t,j] = ethnic_corrected_count[t,j]*
                                   to_per_100000/
                                   pop_e[t,j];
                                   
      ethnic_observed_rates_normapprox[t,j] = theta[t,j]*to_per_100000;
                                 
    }
  }

  // Observed rates - by ethnicity, sex, age

  for (t in 1:n_types) {  
    for (j in 1:(n_ethnicity-1)) {
      for (k in 1:n_sex) {
        for (l in 1:n_age) {
          
            tmp = 0;
            
            // Sum out medication classification
            for (i in 1:n_pdclass) {
              tmp  = tmp + calibration[i,j,k,l] * counts_cesa[t,i,j,k,l];
            }
                                 
            ethnic_rates_sex_age[t,j,k,l] = tmp*to_per_100000*
                     correction_unmedicated[t]*
                     correction_missing_nhi[t]*
                     correction_unknown_ethnicity[t]/
                     pop_esa[t,j,k,l];
                     
            ethnic_rates_sex_age_normapprox[t,j,k,l] = theta_as[t,j,k,l]*to_per_100000;
          
        }
      }
    }
  }

  // sex standardised rates - by ethnicity, age

 for (t in 1:n_types) {    
    for (j in 1:(n_ethnicity-1)) {
        for (l in 1:n_age) {
               
            tmp = (ethnic_rates_sex_age_normapprox[t,j,1,l]*pop_nzs[1,l]+
               ethnic_rates_sex_age_normapprox[t,j,2,l]*pop_nzs[2,l])/
               (pop_nzs[1,l]+pop_nzs[2,l]);
          
          ethnic_rates_age[t,j,l] = tmp;
          
        }
    }
 }
    
  
  // Age-sex-standardised rates (to overall NZ) - by ethnicity
  
  for (t in 1:n_types) { 
    for (j in 1:(n_ethnicity-1)) {
        
        tmp = 0;
        tmp2 = 0;
        
        for (k in 1:n_sex) {
          
          tmp3 = 0;
          
          for (l in 1:n_age) {
            
              tmp = tmp + ethnic_rates_sex_age_normapprox[t,j,k,l]*pop_nzs[k,l];
              tmp2 = tmp2 + ethnic_rates_sex_age_normapprox[t,j,k,l]*pop_whos[l]/2;
              tmp3 = tmp3 + ethnic_rates_sex_age_normapprox[t,j,k,l]*pop_nzsa[l];
            
          }
          ethnic_standardised_rates_sex[t,j,k] = tmp3;
        }
        
      ethnic_standardised_ratios_sex[t,j] = ethnic_standardised_rates_sex[t,j,2]
                                                /ethnic_standardised_rates_sex[t,j,1];
      
      ethnic_standardised_rates[t,j] = tmp;
      ethnic_standardised_rates_who[t,j] = tmp2;
        
      }
  }
  
  // Mean duration
  
  for (j in 1:(n_ethnicity-1)) {
    
    mean_standardised_duration[j] = ethnic_standardised_rates[2,j]/
                          ethnic_standardised_rates[1,j];
                          
    mean_observed_duration[j] = ethnic_observed_rates[2,j]/
                          ethnic_observed_rates[1,j];
  }
  
  
}
  
