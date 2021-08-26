##########################################################################################################
# required libraries

library(mice)
library(tidyverse)
library(ROCR)
library(pROC)
library(caret)

##########################################################################################################

# LOADING DATA AND PREPROCESSING VPS DATASET
load("VPS2017.RData")

source("helper_functions.R")

new_col_names <- VPS %>% names %>% tolower
names(VPS) <- new_col_names
VPS <- VPS %>% mutate(worst.glascow.coma.score = worst.glasgow.coma.score, worst.glascow.coma.score.unknown = worst.glasgow.coma.score.unknown)
VPS <- VPS_preprocessing(VPS)
npat=nrow(VPS)

# Adding SIRS and qSOFA criteria
source('add_sirs_criteria.R')
source('add_qPELOD2_criteria.R')
source('add_PALS_criteria.R')

# Selecting non-null data instances

VPS<-VPS[!(is.na(VPS$abnormal_HR)&is.na(VPS$abnormal_RR)&
             is.na(VPS$abnormal_Temp)&is.na(VPS$abnormal_Leukos)),]

print(paste0('Exclusion dropped ',npat-nrow(VPS),' patients'))

VPS_simp = VPS[,c("outcome","age.in.months",
                  "abnormal_HR","abnormal_RR","abnormal_Temp","abnormal_Leukos",
                  "abnormal_Mentation","abnormal_SBP",
                  "pim.2.score","prism.3.score")] %>% mutate(obs_id = seq(1:nrow(.)))


# creating complete population leukos and mentation datasets and respective reference models

VPS_simp_non_missing_mentation <-  VPS_simp %>% select(-abnormal_Leukos) %>% drop_na()
VPS_simp_non_missing_Leukos <- VPS_simp %>% select(-abnormal_Mentation) %>% drop_na()

complete_case_mentation_model <- VPS_simp_non_missing_mentation %>% 
  glm(outcome ~ abnormal_Mentation + abnormal_SBP + abnormal_RR, data = . , family = binomial(link = 'logit'))
complete_case_Leukos_model <- VPS_simp_non_missing_Leukos %>% 
  glm(outcome ~ abnormal_Leukos + abnormal_HR + abnormal_RR + abnormal_Temp, data = . , family = binomial(link = 'logit'))

##########################################################################################################


###########################################################################################################################
# Leukos model

simulate_abnormal_weighted_Leukos <- function(the_complete_df, runs=10, m = 5, prop_missing = 0.5, weights_vector = c(0.5,0.5), missing_variable = "abnormal_Temp") {
  res <- array(NA, dim = c(4, runs, 3,5))
  dimnames(res) <- list(c("logreg", "missing as normal", "missing as abnormal", "missing discarded"),
                        as.character(1:runs),
                        c("estimate", "2.5 %","97.5 %"),
                        c("Intercept", "abnormal_LeukosTRUE", "abnormal_HRTRUE","abnormal_RRTRUE","abnormal_TempTRUE"))
  
  for(run in 1:runs) 
  {
    print(run)
    the_data <- the_complete_df
    
    # selecting weights_Vector% of rows for imputing missingness
    
    the_data <- the_data %>% mutate(observation_index = 1:nrow(.), 
                                    missing_weight = ifelse(.[[missing_variable]],weights_vector[1], weights_vector[2]))
    
    new_df <- the_data %>% sample_frac(size = prop_missing,weight = missing_weight, replace = F)
    
    data_Leukos <- the_data
    
    
    # inducing missingness
    data_Leukos[data_Leukos$observation_index %in% new_df$observation_index, "abnormal_Leukos"] <- NA
    
    
    ################################
    # missing discarded model
    
    missing_discarded_Leukos_model <- glm(outcome ~ abnormal_Leukos + abnormal_HR + abnormal_RR + abnormal_Temp, 
                                          data = data_Leukos %>% drop_na(),
                                          family = binomial(link = 'logit'))
    
    res[4, run,,] <- cbind(missing_discarded_Leukos_model$coefficients, confint(missing_discarded_Leukos_model)) %>% t

    ################################
    
    ################################
    # MICE model
    
    imp_Leukos <- mice(data_Leukos %>%  select(-obs_id ,-observation_index, -missing_weight), m = m,print = FALSE)
    fit <- with(imp_Leukos, glm(outcome ~ abnormal_Leukos + abnormal_HR + abnormal_RR + abnormal_Temp, family = binomial(link = 'logit')))
    tab <- summary(pool(fit), "all", conf.int = TRUE)
    res[1, run,, ] <- tab[, c("estimate", "2.5 %", "97.5 %")] %>% t
    
    ################################    
    
    ################################    
    # missing as normal model
    
    missing_as_normal_Leukos_model <- glm(outcome ~ abnormal_Leukos + abnormal_HR + abnormal_RR + abnormal_Temp, 
                                          data = data_Leukos %>% mutate_if(is.logical, replace_na, replace = FALSE),
                                          family = binomial(link = 'logit'))
    
    res[2, run,,] <- cbind(missing_as_normal_Leukos_model$coefficients, confint(missing_as_normal_Leukos_model)) %>% t

    ################################   
    
    ################################   
    # missing as abnormal model
    
    missing_as_abnormal_Leukos_model <- glm(outcome ~ abnormal_Leukos + abnormal_HR + abnormal_RR + abnormal_Temp, 
                                            data = data_Leukos %>% mutate_if(is.logical, replace_na, replace = TRUE),
                                            family = binomial(link = 'logit'))
    
    res[3, run,,] <- c(missing_as_abnormal_Leukos_model$coefficients, confint(missing_as_abnormal_Leukos_model)) %>% t
    ################################ 
    
    
  }
  return(list(Leukos = res))
}

prop_missing_list <- c(0.1,0.3,0.5,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95), c(0.1,0.9), c(0.2,0.8), c(0.3,0.7), c(0.4,0.6), c(0.5,0.5),
                            c(0.6,0.4), c(0.7,0.3), c(0.8,0.2), c(0.9,0.1), c(0.95,0.05))

leukos_coeffs <- lapply(prop_missing_list, function(x){
  pmap(list(rep(x, length(weights_vector_list)), weights_vector_list), function(a,b){
    simulate_abnormal_weighted_Leukos(VPS_simp_non_missing_Leukos,
                                      runs = 150, m=5, prop_missing = a, weights_vector = b, missing_variable = "abnormal_Temp")
  })
}) 
saveRDS(leukos_coeffs, "Leukos_Coefficients.RDS")

###########################################################################################################################


###########################################################################################################################
# Mentation model

simulate_abnormal_weighted_Mentation <- function(the_complete_df, runs=10, m = 5, prop_missing = 0.5, weights_vector = c(0.5,0.5), missing_variable = "abnormal_SBP") {
  res <- array(NA, dim = c(4, runs, 3,4))
  dimnames(res) <- list(c("logreg", "missing as normal", "missing as abnormal", "missing discarded"),
                        as.character(1:runs),
                        c("estimate", "2.5 %","97.5 %"),
                        c("Intercept", "abnormal_MentationTRUE","abnormal_SBPTRUE","abnormal_RRTRUE"))
  
  for(run in 1:runs) 
  {
    print(run)
    
    the_data <- the_complete_df
    
    # selecting weights_Vector% of rows for imputing missingness
    
    the_data <- the_data %>% mutate(observation_index = 1:nrow(.), 
                                    missing_weight = ifelse(.[[missing_variable]],weights_vector[1], weights_vector[2]))
    
    new_df <- the_data %>% sample_frac(size = prop_missing,weight = missing_weight, replace = F)
    data_Mentation <- the_data
    
    # inducing missingness
    
    data_Mentation[data_Mentation$observation_index %in% new_df$observation_index, "abnormal_Mentation"] <- NA
    
    ################################
    # missing discarded model
    
    missing_discarded_mentation_model <- glm(outcome ~ abnormal_Mentation + abnormal_SBP + abnormal_RR, 
                                             data = data_Mentation %>% drop_na(),
                                             family = binomial(link = 'logit'))
    
    res[4, run,,] <- cbind(missing_discarded_mentation_model$coefficients, confint(missing_discarded_mentation_model)) %>% t

    ################################
    
    ################################
    # MICE model
    imp_Mentation <- mice(data_Mentation %>%  select(-obs_id ,-observation_index, -missing_weight), m = m,print = FALSE)
    fit <- with(imp_Mentation, glm(outcome ~ abnormal_Mentation + abnormal_SBP + abnormal_RR, family = binomial(link = 'logit')))
    tab <- summary(pool(fit), "all", conf.int = TRUE)
    res[1, run,, ] <- tab[, c("estimate", "2.5 %", "97.5 %")] %>% t
    ################################
    
    ################################
    # missing as normal model
    missing_as_normal_mentation_model <- glm(outcome ~ abnormal_Mentation + abnormal_SBP + abnormal_RR, 
                                             data = data_Mentation %>% mutate_if(is.logical, replace_na, replace = FALSE),
                                             family = binomial(link = 'logit'))
    
    res[2, run,,] <- cbind(missing_as_normal_mentation_model$coefficients, confint(missing_as_normal_mentation_model)) %>% t
    ################################
    
    ################################
    # missing as abnormal model
    missing_as_abnormal_mentation_model <- glm(outcome ~ abnormal_Mentation + abnormal_SBP + abnormal_RR,  
                                               data = data_Mentation %>% mutate_if(is.logical, replace_na, replace = TRUE),
                                               family = binomial(link = 'logit'))
    
    res[3, run,,] <- c(missing_as_abnormal_mentation_model$coefficients, confint(missing_as_abnormal_mentation_model)) %>% t
    ################################
    
  }
  return(list(Mentation = res))
}

###########################################################################################################################

prop_missing_list <- c(0.1,0.3,0.5,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95), c(0.1,0.9), c(0.2,0.8), c(0.3,0.7), c(0.4,0.6), c(0.5,0.5),
                            c(0.6,0.4), c(0.7,0.3), c(0.8,0.2), c(0.9,0.1), c(0.95,0.05))

mentation_coeffs <- lapply(prop_missing_list, function(x){
  pmap(list(rep(x, length(weights_vector_list)), weights_vector_list), function(a,b){
    simulate_abnormal_weighted_Mentation(VPS_simp_non_missing_mentation,
                                         runs = 150, m=5, prop_missing = a, weights_vector = b, missing_variable = "abnormal_SBP")
  })
}) 
saveRDS(mentation_coeffs, "Mentation_Coefficients.RDS")
