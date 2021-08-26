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


# creating complete population leukos dataset

VPS_simp_non_missing_Leukos <- VPS_simp %>% select(-abnormal_Mentation) %>% drop_na()

##########################################################################################################


##########################################################################################################

# Inducing and imputing datasets under biased sampling (without replacement) scheme using Abnormal Temperature
# Applying 5-fold Cross validation to train the models and compute AUC values including Confidence Intervals

simulate_abnormal_weighted_Leukos <- function(the_complete_df, runs=10, m = 5, prop_missing = 0.5, weights_vector = c(0.5,0.5), missing_variable = "abnormal_Temp") {
  
  # adding results in "res" list
  res <- array(NA, dim = c(3, runs, 3))
  dimnames(res) <- list(c("logreg", "missing as normal", "missing discarded"),
                        as.character(1:runs),
                        c("estimate", "2.5 %","97.5 %"))
  
  for(run in 1:runs) 
  {

    the_data <- the_complete_df
    
    # selecting weights_Vector% of rows for imputing missingness
    
    the_data <- the_data %>% mutate(observation_index = 1:nrow(.), 
                                    missing_weight = ifelse(.[[missing_variable]],weights_vector[1], weights_vector[2]))
    
    new_df <- the_data %>% sample_frac(size = prop_missing,weight = missing_weight, replace = F)
    
    data_Leukos <- the_data
    
    # inducing missingness
    data_Leukos[data_Leukos$observation_index %in% new_df$observation_index, "abnormal_Leukos"] <- NA
    
    
    ##############################
    
    # missing discarded model
    
    data_disc_leukos = data_Leukos %>% drop_na()
    outcomes <- data_disc_leukos[,'outcome']
    
    data_disc_leukos <- data_disc_leukos[,c("abnormal_HR","abnormal_RR","abnormal_Temp","abnormal_Leukos")]
    
    outcomes = factor(outcomes)
    
    disc_model <- train(y=outcomes,x=data_disc_leukos,method="glm",family="binomial",
                        trControl = trainControl(method="cv",number = 5,savePredictions=TRUE,classProbs=TRUE,summaryFunction=twoClassSummary))
    
    auc_cis = ci((disc_model$pred)$obs,(disc_model$pred)$Died)
    res[3,run,] = cbind(auc_cis[2],auc_cis[1],auc_cis[3])
 
    remove(disc_model)
    ##############################
    
    ##############################
    
    # MICE model
    
    imp_Leukos <- mice(data_Leukos %>%  select(-obs_id ,-observation_index, -missing_weight), m = m,print = FALSE)
    complete_imputed <- data.frame(complete(imp_Leukos, "long", inc = FALSE))
    outcomes <- complete_imputed[,'outcome']
    outcomes = factor(outcomes)
    complete_imputed <- complete_imputed[,c("abnormal_HR","abnormal_RR","abnormal_Temp","abnormal_Leukos")]
    imp_model <- train(y=outcomes,x=complete_imputed,method="glm",family="binomial",
                       trControl = trainControl(method="cv",number = 5,savePredictions=TRUE,classProbs=TRUE,summaryFunction=twoClassSummary))
    auc_cis = ci((imp_model$pred)$obs,(imp_model$pred)$Died)
    res[1,run,] = cbind(auc_cis[2],auc_cis[1],auc_cis[3])

    remove(imp_model)
    
    ############################
    # Missing as normal model
    
    data_norm_leukos = data_Leukos %>% mutate_if(is.logical, replace_na, replace = FALSE)
    outcomes <- data_norm_leukos[,'outcome']
    outcomes = factor(outcomes)
    
    data_norm_leukos <- data_norm_leukos[,c("abnormal_HR","abnormal_RR","abnormal_Temp","abnormal_Leukos")]
    
    
    norm_model <- train(y=outcomes,x=data_norm_leukos,method="glm",family="binomial",
                        trControl = trainControl(method="cv",number = 5,savePredictions=TRUE,classProbs=TRUE,summaryFunction=twoClassSummary))
    
    auc_cis = ci((norm_model$pred)$obs,(norm_model$pred)$Died)
    
    
    
    res[2,run,] = cbind(auc_cis[2],auc_cis[1],auc_cis[3])

    remove(norm_model)
    
    ##############################
    
  }
  return(list(Leukos = res))
}

prop_missing_list <- c(0.1,0.3,0.5,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95), c(0.1,0.9), c(0.2,0.8), c(0.3,0.7), c(0.4,0.6), c(0.5,0.5),
                            c(0.6,0.4), c(0.7,0.3), c(0.8,0.2), c(0.9,0.1), c(0.95,0.05))

leukos_cis <- lapply(prop_missing_list, function(x){
  pmap(list(rep(x, length(weights_vector_list)), weights_vector_list), function(a,b){
    simulate_abnormal_weighted_Leukos(VPS_simp_non_missing_Leukos,
                                      runs = 150, m=5, prop_missing = a, weights_vector = b)
  })
}) 

saveRDS(leukos_cis, "AUC_CIs_150Runs_Leukos.RDS")
