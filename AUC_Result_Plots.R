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

outcomes <- VPS_simp_non_missing_Leukos[,'outcome']

VPS_simp_non_missing_Leukos <- VPS_simp_non_missing_Leukos[,c('abnormal_Leukos','abnormal_HR','abnormal_RR','abnormal_Temp')]

outcomes = factor(outcomes)

auc_mean = 0
for(i in 1:150)
{
  print(i)
  pop_model <- train(y=outcomes,x=VPS_simp_non_missing_Leukos,method="glm",family="binomial",
                     trControl = trainControl(method="cv",number = 5,savePredictions=TRUE,classProbs=TRUE,summaryFunction=twoClassSummary))
  
  auc_cis = ci((pop_model$pred)$obs,(pop_model$pred)$Died)
  auc_mean = auc_mean + auc_cis[2]
}

auc_leukos = auc_mean/150

##############################################################
# result generation for Leukos model

results_table_maker_leukos <- function(list_of_res_tables, list_of_true_parameters){
  map2(list_of_res_tables, list_of_true_parameters, function(x,y){
    CR <- rowMeans(x[,,"2.5 %"] < y & y < x[,,"97.5 %"])
    AUC <- rowMeans(x[,,'estimate'])
    return(data.frame(AUC,CR))
  })
}


df_maker_leukos <- function(RDS_file, a_prop_missing_list, a_weights_vector_list, sampling_criteria = "Temperature based sampling",
                     m = 5, complete_case_model_coefficients){
  RDS_list <- readRDS(RDS_file)
  map2(a_prop_missing_list,1:length(a_prop_missing_list), function(a,b){
    lapply(1:length(a_weights_vector_list), function(k){results_table_maker_leukos(lapply(1:length(complete_case_model_coefficients), function(y){RDS_list[[b]][[k]][['Leukos']]}) %>%
                                                                              setNames(complete_case_model_coefficients %>% names),
                                                                            lapply(1:length(complete_case_model_coefficients),function(y){complete_case_model_coefficients}))}) %>%
      lapply(function(x){lapply(x, function(y){as.data.frame(y) %>% rownames_to_column(var = "Imputation Method")})}) %>%
      lapply(function(x){do.call("rbind", x)}) %>% map2(a_weights_vector_list, function(x,y){x %>% mutate(a_weights_vector_list = rep(1, nrow(.)))}) %>%
      do.call("rbind", .) %>% gather(key = "Criteria", value = "Calculation", CR, AUC) %>%
      mutate(m = rep(m, nrow(.)),prop_missing = rep(a, nrow(.)), simulation = rep(sampling_criteria, nrow(.)))}) %>%
    do.call("rbind",.)
}

prop_missing_list <- c(0.1,0.3,0.5,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95),c(0.1,0.9),c(0.2,0.8), c(0.3,0.7),c(0.4,0.6),c(0.5,0.5),c(0.6,0.4),c(0.7,0.3),c(0.8,0.2),c(0.9,0.1),c(0.95,0.05))

sbp_mentation_1_df <- df_maker_leukos("AUC_CIs_150Runs_Leukos.RDS", 
                               prop_missing_list, weights_vector_list, sampling_criteria = "Abnormal Temperature Sampling Weight", m=5,auc_leukos)      


for (i in 1:length(c(0.1,0.3,0.5,0.7,0.9)))
{
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='logreg'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='missing as normal'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='missing as abnormal'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='missing discarded'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
}

sbp_mentation_1_df['prop_missing'] = lapply(sbp_mentation_1_df['prop_missing'], as.numeric)
sbp_mentation_1_df['Calculation'] = lapply(sbp_mentation_1_df['Calculation'], as.numeric)
sbp_mentation_1_df['a_weights_vector_list'] = lapply(sbp_mentation_1_df['a_weights_vector_list'], as.numeric)

Leukos_final_df <- (sbp_mentation_1_df) %>% 
  mutate(`Imputation Method` = ifelse(`Imputation Method` == "logreg", "MICE", `Imputation Method`))

####################

###
# MENTATION MODEL 
####################

VPS_simp_non_missing_mentation <- VPS_simp %>% drop_na(abnormal_Mentation,abnormal_SBP,abnormal_RR) %>% select(-abnormal_Leukos) %>% drop_na()

outcomes <- VPS_simp_non_missing_mentation[,'outcome']

VPS_simp_non_missing_mentation <- VPS_simp_non_missing_mentation[,c('abnormal_Mentation','abnormal_SBP','abnormal_RR')]

outcomes = factor(outcomes)
auc_mean = 0

for(i in 1:150)
{
  print(i)
  pop_model <- train(y=outcomes,x=VPS_simp_non_missing_mentation,method="glm",family="binomial",
                     trControl = trainControl(method="cv",number = 5,savePredictions=TRUE,classProbs=TRUE,summaryFunction=twoClassSummary))

  auc_cis = ci((pop_model$pred)$obs,(pop_model$pred)$Died)
  auc_mean = auc_mean + auc_cis[2]
}

auc_mentation = auc_mean/150

print(auc_mentation)

###############################

results_table_maker <- function(list_of_res_tables, list_of_true_parameters){
  map2(list_of_res_tables, list_of_true_parameters, function(x,y){
    CR <- rowMeans(x[,,"2.5 %"] < y & y < x[,,"97.5 %"])
    AUC <- rowMeans(x[,,'estimate'])
    return(data.frame(AUC,CR))
  })
}



df_maker <- function(RDS_file, a_prop_missing_list, a_weights_vector_list, sampling_criteria = "Temperature based sampling",
                     m = 5, complete_case_model_coefficients){
  RDS_list <- readRDS(RDS_file)
  map2(a_prop_missing_list,1:length(a_prop_missing_list), function(a,b){
    lapply(1:length(a_weights_vector_list), function(k){results_table_maker(lapply(1:length(complete_case_model_coefficients), function(y){RDS_list[[b]][[k]][['Mentation']]}) %>%
                                                                              setNames(complete_case_model_coefficients %>% names),
                                                                            lapply(1:length(complete_case_model_coefficients),function(y){complete_case_model_coefficients}))}) %>%
      lapply(function(x){lapply(x, function(y){as.data.frame(y) %>% rownames_to_column(var = "Imputation Method")})}) %>%
      lapply(function(x){do.call("rbind", x)}) %>% map2(a_weights_vector_list, function(x,y){x %>% mutate(a_weights_vector_list = rep(1, nrow(.)))}) %>%
      do.call("rbind", .) %>% gather(key = "Criteria", value = "Calculation", AUC,CR) %>%
      mutate(m = rep(m, nrow(.)),prop_missing = rep(a, nrow(.)), simulation = rep(sampling_criteria, nrow(.)))}) %>%
    do.call("rbind",.)
}


prop_missing_list <- c(0.1,0.3,0.5,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95),c(0.1,0.9),c(0.2,0.8), c(0.3,0.7),c(0.4,0.6),c(0.5,0.5),c(0.6,0.4),c(0.7,0.3),c(0.8,0.2),c(0.9,0.1),c(0.95,0.05))

sbp_mentation_1_df <- df_maker("AUC_CIs_150Runs_Mentation.RDS", 
                               prop_missing_list, weights_vector_list, sampling_criteria = "SBP based sampling", m=5, auc_mentation)      

x = sbp_mentation_1_df
for (i in 1:length(c(0.1,0.3,0.5,0.7,0.9)))
{
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='logreg'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='missing as normal'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='missing as abnormal'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
  sbp_mentation_1_df[((sbp_mentation_1_df['prop_missing']==prop_missing_list[i]) & sbp_mentation_1_df['Imputation Method']=='missing discarded'),'a_weights_vector_list']= c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
}

sbp_mentation_1_df['prop_missing'] = lapply(sbp_mentation_1_df['prop_missing'], as.numeric)
sbp_mentation_1_df['Calculation'] = lapply(sbp_mentation_1_df['Calculation'], as.numeric)
sbp_mentation_1_df['a_weights_vector_list'] = lapply(sbp_mentation_1_df['a_weights_vector_list'], as.numeric)


####################

Mentation_final_df <- (sbp_mentation_1_df) %>% 
  mutate(`Imputation Method` = ifelse(`Imputation Method` == "logreg", "MICE", `Imputation Method`))

#######################################################

# to save the final plot
#######################################################


final_df = rbind(Mentation_final_df,Leukos_final_df)

final_df[final_df['simulation']=='Abnormal Temperature Sampling Weight','simulation'] = 'pSIRS Model'
final_df[final_df['simulation']=='SBP based sampling','simulation'] = 'qSOFA Model'

final_auc_plot <- final_df %>%
  filter(Criteria=='AUC'& `Imputation Method` != "missing as abnormal" & `Criteria` != "complete dataset") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing"),
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_string) %>% sort))

hline.data <- data.frame(z = c(auc_leukos,auc_mentation), simulation = c("pSIRS Model", "qSOFA Model")) 

final_plot = ggplot(final_auc_plot,aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`))+ 
  facet_grid(prop_missing_factor~simulation)+ 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + geom_line(position = position_dodge(width = 0.05)) +
  geom_hline(aes(yintercept = z,linetype="complete dataset"), hline.data)+
  scale_linetype_manual(name = " ",values = 1) +
  ylim(0, 1) +
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) +
  labs(y="AUC",x = "Abnormal Temperature Sampling Weight        Abnormal SBP Sampling Weight", colour = "Imputation Method") +
  theme_bw() +
  theme(legend.title=element_text(face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(),axis.text.x = element_text(angle=90, hjust=1))+
  ggtitle("pSIRS and qSOFA Models - Area under Curve (AUC)") +
  theme(plot.title = element_text(hjust = 0.5))


final_plot

library(viridis)
ggsave('AUC_plot.pdf',  final_plot+scale_color_viridis(discrete = TRUE), width = 8, height = 8,dpi=100)



########################################################################
# To compute p-values
########################################################################

props_req = c(0.1,0.3,0.5,0.7,0.9)
weights_vector = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
for (i in 1:length(props_req))
{
  for (j in 1:length(weights_vector))
  {
    Mentation_final_df = rbind(Mentation_final_df,c('complete dataset',weights_vector[j],'AUC',auc_mentation,5,props_req[i],'SBP based sampling'))
  }
}

props_req = c(0.1,0.3,0.5,0.7,0.9)
weights_vector = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
for (i in 1:length(props_req))
{
  for (j in 1:length(weights_vector))
  {
    Leukos_final_df = rbind(Leukos_final_df,c('complete dataset',weights_vector[j],'AUC',auc_leukos,5,props_req[i],'Abnormal Temperature Sampling Weight'))
  }
}

##########################

Leukos_final_df <- Leukos_final_df %>% filter(Criteria=='AUC')
Mentation_final_df <- Mentation_final_df %>% filter(Criteria=='AUC')

Leukos_final_df['prop_missing'] = lapply(Leukos_final_df['prop_missing'], as.numeric)
Leukos_final_df['Calculation'] = lapply(Leukos_final_df['Calculation'], as.numeric)
Leukos_final_df['a_weights_vector_list'] = lapply(Leukos_final_df['a_weights_vector_list'], as.numeric)

Mentation_final_df['prop_missing'] = lapply(Mentation_final_df['prop_missing'], as.numeric)
Mentation_final_df['Calculation'] = lapply(Mentation_final_df['Calculation'], as.numeric)
Mentation_final_df['a_weights_vector_list'] = lapply(Mentation_final_df['a_weights_vector_list'], as.numeric)


paired_t_test_calculator <- function(a_df,actual_method, alternative_method = "missing discarded", the_criteria = "RMSE", the_simulation = "Temperature based sampling"){
  a_df %>% filter(simulation == the_simulation  & 
                    (`Imputation Method` == actual_method | `Imputation Method` == alternative_method)) %>% 
    group_by(prop_missing,simulation, Criteria) %>% nest() %>%
    mutate(
      t_test_df = lapply(data, function(x) {
        spread(x, `Imputation Method`, Calculation)})
    ) %>% mutate(t_test_p_value = lapply(t_test_df, function(y){
      t.test(y[[actual_method]], y[[alternative_method]], paired = T, alternative = "two.sided")$p.value
    }) %>% unlist) %>% 
    mutate(Average_difference_from_MICE = lapply(t_test_df, function(z){
      mean(z[[alternative_method]] - z[[actual_method]], na.rm = T)
    }) %>% unlist) %>% filter(Criteria == the_criteria) %>% arrange(prop_missing)
}

Leukos_RMSE_df <- paired_t_test_calculator(Leukos_final_df, actual_method="MICE",alternative_method = "complete dataset", the_criteria = "AUC", the_simulation = "Abnormal Temperature Sampling Weight")
Mentation_RMSE_df <- paired_t_test_calculator(Mentation_final_df, actual_method="MICE",alternative_method = "complete dataset", the_criteria = "AUC", the_simulation = "SBP based sampling")

disc_Leukos_RMSE_df <- paired_t_test_calculator(Leukos_final_df, actual_method = 'missing discarded',alternative_method = "complete dataset", the_criteria = "AUC", the_simulation = "Abnormal Temperature Sampling Weight")
disc_Mentation_RMSE_df <- paired_t_test_calculator(Mentation_final_df, actual_method = 'missing discarded',alternative_method = "complete dataset", the_criteria = "AUC", the_simulation = "SBP based sampling")

normal_Leukos_RMSE_df <- paired_t_test_calculator(Leukos_final_df, actual_method = 'missing as normal',alternative_method = "complete dataset", the_criteria = "AUC", the_simulation = "Abnormal Temperature Sampling Weight")
normal_Mentation_RMSE_df <- paired_t_test_calculator(Mentation_final_df, actual_method = 'missing as normal',alternative_method = "complete dataset", the_criteria = "AUC", the_simulation = "SBP based sampling")



(Leukos_RMSE_df %>% filter( t_test_p_value < 0.01) %>% nrow) / (Leukos_RMSE_df %>%  nrow)
(Mentation_RMSE_df %>% filter( t_test_p_value < 0.01) %>% nrow)/(Mentation_RMSE_df %>%  nrow)

(disc_Leukos_RMSE_df %>% filter( t_test_p_value < 0.01) %>% nrow) / (disc_Leukos_RMSE_df %>%  nrow)
(disc_Mentation_RMSE_df %>% filter( t_test_p_value < 0.01) %>% nrow)/(disc_Mentation_RMSE_df %>%  nrow)

(normal_Leukos_RMSE_df %>% filter( t_test_p_value < 0.01) %>% nrow) / (normal_Leukos_RMSE_df %>%  nrow)
(normal_Mentation_RMSE_df %>% filter( t_test_p_value < 0.01) %>% nrow)/(normal_Mentation_RMSE_df %>%  nrow)
