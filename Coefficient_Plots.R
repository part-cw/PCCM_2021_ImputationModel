##########################################################################################################
# required libraries

library(mice)
library(tidyverse)
library(ROCR)
library(pROC)
library(caret)
library(abind)
library(pheatmap)
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

############################################################################################################################

##############################################
# results table maker from the coefficients

results_table_maker <- function(list_of_res_tables, list_of_true_parameters){
  map2(list_of_res_tables, list_of_true_parameters, function(x,y){
    print(x)
    print(dim(x))
    RB <- rowMeans(x[,, "estimate"]) - y
    PB <- 100 * abs((rowMeans(x[,, "estimate"]) - y)/ y)
    CR <- rowMeans(x[,, "2.5 %"] < y & y < x[,, "97.5 %"])
    AW <- rowMeans(x[,, "97.5 %"] - x[,, "2.5 %"])
    RMSE <- sqrt(rowMeans((x[,, "estimate"] - y)^2))
    return(data.frame(RB, PB, CR, AW, RMSE))
  })
}
##############################################


###############################################
# results table

df_maker <- function(RDS_file, a_prop_missing_list, a_weights_vector_list, sampling_criteria = "Temperature based sampling", 
                     m = 5, complete_case_model_coefficients){
  RDS_list <- readRDS(RDS_file)
  map2(a_prop_missing_list,1:length(a_prop_missing_list), function(a,b){
    lapply(1:length(a_weights_vector_list), function(k){results_table_maker(lapply(1:length(complete_case_model_coefficients), function(y){RDS_list[[b]][[k]][[1]][,,,y]}) %>% 
                                                                              setNames(complete_case_model_coefficients %>% names),
                                                                            lapply(1:length(complete_case_model_coefficients),function(y){complete_case_model_coefficients[[y]]}))}) %>% 
      lapply(function(x){lapply(x, function(y){as.data.frame(y) %>% rownames_to_column(var = "Imputation Method")})}) %>% 
      lapply(function(x){map2(x,names(x), function(c,d){c %>% mutate(Imputed_Variable = rep(d, nrow(c)))})}) %>% 
      lapply(function(x){do.call("rbind", x)}) %>% map2(a_weights_vector_list, function(x,y){x %>% mutate(a_weights_vector_list = rep(y[[1]], nrow(.)))}) %>% 
      do.call("rbind", .) %>% gather(key = "Criteria", value = "Calculation", RB, PB, CR, AW, RMSE) %>% 
      mutate(m = rep(m, nrow(.)),prop_missing = rep(a, nrow(.)), simulation = rep(sampling_criteria, nrow(.)))}) %>% 
    do.call("rbind",.)
}


###############
#Mentation results


prop_missing_list <- c(0.1,0.3,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95),c(0.1,0.9),c(0.2,0.8), c(0.3,0.7),c(0.4,0.6),c(0.5,0.5),c(0.6,0.4),c(0.7,0.3),c(0.8,0.2),c(0.9,0.1),c(0.95,0.05))


sbp_leukos_1_df <- df_maker("Leukos_Coefficients.RDS", 
              prop_missing_list, weights_vector_list, m=5, complete_case_mentation_model$coefficients)      


###############
# Leukos results

prop_missing_list <- c(0.1,0.3,0.7,0.9)
weights_vector_list <- list(c(0.05,0.95),c(0.1,0.9),c(0.2,0.8), c(0.3,0.7),c(0.4,0.6),c(0.5,0.5),c(0.6,0.4),c(0.7,0.3),c(0.8,0.2),c(0.9,0.1),c(0.95,0.05))

sbp_mentation_1_df <- df_maker("Mentation_Coefficients.RDS",
                            prop_missing_list, weights_vector_list, sampling_criteria = "SBP based sampling", m=5, complete_case_Leukos_model$coefficients)

#################

#################

# Leukos and Mentation dataframe creation for plotting

Leukos_final_df <- sbp_leukos_1_df %>%
  mutate(`Imputation Method` = ifelse(`Imputation Method` == "logreg", "MICE", `Imputation Method`)) %>% mutate(Imputed_Variable = str_replace(Imputed_Variable, "TRUE", "")) %>%
  mutate(Imputed_Variable = str_replace(Imputed_Variable, "abnormal_", "Abnormal "))

Mentation_final_df <- sbp_mentation_1_df %>% 
  mutate(`Imputation Method` = ifelse(`Imputation Method` == "logreg", "MICE", `Imputation Method`)) %>% mutate(Imputed_Variable = str_replace(Imputed_Variable, "TRUE", "")) %>% 
  mutate(Imputed_Variable = str_replace(Imputed_Variable, "abnormal_", "Abnormal "))

#################

# percentage bias plots

Leukos_m_5_PB_plot <- Leukos_final_df %>%
  filter(prop_missing == 0.5 & simulation == "Temperature based sampling" & `Imputation Method` != "missing as abnormal") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing") ,
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_factor))) %>%
  filter(Criteria == "PB") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`)) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + geom_line( position = position_dodge(width = 0.05)) +
  facet_grid(prop_missing_factor~Imputed_Variable, scales = "free_y") + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) +
  labs(y = "Percent Bias",x = "Abnormal Temperature Sampling Weight", title = "SIRS Model Coefficients - Percent Bias", colour = "Imputation Method") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

Leukos_m_5_PB_plot <- Leukos_m_5_PB_plot %>% plot_theme(y_angle = 0, x_angle = 90) 

Mentation_m_5_PB_plot <- Mentation_final_df %>% 
  filter(prop_missing == 0.5 & simulation == "SBP based sampling" & `Imputation Method` != "missing as abnormal") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing") ,
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_factor))) %>%
  filter(Criteria == "PB") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`)) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + geom_line( position = position_dodge(width = 0.05)) +
  facet_grid(prop_missing_factor~Imputed_Variable, scales = "free_y") + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) +
  labs(y = "Percent Bias",x = "Abnormal SBP Sampling Weight", title = "qSOFA Model Coefficients - Percent Bias", colour = "Imputation Method") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

Mentation_m_5_PB_plot <- Mentation_m_5_PB_plot %>% plot_theme(y_angle=0, x_angle=90)

#################


#################
# RMSE plot

Leukos_m_5_RMSE_plot <- Leukos_final_df %>%
  filter(simulation == "Temperature based sampling" & `Imputation Method` != "missing as abnormal") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing") ,
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_factor) %>% sort)) %>%
  filter(Criteria == "RMSE") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`)) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + 
  geom_line( position = position_dodge(width = 0.05)) +
  facet_grid(prop_missing_factor~Imputed_Variable, scales = "free_y") + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) +
  labs(y = "RMSE",x = "Abnormal Temperature Sampling Weight", title = "SIRS Model Coefficients - RMSE", colour = "Imputation Method") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + ylim(0, 0.7) 

#scale_y_continuous(breaks = seq(0,200,5))

Leukos_m_5_RMSE_plot <- Leukos_m_5_RMSE_plot %>% plot_theme(y_angle = 0, x_angle= 90) 


Mentation_m_5_RMSE_plot <- Mentation_final_df %>%
  filter(simulation == "SBP based sampling" & `Imputation Method` != "missing as abnormal") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing"),
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_string) %>% sort)) %>%
  filter(Criteria == "RMSE") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`)) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + geom_line( position = position_dodge(width = 0.05)) +
  facet_grid(prop_missing_factor~Imputed_Variable, scales = "free_y") + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) +
  labs(y = "RMSE",x = "Abnormal SBP Sampling Weight", title = "qSOFA Model Coefficients - RMSE", colour = "Imputation Method") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

#scale_y_continuous(breaks = seq(0,200,5))

Mentation_m_5_RMSE_plot <- Mentation_m_5_RMSE_plot %>% plot_theme(y_angle = 0, x_angle=90) 

#################
# Coverage rate

Mentation_m_5_CR_plot <- Mentation_final_df %>%
  filter(prop_missing == 0.5 & simulation == "SBP based sampling" & `Imputation Method` != "missing as abnormal") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing") ,
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_factor))) %>%
  filter(Criteria == "CR") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`)) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + geom_line( position = position_dodge(width = 0.05)) +
  geom_hline(yintercept = 0.95) +
  facet_grid(prop_missing_factor~Imputed_Variable, scales = "free_y") + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) + 
  scale_y_continuous(breaks = c(0.25,0.50,0.75,0.95,1))+
  labs(y = "Coverage Rate",x = "Abnormal SBP Sampling Weight", title = "qSOFA Model Coefficients - Coverage Rate", colour = "Imputation Method") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

#scale_y_continuous(breaks = seq(0,200,5))

Mentation_m_5_CR_plot <- Mentation_m_5_CR_plot %>% plot_theme(y_angle = 0) 

Leukos_m_5_CR_plot <- Leukos_final_df %>%
  filter(prop_missing == 0.5 & simulation == "Temperature based sampling" & `Imputation Method` != "missing as abnormal") %>%
  mutate(prop_missing_percentage = as.integer(prop_missing * 100),
         prop_missing_string = str_c(prop_missing_percentage, "% missing") ,
         prop_missing_factor = factor(prop_missing_string),
         prop_missing_factor = factor(prop_missing_factor, unique(prop_missing_factor))) %>%
  filter(Criteria == "CR") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method`)) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + geom_line( position = position_dodge(width = 0.05)) +
  geom_hline(yintercept = 0.95) +
  facet_grid(prop_missing_factor~Imputed_Variable, scales = "free_y") + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95) ) + 
  scale_y_continuous(breaks = c(0.25,0.50,0.75,0.95,1))+
  labs(y = "Coverage Rate",x = "Abnormal Temperature Sampling Weight", title = "SIRS Model Coefficients - Coverage Rate", colour = "Imputation Method") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

#scale_y_continuous(breaks = seq(0,200,5))

Leukos_m_5_CR_plot <- Leukos_m_5_CR_plot %>% plot_theme(y_angle = 0) 


combined_m_5_CR_plot <- final_plot_df %>%
  filter(Imputed_Variable == "Mentation" & prop_missing == 0.5) %>% 
  mutate(Number_of_induced_missing = as.integer(prop_missing * nrow(VPS_simp_complete_cases)),
         Number_of_induced_missing_string = str_c(Number_of_induced_missing, " missing") ,
         Number_of_induced_missing_factor = factor(Number_of_induced_missing_string),
         Number_of_induced_missing_factor = factor(Number_of_induced_missing_factor, unique(Number_of_induced_missing_factor))) %>%
  filter(Criteria == "CR") %>%  ggplot(aes(x = a_weights_vector_list, y = Calculation, colour = `Imputation Method` )) + 
  geom_point(size= 2, position = position_dodge(width = 0.05)) + facet_grid(Number_of_induced_missing_factor~simulation) +
  geom_hline(yintercept = 0.95) +
  scale_y_continuous(breaks = seq(0,1.5,0.1)) + scale_x_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9) ) +
  labs(y = "Coverage Rate",x = "Abnormal Sampling Weight", title = "Mentation Coefficient Coverage Rate", colour = "Imputation Method") +
  theme_bw() + theme()

combined_m_5_CR_plot

#################

# Unit ID plot

unit_id_plot <- VPS %>% mutate(WBC_high_and_low_unknown = low.white.blood.cell.count.unknown & high.white.blood.cell.count.unknown) %>%
  group_by(unit.id) %>% 
  summarize(`WBC` = mean(WBC_high_and_low_unknown),`Glascow Coma Score` = mean(worst.glascow.coma.score.unknown), patients_at_unit = n(), 
            mean_length_of_stay = mean(physical.length.of.stay..days.), percentage_abnormal_glascow_coma_score = mean(abnormal_Mentation, na.rm = T),
            percentage_abnormal_WBC = mean(abnormal_Leukos, na.rm = T)) %>%
  gather(key = "Variable", value = "Proportion_Missing", WBC, `Glascow Coma Score` ) %>% 
  mutate(percentage_abnormal = ifelse(Variable == "WBC", percentage_abnormal_WBC, percentage_abnormal_glascow_coma_score)) %>% 
  ggplot(aes(x = patients_at_unit, y = Proportion_Missing, colour = percentage_abnormal)) + geom_point() +
  facet_wrap(~Variable) + theme_bw() + 
  labs(y = "Proportion Missing", x = "Patients at Unit", colour ="Proportion abnormal")

unit_id_plot <- plot_theme(unit_id_plot, y_angle=0, x_angle=90) 
#################

#################
# Heatmap plots

paired_t_test_calculator <- function(a_df, alternative_method = "missing discarded", the_criteria = "RMSE", the_simulation = "Temperature based sampling"){
  a_df %>% filter(simulation == the_simulation  & 
                               (`Imputation Method` == "MICE" | `Imputation Method` == alternative_method)) %>% 
    group_by(a_weights_vector_list, prop_missing, simulation, Criteria) %>% nest() %>%
    mutate(
      t_test_df = lapply(data, function(x) {
        spread(x, `Imputation Method`, Calculation)})
    ) %>% mutate(t_test_p_value = lapply(t_test_df, function(y){
      t.test(y$MICE, y[[alternative_method]], paired = T, alternative = "less")$p.value
    }) %>% unlist) %>% 
    mutate(Average_difference_from_MICE = lapply(t_test_df, function(z){
      mean(z[[alternative_method]] - z$MICE, na.rm = T)
    }) %>% unlist) %>% filter(Criteria == the_criteria) %>% arrange(prop_missing, a_weights_vector_list)
}


Leukos_RMSE_df <- paired_t_test_calculator(Leukos_final_df, alternative_method = "missing discarded", the_criteria = "RMSE", the_simulation = "Temperature based sampling")
Mentation_RMSE_df <- paired_t_test_calculator(Mentation_final_df, alternative_method = "missing discarded", the_criteria = "RMSE", the_simulation = "SBP based sampling")

(Leukos_RMSE_df %>% filter(prop_missing < 0.5 & t_test_p_value < 0.05) %>% nrow) / (Leukos_RMSE_df %>% filter(prop_missing < 0.5) %>% nrow)
(Mentation_RMSE_df %>% filter(prop_missing < 0.5 & t_test_p_value < 0.05) %>% nrow)/(Mentation_RMSE_df %>% filter(prop_missing < 0.5) %>% nrow)




Leukos_RMSE_p_value_matrix <- Leukos_RMSE_df$t_test_p_value %>% matrix(ncol = Leukos_RMSE_df$prop_missing %>% unique %>% length)
rownames(Leukos_RMSE_p_value_matrix) <- Leukos_RMSE_df$a_weights_vector_list %>% unique %>% str_c("sampling weight=", .)
colnames(Leukos_RMSE_p_value_matrix) <- (Leukos_RMSE_df$prop_missing %>% unique * 100) %>% str_c(., "% missing")
psirs_heatmap = pheatmap(Leukos_RMSE_p_value_matrix, cluster_rows = F, cluster_cols = F, show_colnames = T , main = "pSIRS Model - MICE vs Missing Discarded RMSE p-values", 
         legend_breaks = c(0.05, 0.10, 0.15, max(Leukos_RMSE_p_value_matrix)),legend_labels = c("0.05", "0.10", "0.15", "p-value\n"), legend = T)

# ggsave()

Mentation_RMSE_p_value_matrix <- Mentation_RMSE_df$t_test_p_value %>% matrix(ncol = Mentation_RMSE_df$prop_missing %>% unique %>% length)
rownames(Mentation_RMSE_p_value_matrix) <- Mentation_RMSE_df$a_weights_vector_list %>% unique %>% str_c("sampling weight=", .)
colnames(Mentation_RMSE_p_value_matrix) <- (Mentation_RMSE_df$prop_missing %>% unique * 100) %>% str_c(., "% missing")
pheatmap(Mentation_RMSE_p_value_matrix, cluster_rows = F, cluster_cols = F, show_colnames = T , main = "qSOFA Model - MICE vs Missing Discarded RMSE p-values",
         legend_breaks = c(0.05, 0.10, 0.15, max(Mentation_RMSE_p_value_matrix)),legend_labels = c("0.05", "0.10", "0.15", "p-value\n"), legend = T)


Leukos_RMSE_df <- paired_t_test_calculator(Leukos_final_df, alternative_method = "missing as normal", the_criteria = "RMSE", the_simulation = "Temperature based sampling")
Mentation_RMSE_df <- paired_t_test_calculator(Mentation_final_df, alternative_method = "missing as normal", the_criteria = "RMSE", the_simulation = "SBP based sampling")

(Leukos_RMSE_df %>% filter(prop_missing < 0.5 & t_test_p_value < 0.05) %>% nrow) / (Leukos_RMSE_df %>% filter(prop_missing < 0.5) %>% nrow)
(Mentation_RMSE_df %>% filter(prop_missing < 0.5 & t_test_p_value < 0.05) %>% nrow)/(Mentation_RMSE_df %>% filter(prop_missing < 0.5) %>% nrow)


Leukos_RMSE_p_value_matrix <- Leukos_RMSE_df$t_test_p_value %>% matrix(ncol = Leukos_RMSE_df$prop_missing %>% unique %>% length)
rownames(Leukos_RMSE_p_value_matrix) <- Leukos_RMSE_df$a_weights_vector_list %>% unique %>% str_c("sampling weight=", .)
colnames(Leukos_RMSE_p_value_matrix) <- (Leukos_RMSE_df$prop_missing %>% unique * 100) %>% str_c(., "% missing")
psirs_heatmap = pheatmap(Leukos_RMSE_p_value_matrix, cluster_rows = F, cluster_cols = F, show_colnames = T , main = "pSIRS Model - MICE vs Missing as normal RMSE p-values", 
                         legend_breaks = c(0.05, 0.10, 0.15, max(Leukos_RMSE_p_value_matrix)),legend_labels = c("0.05", "0.10", "0.15", "p-value\n"), legend = T)

# ggsave()

Mentation_RMSE_p_value_matrix <- Mentation_RMSE_df$t_test_p_value %>% matrix(ncol = Mentation_RMSE_df$prop_missing %>% unique %>% length)
rownames(Mentation_RMSE_p_value_matrix) <- Mentation_RMSE_df$a_weights_vector_list %>% unique %>% str_c("sampling weight=", .)
colnames(Mentation_RMSE_p_value_matrix) <- (Mentation_RMSE_df$prop_missing %>% unique * 100) %>% str_c(., "% missing")
pheatmap(Mentation_RMSE_p_value_matrix, cluster_rows = F, cluster_cols = F, show_colnames = T , main = "qSOFA Model - MICE vs Missing as normal RMSE p-values",
         legend_breaks = c(0.05, 0.10, 0.15, max(Mentation_RMSE_p_value_matrix)),legend_labels = c("0.05", "0.10", "0.15", "p-value\n"), legend = T)

#################

## FRECHET DISTANCE 
library(SimilarityMeasures)

aaa <- Mentation_final_df %>% 
  filter(simulation == "SBP based sampling" & Imputed_Variable == "Abnormal SBP" & prop_missing == 0.9 & Criteria == "RMSE" & `Imputation Method` == "MICE") %>% 
  select(a_weights_vector_list, Calculation) %>% as.matrix() 

bbb <- Mentation_final_df %>% 
  filter(simulation == "SBP based sampling" & Imputed_Variable == "Abnormal SBP" & prop_missing == 0.9 & Criteria == "RMSE" & `Imputation Method` == "missing discarded") %>% 
  select(a_weights_vector_list, Calculation) %>% as.matrix() 

ccc <- Mentation_final_df %>% 
  filter(simulation == "SBP based sampling" & Imputed_Variable == "Abnormal sBP" & prop_missing == 0.9 & Criteria == "RMSE" & `Imputation Method` == "missing as normal") %>% 
  select(a_weights_vector_list, Calculation) %>% as.matrix() 


Frechet(aaa,ccc)
bbb

Frechet_distance_calculator <- function(the_df,the_simulation, the_imputed_variable,the_prop_missing, the_criteria, the_imputation_method){
  mice_matrix <- the_df %>% 
    filter(simulation == the_simulation & Imputed_Variable == the_imputed_variable & prop_missing == the_prop_missing & 
             Criteria == the_criteria & `Imputation Method` == "MICE") %>% 
    select(a_weights_vector_list, Calculation) %>% as.matrix() 
  
  second_matrix <- the_df %>% 
   filter(simulation == the_simulation & Imputed_Variable == the_imputed_variable & prop_missing == the_prop_missing & 
             Criteria == the_criteria & `Imputation Method` == the_imputation_method) %>% 
    select(a_weights_vector_list, Calculation) %>% drop_na() %>% as.matrix()
  
  Frechet(mice_matrix, second_matrix) %>% return()
#return(second_matrix)
  }

Frechet_distance_calculator_same_method(the_df = Mentation_final_df,the_simulation = "SBP based sampling", the_imputed_variable = "Abnormal SBP", 
                                        the_prop_missing = c(0.1,0.9), the_criteria = "RMSE", the_imputation_method = "missing discarded")


Mentation_Frechet_df <-
  lapply(Mentation_final_df$Imputed_Variable %>% unique, function(x){
    lapply(c("missing discarded", "missing as normal", "missing as abnormal"), function(y){
      lapply(c("RMSE", "PB"), function(z){
        lapply(c(0.1, 0.3, 0.5, 0.7, 0.9), function(a){
          the_result <-Frechet_distance_calculator(the_df = Mentation_final_df, the_simulation = "SBP based sampling", the_imputed_variable = x,
                                      the_prop_missing = a, the_criteria = z, the_imputation_method = y) 
          return(data.frame(Imputed_Variable = x, Imputation_method = y, Criteria = z, prop_missing = a, Frechet_distance = the_result))
        }) %>% do.call("rbind",.) 
      }) %>% do.call("rbind", .)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .)

rownames(Mentation_Frechet_df) <- NULL 
Mentation_Frechet_df<-remove_rownames(Mentation_Frechet_df)

Leukos_Frechet_df <-
  lapply(Leukos_final_df$Imputed_Variable %>% unique, function(x){
    lapply(c("missing discarded", "missing as normal", "missing as abnormal"), function(y){
      lapply(c("RMSE", "PB"), function(z){
        lapply(c(0.1, 0.3, 0.5, 0.7, 0.9), function(a){
          the_result <-Frechet_distance_calculator(the_df = Leukos_final_df, the_simulation = "Temperature based sampling", the_imputed_variable = x,
                                                   the_prop_missing = a, the_criteria = z, the_imputation_method = y) 
          return(data.frame(Imputed_Variable = x, Imputation_method = y, Criteria = z, prop_missing = a, Frechet_distance = the_result))
        }) %>% do.call("rbind",.) 
      }) %>% do.call("rbind", .)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .)
rownames(Leukos_Frechet_df) <- NULL 
Leukos_Frechet_df<-remove_rownames(Leukos_Frechet_df)

Mentation_Frechet_table <- Mentation_Frechet_df %>% filter(prop_missing != 0.9) %>%  group_by(Imputed_Variable,Imputation_method, Criteria) %>% 
  summarize(the_mean = mean(Frechet_distance), the_sd = sd(Frechet_distance), the_min = min(Frechet_distance), the_max = max(Frechet_distance), 
            the_median = median(Frechet_distance), the_IQR = IQR(Frechet_distance)) %>% filter(Criteria == "RMSE") %>% 
  mutate(results = str_c(the_mean %>% round(3), " (" , the_sd %>% round(3),") ", "[", the_min %>% round(3), "-", the_max %>% round(3), "]")) %>% 
  select(Imputed_Variable, Imputation_method, results) %>% spread(Imputation_method, results) %>% remove_rownames() %>% column_to_rownames("Imputed_Variable") %>% 
  select(`Discarded Curve` = `missing discarded`, `Normal Curve` = `missing as normal`, `Abnormal Curve` = `missing as abnormal` )
rownames(Mentation_Frechet_table) <-  rownames(Mentation_Frechet_table) %>% str_c(" mean (sd) [min-max]")

Leukos_Frechet_table <- Leukos_Frechet_df %>% filter(prop_missing != 0.9) %>% group_by(Imputed_Variable,Imputation_method, Criteria) %>% 
  summarize(the_mean = mean(Frechet_distance), the_sd = sd(Frechet_distance), the_min = min(Frechet_distance), the_max = max(Frechet_distance), 
            the_median = median(Frechet_distance), the_IQR = IQR(Frechet_distance)) %>% filter(Criteria == "RMSE") %>% 
  mutate(results = str_c(the_mean %>% round(3), " (" , the_sd %>% round(3),") ", "[", the_min %>% round(3), "-", the_max %>% round(3), "]")) %>% 
  select(Imputed_Variable, Imputation_method, results) %>% spread(Imputation_method, results) %>% remove_rownames() %>% column_to_rownames("Imputed_Variable") %>% 
  select(`Discarded Curve` = `missing discarded`, `Normal Curve` = `missing as normal`, `Abnormal Curve` = `missing as abnormal` )

rownames(Leukos_Frechet_table) <-  rownames(Leukos_Frechet_table) %>% str_c(" mean (sd) [min-max]")

rbind(Leukos_Frechet_table, Mentation_Frechet_table )%>% 
  htmlTable(rgroup = c("SIRS Coefficients","qSOFA Coefficients"),
            n.rgroup = c(5,nrow(.) - 5),
            cgroup = c("Ferchet Distance from MICE RMSE Curve"), n.cgroup = c(3))


Frechet_distance_calculator_same_method <- function(the_df,the_simulation, the_imputed_variable,the_prop_missing_vector, the_criteria, the_imputation_method){

  first_matrix <- the_df %>% 
    filter(simulation == the_simulation & Imputed_Variable == the_imputed_variable & prop_missing == the_prop_missing_vector[1] & 
             Criteria == the_criteria & `Imputation Method` == the_imputation_method) %>% 
    select(a_weights_vector_list, Calculation) %>% drop_na() %>% as.matrix()
  
 
  second_matrix <- the_df %>% 
    filter(simulation == the_simulation & Imputed_Variable == the_imputed_variable & prop_missing == the_prop_missing_vector[2] & 
             Criteria == the_criteria & `Imputation Method` == the_imputation_method) %>% 
    select(a_weights_vector_list, Calculation) %>% drop_na() %>% as.matrix()
  
  Frechet(first_matrix, second_matrix) %>% return()
  #return(second_matrix)
}

Mentation_Frechet_df_same_method <-
  lapply(Mentation_final_df$Imputed_Variable %>% unique, function(x){
    lapply(c("MICE","missing discarded", "missing as normal", "missing as abnormal"), function(y){
      lapply(c("RMSE", "PB"), function(z){
        lapply(list(c(0.1,0.3), c(0.1,0.5), c(0.1,0.7), c(0.1,0.9), c(0.3,0.5), c(0.3,0.7), c(0.3, 0.9), c(0.5, 0.7), c(0.5,0.9), c(0.7,0.9)), function(a){
          the_result <-Frechet_distance_calculator_same_method(the_df = Mentation_final_df, the_simulation = "SBP based sampling", the_imputed_variable = x,
                                                   the_prop_missing = a, the_criteria = z, the_imputation_method = y) 
          return(data.frame(Imputed_Variable = x, Imputation_method = y, Criteria = z, prop_missing = a[2], Frechet_distance = the_result))
        }) %>% do.call("rbind",.) 
      }) %>% do.call("rbind", .)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .)

rownames(Mentation_Frechet_df_same_method) <- NULL

Leukos_Frechet_df_same_method <-
  lapply(Leukos_final_df$Imputed_Variable %>% unique, function(x){
    lapply(c("MICE","missing discarded", "missing as normal", "missing as abnormal"), function(y){
      lapply(c("RMSE", "PB"), function(z){
        lapply(list(c(0.1,0.3), c(0.1,0.5), c(0.1,0.7), c(0.1,0.9), c(0.3,0.5), c(0.3,0.7), c(0.3, 0.9), c(0.5, 0.7), c(0.5,0.9), c(0.7,0.9)), function(a){
          the_result <-Frechet_distance_calculator_same_method(the_df = Leukos_final_df, the_simulation = "Temperature based sampling", the_imputed_variable = x,
                                                               the_prop_missing = a, the_criteria = z, the_imputation_method = y) 
          return(data.frame(Imputed_Variable = x, Imputation_method = y, Criteria = z, prop_missing = a[2], Frechet_distance = the_result))
        }) %>% do.call("rbind",.) 
      }) %>% do.call("rbind", .)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .)

rownames(Leukos_Frechet_df_same_method) <- NULL

Mentation_Frechet_same_method_table <- Mentation_Frechet_df_same_method %>%  group_by(Imputed_Variable,Imputation_method, Criteria) %>% 
  summarize(the_mean = mean(Frechet_distance), the_sd = sd(Frechet_distance), the_min = min(Frechet_distance), the_max = max(Frechet_distance), 
            the_median = median(Frechet_distance), the_IQR = IQR(Frechet_distance)) %>% filter(Criteria == "RMSE") %>% 
  mutate(results = str_c(the_mean %>% round(2), " (" , the_sd %>% round(3),") ", "[", the_min %>% round(2), "-", the_max %>% round(2), "]")) %>% 
  select(Imputed_Variable, Imputation_method, results) %>% spread(Imputation_method, results) %>% remove_rownames() %>% column_to_rownames("Imputed_Variable") %>% 
  select(`MICE Curves` = MICE,`Discarded Curves` = `missing discarded`, `Normal Curves` = `missing as normal`, `Abnormal Curves` = `missing as abnormal` )
rownames(Mentation_Frechet_same_method_table) <-  rownames(Mentation_Frechet_table) %>% str_c(" mean (sd) [min-max]")

Leukos_Frechet_same_method_table <- Leukos_Frechet_df_same_method %>% group_by(Imputed_Variable,Imputation_method, Criteria) %>% 
  summarize(the_mean = mean(Frechet_distance), the_sd = sd(Frechet_distance), the_min = min(Frechet_distance), the_max = max(Frechet_distance), 
            the_median = median(Frechet_distance), the_IQR = IQR(Frechet_distance)) %>% filter(Criteria == "RMSE") %>% 
  mutate(results = str_c(the_mean %>% round(3), " (" , the_sd %>% round(2),") ", "[", the_min %>% round(2), "-", the_max %>% round(2), "]")) %>% 
  select(Imputed_Variable, Imputation_method, results) %>% spread(Imputation_method, results) %>% remove_rownames() %>% column_to_rownames("Imputed_Variable") %>% 
  select(`MICE Curves` = MICE,`Discarded Curves` = `missing discarded`, `Normal Curves` = `missing as normal`, `Abnormal Curves` = `missing as abnormal`)

rownames(Leukos_Frechet_same_method_table) <-  rownames(Leukos_Frechet_table) %>% str_c(" mean (sd) [min-max]")

rbind(Leukos_Frechet_same_method_table, Mentation_Frechet_same_method_table )%>% 
  htmlTable(rgroup = c("SIRS Coefficients","qSOFA Coefficients"),
            n.rgroup = c(5,nrow(.) - 5),
            cgroup = c("Ferchet Distance Between RMSE Curves"), n.cgroup = c(4))



