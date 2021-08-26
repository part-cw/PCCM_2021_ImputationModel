# PCCM_2021_ImputationModel
Title: Prediction model performance with different imputation strategies – a simulation study using a North American intensive care unit registry
Authors: Jonathan Steif, Rollin Brant, Rama Sreepada, Nicholas West, Srinivas Murthy, Matthias Görges

1. helper_functions.R, add_PALS_criteria.R, add_qPELOD2_criteria.R, add_sirs_criteria.R: these files are used to pre-process the initial dataset to create new features. These files are called from other files, and need not be run explicitly by the user.
2. Leukos_Mentation_Coefficients.R : this file is used to build pSIRS and qSOFA models using 150-run, and compute the coefficients for missing as normal, missing discarded, and MICE imputation strategies.
3. Coefficient_Plots.R: this file is used to pool the coefficients from 150 runs, and generate percent bias, RMSE plots, and Frechet distance table, and heatmaps provided with the article.
4. AUC_WithCIs_Code_Leukos.R, AUC_WithCIs_Code_Mentation.R : These files are used to compute AUC values using 150-run, 5-fold cross-validation on pSIRS and qSOFA models for missing as normal, missing discarded, and MICE imputation strategies.
5. AUC_Result_Plots.R : This file is used to pool AUC results from 150-run, and generate AUC plots for pSIRS and qSOFA models. 
