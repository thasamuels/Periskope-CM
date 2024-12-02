# Cryptococcal meningitis prognostic model - Periskope-CM
# 12: Save Outputs
# Description: Saves relevant outputs for Rmd file. 
# Started: 23/11/23
# Author: Tom Samuels, University College London

library(tidyverse)
library(ggpubr)
library(janitor)
library(rms)
library(naniar)
library(pROC)
library(flextable)
library(gtsummary)
library(Amelia)
library(metafor)
library(rmda)
library(mice)
library(viridis)
library(readr)
library(gt)
library(Hmisc)
library(corrplot)
library(patchwork)
library(webshot2)
library(finalfit)
library(survminer)
library(interactionRCS)
library(gtools)

# Source from relevant R code to create objects for saving (takes a while to run). 
source("11 - Treatment Effect Analysis.R")
source("02 - Descriptive Analysis.R")
source("04 - MI variable selection.R")
source("05 - IECV.R")
source("06 - Held-out Validation.R")
source('07 - DCA and Benchmarking.R')
source("08 - Sensitivity Analyses.R")
source("09 - XGBoost Machine Learning.R")
source("10 - Ten-week Regression Model.R")

# Pull through objects from XGBoost analysis
calibration_val_x <- readRDS("XGBoost data/calibration_val_x.rds")
marginals_plot_x <- readRDS("XGBoost data/marginals_plot.rds")
validation_summary_pooled_x <- readRDS("XGBoost data/validation_summary_pooled_x.rds")
prediction_density_x_val <- readRDS("XGBoost data/prediction_density_XGB.rds")
iecv_master_pooled_x <- readRDS("XGBoost data/iecv_master_pooled_x.rds")

# Create validation results dataframe for XGBoost using summary
validation_results_x <- validation_summary_pooled_x %>%
  mutate(
    Model = 'XGBoost',
    AUROC = paste0(round(auc, 2), " (", round(auc.lower, 2), " - ", round(auc.upper, 2), ")"),
    Slope = paste0(round(slope, 2), " (", round(slope.lower, 2), " - ", round(slope.upper, 2), ")"),
    CITL = paste0(round(citl, 2), " (", round(citl.lower, 2), " - ", round(citl.upper, 2), ")")
  ) %>%
  select(Model, AUROC, Slope, CITL)

# Save objects for Rmd file. 
save(
  # Summary Tables
  tbl_summary_1, #with missing (n)
  tbl_summary_country, #distribution of variables by country
  tbl_summary_devval, #distribution of variables by development vs validation cohort
  
  # Main models
  tbl_var_sel_b_detail, #Variables selected into MI Basic Model by MI dataset
  basic_plot_assoc_mi,  #Plot associations for the MI Basic Model
  basic_regress_table, # Table of coefficients for Basic MI Model
  basic_mi_model, #Basic Model Object
  tbl_var_sel_r_detail, #Variables selected into MI Research Model by MI dataset
  table_var_sel_complete, #Combined variable selection table
  research_plot_assoc_mi, #Plot associations for the MI Research Model
  research_regress_table, # Table of log(OR) for Research MI Model
  research_mi_model, #Model object
  predictions_density, #Density Plot of both models predictions in all MI data
  
  # IECV
  iecv_master_pooled_b, #IECV output data for Basic Model
  iecv_master_pooled_r, #IECV output data for Research Model
  model.cstat.b, model.CITL.b, model.CS.b, #Values for IECV forest plot Basic Model
  model.cstat.r, model.CITL.r, model.CS.r, #Values for IECV Forest Research Plot
  calibration_iecv_b, calibration_iecv_r, # IECV Calibration plots
  calibration_iecv_recal_b, calibration_iecv_recal_r, #IECV Recalibration plots
  iecv_prediction_density, #Density Plot of Prediction values by MI dataset for Research model
  
  #Validation
  validation_summary_pooled_b, validation_summary_pooled_r, #Summary validation data for B and R Model
  validation_table, combined_valdata, # Summary Table of Validation stats for R and B Models
  calibration_val_b, calibration_val_r, #Pooled Validation Calibration Plots for R and B Model
  calibration_val_each_mi_b, calibration_val_each_mi, #Cal plot by MI dataset for R and B
  calibration_b_recal, calibration_r_recal, #Recalibration Validation Plots for R and B
  iecv_prediction_density_val, combined_devval_density, #Validation and Combined Density Predictions by MI
  prediction_density_research_val, prediction_density_basic_val, #Density predictions for research and basic model, pooled MI sets
  dca_plot_val_combined, # Decision curve analysis plot
  
  #Treatment effect results based on main model
  risk_treat_effect_plot_combined,
  ambition_treat_effect_summary, ambition_treat_effect_summary_gt,
  acta_treat_effect_summary, acta_treat_effect_summary_gt,
  basic_group_boxplot, research_group_boxplot,
  tercile_thresholds_b, tercile_thresholds_r,
  box_dens_arranged_nonrotated,
  
  # Sensitivity analyses
  iecv_singlefactor_auc_table, single_predictor_table, # single factor table AUC (pooled and IECV)
  validation_table_10w, calibration_val_b_10w, calibration_val_r_10w, #10 week endpoint for 2 week model
  
  #Ten week model
  tbl_var_sel_b_detail_10, basic_plot_assoc_mi_10, iecv_master_pooled_b10, model.cstat.b.10, 
  model.CITL.b.10, model.CS.b.10, calibration_iecv_b10,
  tbl_var_sel_r_detail_10, research_plot_assoc_mi_10, iecv_master_pooled_r10, model.cstat.r.10, 
  table_var_sel_complete_10week,
  model.CITL.r.10, model.CS.r.10, calibration_iecv_r10, calibration_iecv_recal_r10,
  calibration_val_r_10week, calibration_val_b_10week, 
  research_regress_table_10w, basic_regress_table_10w,
  
  # Sample sizes and other values
  samplesize_total, samplesize_dev, samplesize_val, samplesize_acta, samplesize_ambition,
  median_research_prediction, median_basic_prediction, 
  acta_death_2week, acta_death_10week, ambition_death_2week, ambition_death_10week,
  total_death_2week, total_death_10week, dev_death_2week, dev_death_10week, val_death_2week, val_death_10week,
  flowchart_data, median_r_pred_death, median_r_pred_alive, median_b_pred_death, median_b_pred_alive,
  research_pred_stats, basic_pred_stats, #stats for all-mi predictions for both main models, in df format
  basic_cstat_i2, basic_citl_i2, basic_slope_i2, # measures of heterogeneity for IECV
  research_cstat_i2, research_citl_i2, research_slope_i2,
  basic10w_cstat_i2, basic10w_citl_i2, basic10w_slope_i2,
  research10w_cstat_i2, research10w_citl_i2, research10w_slope_i2,
  x_cstat_i2, x_citl_i2, x_slope_i2,
  
  # Results dataframes for combine
  validation_results_main, #validation_results_treat, 
  validation_results_mainmodel_10w,
  validation_results_10weekmodel,
  
  # Delong analysis and stratified validation
  delong_table, validation_strat,
  
  # XGBoost model
  validation_summary_pooled_x, marginals_plot_x, calibration_val_x, validation_results_x, 
  prediction_density_x_val, iecv_master_pooled_x, xgb_variableimport,
  
  # Comparative analysis with Zhao et al model
  validation_results_zhao, calibration_val_z, 
  
  file = "CM_table_and_plot_objects_231123_v2.rda"
)

# Get Model formulas

(pred_logit_b <- Function(basic_mi_model))
(pred_logit <- Function(research_mi_model))


# Final Model Formulas

Periskope_CM_basic <- function(ecog, gcs_bin, haem_gl, neut, treatment){
  inv.logit(-1.9538382+0.55680611*(gcs_bin=="11-14")+0.99380898*(gcs_bin=="<=10")+1.0270395*(ecog=="Restricted activity")+
              1.7235581*(ecog=="Ambulatory")+1.7816894*(ecog=="Limited self-care")+2.7167113*(ecog=="Bedbound")+1.0192472*(treatment=="AmBd_Flu_1w")+
              0.16194237*(treatment=="AmBd_5FC_2w")+0.5204576*(treatment=="AmBd_Flu_2w")+0.12240689*(treatment=="Flu_5FC_PO")+
              0.3132614* neut-0.0097500838*pmax(neut-1.12,0)^3+0.014065387*pmax(neut-2.5,0)^3-0.0043153033*pmax(neut-5.618,0)^3-0.029004943* 
              haem_gl+7.3953556e-06*pmax(haem_gl-83.000002,0)^3-1.4790711e-05*pmax(haem_gl-110,0)^3+7.3953551e-06*pmax(haem_gl-137,0)^3)
}

Periskope_CM_research <- function(ecog, gcs_bin, haem_gl, neut, treatment, csf_OP, csf_qculture_log){
  inv.logit(-1.7666866+0.58108236*(gcs_bin=="11-14")+1.2728151*(gcs_bin=="<=10")+0.93532768*(ecog=="Restricted activity")+1.6564446*(ecog=="Ambulatory")+
              1.6210341*(ecog=="Limited self-care")+2.549854*(ecog=="Bedbound")+0.88776442*(treatment=="AmBd_Flu_1w")+0.10466669*(treatment=="AmBd_5FC_2w")+
              0.40279503*(treatment=="AmBd_Flu_2w")+0.06123988*(treatment=="Flu_5FC_PO")+0.25507429* 
              neut-0.0050287923*pmax(neut-1.12,0)^3+0.0072544926*pmax(neut-2.5,0)^3-0.0022257003*pmax(neut-5.618,0)^3-0.033843364* 
              haem_gl+7.4246174e-06*pmax(haem_gl-83.000002,0)^3-1.4849234e-05*pmax(haem_gl-110,0)^3+7.4246169e-06*pmax(haem_gl-137,0)^3+0.0039830694* 
              csf_OP+7.3805673e-06*pmax(csf_OP-8,0)^3-1.1873087e-05*pmax(csf_OP-22,0)^3+4.4925192e-06*pmax(csf_OP-45,0)^3-0.21419185* 
              csf_qculture_log+0.0154645*pmax(csf_qculture_log-0.69897,0)^3-0.066804923*pmax(csf_qculture_log-4.9799974,0)^3+0.051340423*pmax(csf_qculture_log-6.2695067,0)^3)
}

