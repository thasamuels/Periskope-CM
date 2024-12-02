# Cryptococcal meningitis prognostic model - Periskope-CM
# 11: Treatment Effect Analysis using Primary model
# Description: Treatment effect analysis, by category and modeled continuously. 
# Started: 5/1/24
# Author: Tom Samuels and Rishi Gupta, University College London

library(tidyverse)
library(ggpubr)
library(janitor)
library(rms)
library(naniar)
library(pROC)
library(flextable)
library(gtsummary)
library(Amelia)
library(rmda)
library(mice)
library(viridis)
library(Hmisc)
library(corrplot)
library(gt)
library(cowplot)
library(patchwork)
library(ggpubr)

#### Note, if run sequentially after other files in CM analysis, the treatment effect plots will error due to a datadist issue. 
# There is no easy fix to this, so best to run this file first. 

#Read in data
all_mi <- readRDS("complete_mi_aregimpute_131123.rds")

# Clear labels
all_mi <- data.frame(all_mi)

# Read in final models
basic_model_2w <- readRDS("basic_mi_model_251023.rds")
research_model_2w <- readRDS("research_mi_model_251023.rds")
basic_model_10w <- readRDS("basic_mi_model_10week_161123.rds")
research_model_10w <- readRDS("research_mi_model_10week_161123.rds")

# Re-train the model for the ACTA data (without treatment as a component variable)

basic_model_notreat <- readRDS('basicmodel_2w_notreat.rds')
research_model_notreat <- readRDS('researchmodel_2w_notreat.rds')

# Set treatment variable to baseline, and transfer real treatment allocation
all_mi <- all_mi %>%
  mutate(treatment_real = treatment) %>%
  mutate(treatment = '1wk AmBd+5FC/Ambisome')

# Create treatment model predictions for MI data
all_mi$treat_pred_b <- predict(basic_model_2w, type = 'fitted', newdata = all_mi) # change here if want to use 10w models
all_mi$treat_pred_r <- predict(research_model_2w, type = 'fitted', newdata = all_mi)

# Add risk strata
tercile_thresholds_b <- quantile(all_mi$treat_pred_b[all_mi$.imp>0], probs = c(1/3, 2/3), na.rm = T)
tercile_thresholds_r <- quantile(all_mi$treat_pred_r[all_mi$.imp>0], probs = c(1/3, 2/3), na.rm = T)

all_mi <- all_mi %>%
  mutate(
    basic_risk_group = cut(treat_pred_b,
                           breaks = c(-Inf, tercile_thresholds_b, Inf), # Ensure all values are covered
                           labels = c("Low Risk", "Medium Risk", "High Risk"),
                           include.lowest = TRUE),
    research_risk_group = cut(treat_pred_r,
                              breaks = c(-Inf, tercile_thresholds_r, Inf), # Ensure all values are covered
                              labels = c("Low Risk", "Medium Risk", "High Risk"),
                              include.lowest = TRUE),
    death_10week = (as.numeric(death_10week) - 1)
  )

# tabyl(all_mi$basic_risk_group)
# tabyl(all_mi$research_risk_group)


#################### Create boxplots and density plots of risk distribution #########################



basic_group_boxplot <- ggplot(data = all_mi[all_mi$.imp>0,],
                              aes(x = basic_risk_group, 
                                  y = treat_pred_b,
                                  fill = factor(basic_risk_group)
                              )) +
  coord_cartesian(ylim = c(0,1)) +
  geom_boxplot(show.legend = F) +
  labs(y = NULL,
       x = 'Tercile Group') +
  scale_fill_manual(values = c("Low Risk" = "#75C141", "Medium Risk" = "#FFA500", 
                               "High Risk" = "#FF6347")) +
  geom_hline(yintercept = c(tercile_thresholds_b[1], tercile_thresholds_b[2]), 
             linetype = "dashed", color = "black", linewidth = 0.3) +
  theme_pubr() +
  theme(axis.title.x = element_text(margin = margin(t = 15)))

basic_group_boxplot

research_group_boxplot <- ggplot(data = all_mi[all_mi$.imp>0,],
                                 aes(x = research_risk_group, 
                                     y = treat_pred_r,
                                     fill = factor(research_risk_group)
                                 )) +
  coord_cartesian(ylim = c(0,1)) +
  geom_boxplot(show.legend = F) +
  labs(y = NULL,
       x = 'Tercile Group') +
  scale_fill_manual(values = c("Low Risk" = "#75C141", "Medium Risk" = "#FFA500", 
                               "High Risk" = "#FF6347")) +
  geom_hline(yintercept = c(tercile_thresholds_r[1], tercile_thresholds_r[2]), 
             linetype = "dashed", color = "black", linewidth = 0.3) +
  theme_pubr() +
  theme(axis.title.x = element_text(margin = margin(t = 15)))

research_group_boxplot


# Extract summary statistics for distributions
# For the mode, need to round up predictions or wont be able to identify correctly as too granular

all_mi <- all_mi %>% 
  mutate(
    pred_r_round = round(treat_pred_r, 3),
    pred_b_round = round(treat_pred_b, 3)
  )

research_pred_stats <- all_mi %>% 
  filter(.imp > 0) %>% 
  summarise(
    mean = round(mean(treat_pred_r), 2),
    median = round(median(treat_pred_r), 2),
    lquartile = round(quantile(treat_pred_r, probs = 0.25), 2),
    uquartile = round(quantile(treat_pred_r, probs = 0.75), 2),
    mode = as.numeric(names(sort(table(pred_r_round), decreasing = T))[1])
  )

basic_pred_stats <- all_mi %>% 
  filter(.imp > 0) %>% 
  summarise(
    mean = round(mean(treat_pred_b), 2),
    median = round(median(treat_pred_b), 2),
    lquartile = round(quantile(treat_pred_b, probs = 0.25), 2),
    uquartile = round(quantile(treat_pred_b, probs = 0.75), 2),
    mode = as.numeric(names(sort(table(pred_b_round), decreasing = T))[1])
  )

### Create density distributions for predictions

prediction_density_research <- ggplot(data = all_mi[all_mi$.imp>0,], 
                                      aes(x = treat_pred_r)) +
  geom_density(
    alpha = 0.6, linewidth = 0.5, show.legend = F, fill = '#587D97') +
  labs(y = "Prediction Density",
       x = 'Predicted Mortality Risk') +
  geom_vline(xintercept = c(tercile_thresholds_r[1], tercile_thresholds_r[2]), 
             linetype = "dashed", color = "black", linewidth = 0.3) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_pubr()

prediction_density_research

prediction_density_basic <- ggplot(data = all_mi[all_mi$.imp>0,], 
                                   aes(x = treat_pred_b)) +
  geom_density(
    alpha = 0.6, linewidth = 0.5, 
    fill = '#A4C0D9',
    show.legend = F) +
  labs(y = "Prediction Density",
       x = 'Predicted Mortality Risk') +
  geom_vline(xintercept = c(tercile_thresholds_b[1], tercile_thresholds_b[2]), 
             linetype = "dashed", color = "black", linewidth = 0.3) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_pubr()

prediction_density_basic

# Non-rotated version
box_dens_r_nonrotated <- plot_grid(research_group_boxplot, prediction_density_research, align = 'hv', rel_widths = c(1,1))
box_dens_b_nonrotated <- plot_grid(basic_group_boxplot, prediction_density_basic, align = 'hv', rel_widths = c(1,1))

box_dens_arranged_nonrotated <- ggarrange(
  basic_group_boxplot, prediction_density_basic, research_group_boxplot, prediction_density_research, 
  nrow = 2, ncol = 2,
  labels = c("A", "B", 'C', 'D'),
  label.x = c(0.15, 0, 0.15, 0)
)
box_dens_arranged_nonrotated <- annotate_figure(box_dens_arranged_nonrotated, left = "Predicted Mortality Risk")
box_dens_arranged_nonrotated


###################### Treatment effect table ############################

# Split Ambition data

ambition  <- all_mi %>% filter(trial == "Ambition") %>%
  mutate(
    treatment = ifelse(trial_arm == 'Ambisome', "Ambisome", "AmBd_1wk") %>%
      factor(levels = c("AmBd_1wk", 'Ambisome')))

# Split ACTA data

acta <- all_mi %>% filter(trial != 'Ambition')

acta_2arms <- acta %>%
  filter(treatment_real == '1wk AmBd+5FC/Ambisome' | treatment_real == 'FLU+5FC PO') %>%
  mutate(
    treatment = ifelse(treatment_real == '1wk AmBd+5FC/Ambisome', "AmBd_1wk", "Oral_regimen"),
    treatment = factor(treatment, levels = c("AmBd_1wk", "Oral_regimen")))


# Function to get risk-stratified treatment effects

get_risk_strata_treat_effects <- function(data, model, risk_strata){
  
  # Absolute mortality summary
  
  data1 <- data
  data1$risk_group <- data1[[risk_strata]]
  
  abs_mort <- data1 %>%
    filter(.imp==0) %>%
    group_by(risk_group, treatment) %>%
    dplyr::summarize(
      deaths = sum(death_10week),
      total = n(),
      pct_death = round(100 * deaths / total, 1)) %>%
    mutate(
      summary = paste0(
        deaths, "/", total, " (", pct_death, "%)")) %>%
    mutate(model = model) %>%
    select(model, risk_group, treatment, summary) %>%
    pivot_wider(
      names_from = treatment, values_from = ("summary"))
  
  # Create a MIDS object
  
  data1_mids <- data1 %>%
    select(.imp, study_id,
           age, sex, gcs_bin, ecog,
           neut, haem_gl,
           treatment,
           csf_OP, csf_qculture_log,
           time_days_numeric,
           risk_group, 
           death_2week, death_10week) %>%
    as.mids(.id = "study_id")
  
  # Linear model for risk difference
  
  glm_contrast_treatmodel <- fit.mult.impute(death_10week ~ risk_group * treatment,
                                             Glm,
                                             #family = gaussian(link = "identity"),
                                             data1_mids)
  
  contrast_riskdiff_treatment <- contrast(glm_contrast_treatmodel,
                                          list(treatment = levels(data1$treatment)[2], risk_group = c('Low Risk', 'Medium Risk', 'High Risk')),
                                          list(treatment = levels(data1$treatment)[1], risk_group = c('Low Risk', 'Medium Risk', 'High Risk')),
                                          type = 'joint')
  
  contrast_treatment_risk_diff <- data.frame(matrix(nrow = 3, ncol = 5))
  contrast_treatment_risk_diff$risk_group <- contrast_riskdiff_treatment$risk_group
  contrast_treatment_risk_diff$estimate <- contrast_riskdiff_treatment$Contrast
  contrast_treatment_risk_diff$lower_ci <- contrast_riskdiff_treatment$Lower
  contrast_treatment_risk_diff$upper_ci <- contrast_riskdiff_treatment$Upper
  contrast_treatment_risk_diff$p.value.mort <- contrast_riskdiff_treatment$Pvalue
  
  mort_diff_formerge <- contrast_treatment_risk_diff %>%
    select(risk_group, estimate, lower_ci, upper_ci, p.value.mort) %>%
    mutate_at(vars(c(estimate, lower_ci, upper_ci)), ~ (.) * 100) %>%
    mutate(
      mort_diff = paste0(round(estimate, 1), "% (", round(lower_ci, 1), "-", round(upper_ci, 1), ")"),
    ) %>%
    select(risk_group, mort_diff, p.value.mort) %>%
    mutate(
      p.value.mort = case_when(
        p.value.mort >= 0.1 ~ round(p.value.mort, 1),
        p.value.mort < 0.1 & p.value.mort >= 0.01 ~ round(p.value.mort, 2),
        p.value.mort < 0.01 & p.value.mort >= 0.001 ~ round(p.value.mort, 3),
        p.value.mort <0.001 ~ p.value.mort
      )
    ) %>%
    mutate(
      p.value.mort = ifelse(p.value.mort <0.001, "<0.001", p.value.mort)
    )
  
  mort_diff_formerge
  
  # Cox model for hazard ratios
  
  cox_treat_risk_group <- fit.mult.impute(Surv(time = time_days_numeric, death_10week) ~ risk_group * treatment,
                                          cph,
                                          data1_mids)
  
  contrast_cox_treatment <- contrast(cox_treat_risk_group,
                                     list(treatment = levels(data1$treatment)[2], risk_group = c('Low Risk', 'Medium Risk', 'High Risk')),
                                     list(treatment = levels(data1$treatment)[1], risk_group = c('Low Risk', 'Medium Risk', 'High Risk')),
                                     type = 'joint')
  
  contrast_treatment_hr <- data.frame(matrix(nrow = 3, ncol = 5))
  contrast_treatment_hr$risk_group <- contrast_cox_treatment$risk_group
  contrast_treatment_hr$estimate <- contrast_cox_treatment$Contrast
  contrast_treatment_hr$lower_ci <- contrast_cox_treatment$Lower
  contrast_treatment_hr$upper_ci <- contrast_cox_treatment$Upper
  contrast_treatment_hr$p_value <- contrast_cox_treatment$Pvalue
  
  treat_group_hr_formerge <- contrast_treatment_hr %>%
    select(risk_group, estimate, lower_ci, upper_ci, p_value) %>%
    rename(p.value.hr = p_value) %>%
    mutate_at(vars(c(estimate, lower_ci, upper_ci)), ~ exp(.)) %>%
    mutate(
      HR = paste0(round(estimate, 2), " (", round(lower_ci, 2), "-", round(upper_ci, 2), ")"),
    ) %>%
    select(risk_group, HR, p.value.hr) %>%
    mutate(
      p.value.hr = case_when(
        p.value.hr >= 0.1 ~ round(p.value.hr, 1),
        p.value.hr < 0.1 & p.value.hr >= 0.01 ~ round(p.value.hr, 2),
        p.value.hr < 0.01 & p.value.hr >= 0.001 ~ round(p.value.hr, 3),
        p.value.hr <0.001 ~ p.value.hr
      )
    ) %>%
    mutate(
      p.value.hr = ifelse(p.value.hr <0.001, "<0.001", p.value.hr)
    )
  
  treat_group_hr_formerge
  
  # Join tables
  
  treatment_effect_table <- abs_mort %>% filter(!is.na(risk_group)) %>%
    left_join(mort_diff_formerge) %>%
    left_join(treat_group_hr_formerge)
  
  return(treatment_effect_table)
  
}

## Run for Ambition

ambition_basic <- get_risk_strata_treat_effects(data=ambition, risk_strata = "basic_risk_group", model = "Basic")
ambition_research <- get_risk_strata_treat_effects(data=ambition, risk_strata = "research_risk_group", model = "Research")
ambition_treat_effect_summary <- bind_rows(ambition_basic, ambition_research)
ambition_treat_effect_summary

## Run for ACTA

acta_basic <- get_risk_strata_treat_effects(data=acta_2arms, risk_strata = "basic_risk_group", model = "Basic")
acta_research <- get_risk_strata_treat_effects(data=acta_2arms, risk_strata = "research_risk_group", model = "Research")
acta_treat_effect_summary <- bind_rows(acta_basic, acta_research)
acta_treat_effect_summary

## Tidy up tables

### ACTA

acta_treat_effect_summary_gt <- acta_treat_effect_summary %>% 
  arrange(model) %>%
  gt(rowname_col = 'risk_group',
     groupname_col = 'model') %>%
  tab_spanner(columns = c(3:4), label = "Deaths") %>%
  tab_spanner(columns = c(5:6), label = "Mortality Difference") %>%
  tab_spanner(columns = c(7:8), label = "Hazard Ratio") %>%
  cols_label(AmBd_1wk = "SOC",
             Oral_regimen = "Oral regimen",
             mort_diff = "Oral v SOC",
             p.value.mort = "p value",
             HR = "Oral v SOC",
             p.value.hr = "p value") %>%
  tab_style(style = cell_text(weight = "bold", align = 'center'), locations = cells_column_labels()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_spanners()) %>%
  tab_style(style = cell_text(weight = 'bold'), locations = cells_row_groups()) %>%
  tab_footnote(footnote = paste0("Risk groups were defined by terciles of predicted risk in the pooled dataset. Thresholds of ", round((tercile_thresholds_b[1] * 100), 1), "% and ", 
                                 round((tercile_thresholds_b[2] * 100), 1), "% delineated Low-Medium and Medium-High risk respectively in the basic treatment model. 
               Thresholds of ", round((tercile_thresholds_r[1] * 100), 1), "% and ", round((tercile_thresholds_r[2] * 100), 1), 
               "% delineated Low-Medium and Medium-High risk respectively in the research treatment model."),
               locations = cells_row_groups()) %>%
  tab_footnote(footnote = "SOC = Standard of Care", locations = cells_column_labels(columns = c(3,5,7))) %>%
  tab_footnote(footnote = "Brackets indicate proportions", 
               locations = cells_column_labels(columns = c(3:4))) %>%
  tab_footnote(footnote = "Brackets indicate 95% confidence intervals", 
               locations = cells_column_labels(columns = c(5,7))) 

acta_treat_effect_summary_gt

### Ambition

ambition_treat_effect_summary_gt <- ambition_treat_effect_summary %>% 
  arrange(model) %>%
  gt(rowname_col = 'risk_group',
     groupname_col = 'model') %>%
  tab_spanner(columns = c(3:4), label = "Deaths") %>%
  tab_spanner(columns = c(5:6), label = "Mortality Difference") %>%
  tab_spanner(columns = c(7:8), label = "Hazard Ratio") %>%
  tab_stubhead(label = 'Ambition') %>%
  cols_label(AmBd_1wk = "SOC",
             Ambisome = "Ambisome",
             mort_diff = "Ambisome v SOC",
             p.value.mort = "p value",
             HR = "Ambisome v SOC",
             p.value.hr = "p value") %>%
  tab_style(style = cell_text(weight = "bold", align = 'center'), locations = cells_column_labels()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_spanners()) %>%
  tab_style(style = cell_text(weight = 'bold'), locations = cells_row_groups()) %>%
  tab_style(style = cell_text(weight = 'bold', align = 'center', v_align = 'middle'), locations = cells_stubhead()) %>%
  tab_footnote(footnote = paste0("Risk groups were defined by terciles of predicted risk in the pooled dataset. Thresholds of ", round((tercile_thresholds_b[1] * 100), 1), "% and ", 
                                 round((tercile_thresholds_b[2] * 100), 1), "% delineated Low-Medium and Medium-High risk respectively in the basic treatment model. 
               Thresholds of ", round((tercile_thresholds_r[1] * 100), 1), "% and ", round((tercile_thresholds_r[2] * 100), 1), 
               "% delineated Low-Medium and Medium-High risk respectively in the research treatment model."),
               locations = cells_row_groups()) %>%
  tab_footnote(footnote = "SOC = Standard of Care", locations = cells_column_labels(columns = c(3,5,7))) %>%
  tab_footnote(footnote = "Brackets indicate proportions", 
               locations = cells_column_labels(columns = c(3:4))) %>%
  tab_footnote(footnote = "Brackets indicate 95% confidence intervals", 
               locations = cells_column_labels(columns = c(5,7))) 

ambition_treat_effect_summary_gt



####################### Treatment Effect Plots modelled with continuous risk #################

# Hazard ratios----

## Function

get_risk_splines_treat_effects <- function(data, model, risk_var, trial_name){
  
  data1 <- data
  data1$pred_risk <- data1[[risk_var]]
  
  data1 <- data1 %>%
    mutate(
      treatment = ifelse(treatment == "AmBd_1wk", 0, 1)
    )
  
  # Create a MIDS object
  data1_mids <- data1 %>%
    select(.imp, study_id,
           age, sex, gcs_bin, ecog,
           neut, haem_gl,
           treatment,
           csf_OP, csf_qculture_log,
           time_days_numeric,
           pred_risk, 
           death_2week, death_10week) %>%
    as.mids(.id = "study_id")
  
  # Run the fit.mult.impute model
  cox_treat_effect_spl <- fit.mult.impute(Surv(time = time_days_numeric, death_10week) ~ treatment * rcs(pred_risk, 3),
                                          cph,
                                          data1_mids)
  
  
  # Define a sequence of predicted risk values for which to calculate hazard ratios
  covariate_values <- seq(0, 0.45, by=0.0005)
  
  # Use the contrast() function to calculate hazard ratios
  hr_list <- lapply(covariate_values, function(x) {
    
    contrasts <- contrast(cox_treat_effect_spl, 
                          list(treatment = 1, pred_risk = x),
                          list(treatment = 0, pred_risk = x))
    #c(HR = 100*contrasts$Contrast, Lower = 100*contrasts$Lower, Upper = 100*contrasts$Upper) # to plot risk differences
    c(HR = exp(contrasts$Contrast), Lower = exp(contrasts$Lower), Upper = exp(contrasts$Upper))
  })
  
  # Convert the list to a data frame
  hr_df <- data.frame(do.call(rbind, hr_list))
  names(hr_df) <- c("HR", "Lower", "Upper")
  hr_df$pred_risk <- covariate_values
  hr_df$model <- model
  hr_df$trial <- trial_name
  
  return(hr_df)
}

## Run for each model and dataset and combine

ambition_basic_spl <- get_risk_splines_treat_effects(data=ambition, risk_var = "treat_pred_b", model = "Basic", trial_name = "Ambition (Ambition regimen vs 1-wk AmBd+5FC)")
ambition_research_spl <- get_risk_splines_treat_effects(data=ambition, risk_var = "treat_pred_r", model = "Research", trial_name = "Ambition (Ambition regimen vs 1-wk AmBd+5FC)")
acta_basic_spl <- get_risk_splines_treat_effects(data=acta_2arms, risk_var = "treat_pred_b", model = "Basic", trial_name = "ACTA (Oral regimen vs 1-wk AmBd+5FC)")
acta_research_spl <- get_risk_splines_treat_effects(data=acta_2arms, risk_var = "treat_pred_r", model = "Research", trial_name = "ACTA (Oral regimen vs 1-wk AmBd+5FC)")
risk_hr_spl_all <- bind_rows(ambition_basic_spl, ambition_research_spl, acta_basic_spl, acta_research_spl)

## Plot

risk_treat_effect_plot_acta_basic <- 
  ggplot(data=risk_hr_spl_all[risk_hr_spl_all$trial == 'ACTA (Oral regimen vs 1-wk AmBd+5FC)' & 
                                risk_hr_spl_all$model == 'Basic' ,]) +
  scale_y_continuous(trans = 'log2') +
  coord_cartesian(xlim = c(0,0.4), ylim = c(0.25,4)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_smooth(aes(x = pred_risk, y = HR,
                  color = ifelse(pred_risk <= tercile_thresholds_b[1], "Low",
                                 ifelse(pred_risk <= tercile_thresholds_b[2], "Medium", "High"))),
              method = loess, se = F
  ) +
  geom_ribbon(aes(x = pred_risk, ymin = Lower, ymax = Upper,
                  fill = ifelse(pred_risk <= tercile_thresholds_b[1], "Low",
                                ifelse(pred_risk <= tercile_thresholds_b[2], "Medium", "High"))
  ),
  alpha = 0.2,
  show.legend = F) +
  scale_color_manual(limits = c("Low", 'Medium', "High"),
                     values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347")) +
  scale_fill_manual(limits = c("Low", 'Medium', "High"),
                    values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347"),
  ) +
  labs(color = "Risk Category",
       y = "", x = NULL) + 
  theme_pubr() +
  facet_wrap(trial~model, scales="free")
risk_treat_effect_plot_acta_basic

risk_treat_effect_plot_acta_research <- 
  ggplot(data=risk_hr_spl_all[risk_hr_spl_all$trial == 'ACTA (Oral regimen vs 1-wk AmBd+5FC)' & 
                                risk_hr_spl_all$model == 'Research' ,]) +
  scale_y_continuous(trans = 'log2') +
  coord_cartesian(xlim = c(0,0.4), ylim = c(0.25,4)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_smooth(aes(x = pred_risk, y = HR,
                  color = ifelse(pred_risk <= tercile_thresholds_r[1], "Low",
                                 ifelse(pred_risk <= tercile_thresholds_r[2], "Medium", "High"))),
              method = loess, se = F
  ) +
  geom_ribbon(aes(x = pred_risk, ymin = Lower, ymax = Upper,
                  fill = ifelse(pred_risk <= tercile_thresholds_r[1], "Low",
                                ifelse(pred_risk <= tercile_thresholds_r[2], "Medium", "High"))
  ),
  alpha = 0.2,
  show.legend = F) +
  scale_color_manual(limits = c("Low", 'Medium', "High"),
                     values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347")) +
  scale_fill_manual(limits = c("Low", 'Medium', "High"),
                    values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347"),
  ) +
  labs(color = "Risk Category",
       y = "", x = NULL) + 
  theme_pubr() +
  facet_wrap(trial~model, scales="free")
risk_treat_effect_plot_acta_research


risk_treat_effect_plot_ambition_basic <- 
  ggplot(data=risk_hr_spl_all[risk_hr_spl_all$trial == "Ambition (Ambition regimen vs 1-wk AmBd+5FC)" & 
                                risk_hr_spl_all$model == 'Basic' ,]) +
  scale_y_continuous(trans = 'log2') +
  coord_cartesian(xlim = c(0,0.4), ylim = c(0.25,4)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_smooth(aes(x = pred_risk, y = HR,
                  color = ifelse(pred_risk <= tercile_thresholds_b[1], "Low",
                                 ifelse(pred_risk <= tercile_thresholds_b[2], "Medium", "High"))),
              method = loess, se = F
  ) +
  geom_ribbon(aes(x = pred_risk, ymin = Lower, ymax = Upper,
                  fill = ifelse(pred_risk <= tercile_thresholds_b[1], "Low",
                                ifelse(pred_risk <= tercile_thresholds_b[2], "Medium", "High"))
  ),
  alpha = 0.2,
  show.legend = F) +
  scale_color_manual(limits = c("Low", 'Medium', "High"),
                     values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347")) +
  scale_fill_manual(limits = c("Low", 'Medium', "High"),
                    values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347"),
  ) +
  labs(color = "Risk Category",
       y = "", x = NULL) + 
  theme_pubr() +
  facet_wrap(trial~model, scales="free")
risk_treat_effect_plot_ambition_basic

risk_treat_effect_plot_ambition_research <- 
  ggplot(data=risk_hr_spl_all[risk_hr_spl_all$trial == "Ambition (Ambition regimen vs 1-wk AmBd+5FC)" & 
                                risk_hr_spl_all$model == 'Research' ,]) +
  scale_y_continuous(trans = 'log2') +
  coord_cartesian(xlim = c(0,0.4), ylim = c(0.25,4)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_smooth(aes(x = pred_risk, y = HR,
                  color = ifelse(pred_risk <= tercile_thresholds_r[1], "Low",
                                 ifelse(pred_risk <= tercile_thresholds_r[2], "Medium", "High"))),
              method = loess, se = F
  ) +
  geom_ribbon(aes(x = pred_risk, ymin = Lower, ymax = Upper,
                  fill = ifelse(pred_risk <= tercile_thresholds_r[1], "Low",
                                ifelse(pred_risk <= tercile_thresholds_r[2], "Medium", "High"))
  ),
  alpha = 0.2,
  show.legend = F) +
  scale_color_manual(limits = c("Low", 'Medium', "High"),
                     values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347")) +
  scale_fill_manual(limits = c("Low", 'Medium', "High"),
                    values = c("Low" = "#75C141", "Medium" = "#FFA500", "High" = "#FF6347"),
  ) +
  labs(color = "Risk Category",
       y = "", x = NULL) + 
  theme_pubr() +
  facet_wrap(trial~model, scales="free")
risk_treat_effect_plot_ambition_research


risk_treat_effect_plot_combined <- ggarrange(
  risk_treat_effect_plot_acta_basic, risk_treat_effect_plot_acta_research,
  risk_treat_effect_plot_ambition_basic, risk_treat_effect_plot_ambition_research,
  nrow = 2, ncol = 2, common.legend = T
)
risk_treat_effect_plot_combined <- annotate_figure(risk_treat_effect_plot_combined,
                                                   bottom = "Predicted mortality",
                                                   left = "Treatment effect \n expressed as HR"
)
risk_treat_effect_plot_combined
