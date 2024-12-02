# Cryptococcal meningitis prognostic model - Periskope-CM
# 04: Variable Selection for CM model (MI aregImpute analysis)
# Description: Variable selection for prognostic model using imputed data, both basic and research versions, using aregImpute.
              # Plot associations and extract model coefficents. 
# Started: 25/10/23
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
library(broom)
library(broom.helpers)

#Read in data

all_mi_aregi <- read_rds("complete_mi_aregimpute_131123.rds")

# Split into dev and val

dev_mi_aregi <- all_mi_aregi %>% filter(trial == "ACTA" | (trial == "Ambition" & country != "Malawi"))
val_mi_aregi <- all_mi_aregi %>% filter(country == "Malawi" & trial == "Ambition")

samplesize_dev <- as.numeric(length(dev_mi_aregi$study_id[dev_mi_aregi$.imp==1]))

label(dev_mi_aregi$death_2week) <- "Two-week Mortality"
label(dev_mi_aregi$treatment) <- "Treatment received"
label(dev_mi_aregi$gcs_bin) <- "GCS score"
label(dev_mi_aregi$ecog) <- "ECOG"
label(dev_mi_aregi$neut) <- "Neutrophil count"
label(dev_mi_aregi$haem_gl) <- "Haemoglobin"
label(dev_mi_aregi$csf_OP) <- "CSF opening pressure"
label(dev_mi_aregi$csf_qculture_log) <- "Log(CSF quantitative culture)"


########################### BASIC MODEL ################################


# Set datadist

attach(dev_mi_aregi)
ddist <- datadist(age, sex, weight, seizure, gcs_bin, ecog,
                  neut, haem_gl, 
                  treatment, 
                  death_2week,
                  q.display = c(0.025, 0.975))
options(datadist = 'ddist')
detach(dev_mi_aregi)
attached <- search()
attached [!grepl("package", attached)]



### Variable selection in each MI dataset using AIC

# Define the function for selecting the variables
var_sel_basic <- function(data, im, f){
  dataset1 <- data %>% filter(.imp == im)
  model <- lrm(f, dataset1, x = T, y = T)
  bw1 <- fastbw(model)
  age <- "age" %in% bw1$names.kept
  sex <- "sex" %in% bw1$names.kept
  weight <- "weight" %in% bw1$names.kept
  seizure <- "seizure" %in% bw1$names.kept
  gcs_bin <- "gcs_bin" %in% bw1$names.kept
  ecog <- "ecog" %in% bw1$names.kept
  neut <- "neut" %in% bw1$names.kept
  haem_gl <- "haem_gl" %in% bw1$names.kept
  treatment <- "treatment" %in% bw1$names.kept
  var_sel_b <- bind_cols('Age' = age, 'Sex' = sex, 'Weight' = weight, 'Seizure' = seizure,
                         'GCS' = gcs_bin, 'ECOG' = ecog, 'Treatment regimen' = treatment,
                         'Neutrophil count' = neut, 'Haemoglobin' = haem_gl)
  var_sel_b$n <- length(model$linear.predictors)
  var_sel_b$imp <- im
  return(var_sel_b)
}

# Define the formula to input into the variable selector function

f <- death_2week ~ 
  rcs(age, 3) +
  sex + rcs(weight, 3) + 
  seizure + ecog + gcs_bin + treatment +
  rcs(neut, 3) + rcs(haem_gl, 3)

# Set the data frame to take the result of the function
var_sel_basic_master <- data.frame()

# Run the function
for (i in c(1:10)) {
  var_sel_b <- var_sel_basic(im = i, f = f, data = dev_mi_aregi)
  var_sel_basic_master <- bind_rows(var_sel_basic_master, var_sel_b)
}

# Tabulate to check which variables selected (>50% imputed datasets)

tbl_var_sel_b <- var_sel_basic_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  pivot_longer(cols = c("Age","Sex", "Weight", "Seizure", "ECOG", "GCS", "Treatment regimen", 
                        "Neutrophil count", "Haemoglobin"),
               names_to = "Variable") %>%
  group_by(Variable) %>%
  summarise(sel_imp = sum(value)) %>%
  ungroup() %>%
  mutate(sel5 = sel_imp > 5) %>%
  filter(sel5 ==T) %>%
  select(-sel5) %>%
  rename('Imputations Selected' = sel_imp)
tbl_var_sel_b

tbl_var_sel_b_detail <- var_sel_basic_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("**Number of MI datasets\nselected (total n=10)**")) %>%
  as_flex_table()
tbl_var_sel_b_detail


#Convert dev_mi_aregi to mids object so can run in the fit.mult.impute function
dev_mi_aregi <- data.frame(dev_mi_aregi)

dev_mi_aregi_mids <- dev_mi_aregi %>%
  select(.imp, study_id,
         age, sex, weight, seizure, gcs_bin, ecog,
         neut, haem_gl,
         treatment,
         death_2week) %>%
  as.mids(.id = "study_id")

# Run the fit.mult.impute with the selected variables to create the model
basic_mi_model <- fit.mult.impute(death_2week ~
                                    gcs_bin + ecog + treatment +
                                    rcs(neut, 3) + rcs(haem_gl, 3),
                                  lrm,
                                  dev_mi_aregi_mids)


##### Examine variable mortality associations in plot

# Build the plot graph manually rather than in rms, so all axis are in the same orientation

basic_predictions_df <- Predict(basic_mi_model, fun = plogis) %>% 
  rownames_to_column(var = 'variable') %>%
  as.data.frame() %>% 
  mutate(variable = case_when(
    substr(variable, 1, 3) == 'gcs' ~ 'gcs',
    substr(variable, 1, 4) == 'ecog' ~ 'ecog',
    substr(variable, 1, 9) == 'treatment' ~ 'treatment',
    substr(variable, 1, 4) == 'neut' ~ 'neut',
    substr(variable, 1, 4) == 'haem' ~ 'haem',
    .default = variable
  )) %>% 
  mutate(
    x_var_type = ifelse(variable == 'neut' | variable == 'haem',
                        'numeric', 'character'),
    x = case_when(
      variable == 'gcs' ~ gcs_bin,
      variable == 'ecog' ~ ecog,
      variable == 'treatment' ~ treatment,
      variable == 'neut' ~ as.character(neut),
      variable == 'haem' ~ as.character(haem_gl),
      .default = 'Failed'
    ),
    xlabel = case_when(
      variable=="ecog" ~ "ECOG",
      variable=="gcs" ~ "GCS",
      variable=="treatment" ~ "Treatment arm",
      variable=="neut" ~ "Neutrophils (x10^9/L)",
      variable=="haem" ~ "Haemoglobin (g/L)",
      .default = NA_character_
    ),
    x = case_when(
      xlabel!="ECOG" ~ x,
      x=="Normal" ~"0",
      x=="Restricted activity" ~"1",
      x=="Ambulatory" ~"2",
      x=="Limited self-care" ~"3",
      x=="Bedbound" ~"4")
  ) %>% 
  mutate(
    x = case_when(
      xlabel!="Treatment arm" ~ x,
      x=="1wk AmBd+5FC/Ambisome" ~"1",
      x=="1wk AmBd+FLU" ~"2",
      x=="2wks AmBd+5FC" ~"3",
      x=="2wks AmBd+FLU" ~"4",
      x=="FLU+5FC PO" ~"5")
  )


get_marg_plot_logistic <- function(dataset){
  
  dataset1 <- dataset
  dataset1$pred <- dataset1$yhat
  
  numerics <- dataset1 %>%  
    filter(x_var_type=="numeric") %>%
    ggplot() +
    geom_line(aes(x=as.numeric(x), y=pred)) +
    geom_ribbon(aes(x=as.numeric(x), ymin=lower, ymax=upper), linetype=2, alpha=0.1) +
    facet_wrap(~xlabel, scales = "free_x", ncol = 2) +
    theme_pubr() +
    xlab("") +
    coord_cartesian(ylim=c(0,0.3)) +
    ylab("")
  
  factors <- dataset1 %>%  
    filter(x_var_type=="character") %>%
    ggplot() +
    geom_point(aes(x=x, y=pred), size=1) +
    ggplot2::geom_errorbar(aes(x=x, ymin=lower, ymax=upper), width=0.01, color="black") +
    facet_wrap(~xlabel, scales = "free_x", ncol=2) +
    theme_pubr() +
    xlab("") +
    coord_cartesian(ylim=c(0,0.3)) +
    ylab("") +
    theme(axis.text.x = element_text(angle = 30, 
                                     vjust = 0.7,
                                     size = 6
    ))
  
  plot_marginals_research_2w <- ggarrange(numerics, factors, 
                                          ncol=1, heights = c(1,2))
  return(plot_marginals_research_2w)
}

basic_plot_assoc_mi <- get_marg_plot_logistic(dataset=basic_predictions_df)

basic_plot_assoc_mi <- annotate_figure(
  basic_plot_assoc_mi,
  bottom = 'Predictor value',
  left = 'Risk of mortality (%)'
)

basic_plot_assoc_mi


# Table of regression coefficients

basic_regress_table <- basic_mi_model %>% 
  tidy_and_attach(exponentiate = F, conf.int = T) %>% 
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(where(is.numeric), ~ round(. , 2))) %>% 
  mutate(ci = paste0(conf.low, ' to ', conf.high),
         p.value = case_when(
           p.value >= 0.1 ~ round(p.value, 1),
           p.value < 0.1 & p.value >= 0.01 ~ round(p.value, 2),
           p.value < 0.01 & p.value >= 0.001 ~ round(p.value, 3),
           p.value <0.001 ~ p.value
         )) %>% 
  mutate(
    p.value = ifelse(p.value <0.001, "<0.001", p.value)
  ) %>% 
  select(term, estimate, ci) %>% 
  rename(Variable = term, Estimate = estimate, `95% Confidence Int.` = ci#, 
         # `p value` = p.value
  ) %>% 
  flextable() %>% 
  add_footer_lines(values = c("Restricted cubic spline knot positions are:",
                              paste0("- Neutrophils = ",round(basic_mi_model$Design$parms$neut,1)[1], ", ",
                                     round(basic_mi_model$Design$parms$neut,1)[2],", ",
                                     round(basic_mi_model$Design$parms$neut,1)[3],";"),
                              paste0("- Haemoglobin = ",round(basic_mi_model$Design$parms$haem_gl,1)[1], ", ",
                                     round(basic_mi_model$Design$parms$haem_gl,1)[2],", ",
                                     round(basic_mi_model$Design$parms$haem_gl,1)[3],".")))

# Rename variables
basic_regress_table[["body"]][["content"]][["data"]][[2]][['txt']] <- 'Glasgow Coma Score: 11-14'
basic_regress_table[["body"]][["content"]][["data"]][[3]][['txt']] <- 'Glasgow Coma Score: ≤10'
basic_regress_table[["body"]][["content"]][["data"]][[4]][['txt']] <- 'ECOG: Restricted activity'
basic_regress_table[["body"]][["content"]][["data"]][[5]][['txt']] <- 'ECOG: Ambulatory'
basic_regress_table[["body"]][["content"]][["data"]][[6]][['txt']] <- 'ECOG: Limited self-care'
basic_regress_table[["body"]][["content"]][["data"]][[7]][['txt']] <- 'ECOG: Bedbound'
basic_regress_table[["body"]][["content"]][["data"]][[8]][['txt']] <- 'Treatment: 1wk AmBd+FLU'
basic_regress_table[["body"]][["content"]][["data"]][[9]][['txt']] <- 'Treatment: 2wk AmBd+5FC'
basic_regress_table[["body"]][["content"]][["data"]][[10]][['txt']] <- 'Treatment: 2wk AmBd+FLU'
basic_regress_table[["body"]][["content"]][["data"]][[11]][['txt']] <- 'Treatment: 5FC+FLU'
basic_regress_table[["body"]][["content"]][["data"]][[12]][['txt']] <- 'Neutrophils'
basic_regress_table[["body"]][["content"]][["data"]][[13]][['txt']] <- 'Neutrophils (spline 1)'
basic_regress_table[["body"]][["content"]][["data"]][[14]][['txt']] <- 'Haemoglobin'
basic_regress_table[["body"]][["content"]][["data"]][[15]][['txt']] <- 'Haemoglobin (spline 1)'

basic_regress_table


# Save the basic model
saveRDS(basic_mi_model, file = "basic_mi_model_251023.rds")



########################## RESEARCH MODEL ###################################

# Reattach the datadist

attach(dev_mi_aregi)
ddist <- datadist(age, sex, weight, seizure, gcs_bin, ecog,
                  neut, haem_gl, 
                  cd4, csf_OP, csf_cellcount, csf_qculture_log,
                  treatment, 
                  death_2week,
                  q.display = c(0.025, 0.975))
options(datadist = 'ddist')
detach(dev_mi_aregi)
attached <- search()
attached [!grepl("package", attached)]


### Variable selection in each MI dataset by AIC

var_sel_research <- function(data, im, f){
  dataset1 <- data %>% filter(.imp == im)
  model <- lrm(f, dataset1, x = T, y = T)
  bw1 <- fastbw(model)
  age <- "age" %in% bw1$names.kept
  sex <- "sex" %in% bw1$names.kept
  weight <- "weight" %in% bw1$names.kept
  seizure <- "seizure" %in% bw1$names.kept
  gcs_bin <- "gcs_bin" %in% bw1$names.kept
  ecog <- "ecog" %in% bw1$names.kept
  neut <- "neut" %in% bw1$names.kept
  haem_gl <- "haem_gl" %in% bw1$names.kept
  treatment <- "treatment" %in% bw1$names.kept
  cd4 <- "cd4" %in% bw1$names.kept
  csf_OP <- "csf_OP" %in% bw1$names.kept
  csf_cellcount <- "csf_cellcount" %in% bw1$names.kept
  csf_qculture_log <- "csf_qculture_log" %in% bw1$names.kept
  var_sel_r <- bind_cols('Age' = age, 'Sex' = sex, 'Weight' = weight, 'Seizure' = seizure,
                         'GCS' = gcs_bin, 'ECOG' = ecog, 'Treatment regimen' = treatment,
                         'Neutrophil count' = neut, 'Haemoglobin' = haem_gl,
                         'CD4 count' = cd4, 'CSF opening pressure' = csf_OP, 
                         'CSF cell count' = csf_cellcount, 'CSF Quantitative Culture (log)' = csf_qculture_log)
  var_sel_r$n <- length(model$linear.predictors)
  var_sel_r$imp <- im
  return(var_sel_r)
}

f_r <- death_2week ~ 
  rcs(age, 3) + sex + rcs(weight, 3) + seizure + gcs_bin + ecog + 
  rcs(neut, 3) + rcs(haem_gl, 3) + 
  treatment + 
  rcs(cd4, 3) + rcs(csf_OP, 3) + rcs(csf_cellcount, 3) + rcs(csf_qculture_log, 3)

var_sel_r_master <- data.frame()

for (i in c(1:10)) {
  var_sel_r <- var_sel_research(im = i, f = f_r, data = dev_mi_aregi)
  var_sel_r_master <- bind_rows(var_sel_r_master, var_sel_r)
}

# Tabulate to check which variables selected (>50% imputed datasets)

tbl_var_sel_r <- var_sel_r_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  pivot_longer(cols = c("Age","Sex", "Weight", "Seizure", "ECOG", "GCS", "Treatment regimen", 
                        "Neutrophil count", "Haemoglobin", "CD4 count", "CSF opening pressure", 
                        "CSF cell count", "CSF Quantitative Culture (log)"),
               names_to = "Variable") %>%
  group_by(Variable) %>%
  summarise(sel_imp = sum(value)) %>%
  ungroup() %>%
  mutate(sel5 = sel_imp > 5) %>%
  filter(sel5 ==T) %>%
  select(-sel5) %>%
  rename('Imputations Selected' = sel_imp)
tbl_var_sel_r

tbl_var_sel_r_detail <- var_sel_r_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("**Number of MI datasets\nselected (total n=10)**")) %>%
  as_flex_table()
tbl_var_sel_r_detail



#Run the Model to obtain intercepts and slopes

dev_mi_aregi_mids <- dev_mi_aregi %>%
  select(.imp, study_id,
         age, sex, weight, seizure, gcs_bin, ecog,
         neut, haem_gl,
         treatment,
         cd4, csf_OP, csf_cellcount, csf_qculture_log,
         death_2week) %>%
  as.mids(.id = "study_id")

# Force treatment even though not selected. CSF OP was selected in the original model, and so included here
research_mi_model <- fit.mult.impute(death_2week ~
                                       gcs_bin + ecog + treatment +
                                       rcs(neut, 3) + rcs(haem_gl, 3) + 
                                       rcs(csf_OP, 3) + 
                                       rcs(csf_qculture_log, 3),
                                     lrm,
                                     dev_mi_aregi_mids)


####### Plot associations with everything in same axis

research_predictions_df <- Predict(research_mi_model, fun = plogis) %>% 
  rownames_to_column(var = 'variable') %>%
  as.data.frame() %>% 
  mutate(variable = case_when(
    substr(variable, 1, 3) == 'gcs' ~ 'gcs',
    substr(variable, 1, 4) == 'ecog' ~ 'ecog',
    substr(variable, 1, 9) == 'treatment' ~ 'treatment',
    substr(variable, 1, 4) == 'neut' ~ 'neut',
    substr(variable, 1, 4) == 'haem' ~ 'haem',
    substr(variable, 1, 6) == 'csf_OP' ~ 'csf_op',
    substr(variable, 1, 9) == 'csf_qcult' ~ 'csf_qcult',
    .default = variable
  )) %>% 
  mutate(
    x_var_type = ifelse(variable == 'neut' | variable == 'haem' |
                          variable == 'csf_op' | variable == 'csf_qcult',
                        'numeric', 'character'),
    x = case_when(
      variable == 'gcs' ~ gcs_bin,
      variable == 'ecog' ~ ecog,
      variable == 'treatment' ~ treatment,
      variable == 'neut' ~ as.character(neut),
      variable == 'haem' ~ as.character(haem_gl),
      variable == 'csf_op' ~ as.character(csf_OP),
      variable == 'csf_qcult' ~ as.character(csf_qculture_log),
      .default = 'Failed'
    ),
    xlabel = case_when(
      variable=="ecog" ~ "ECOG",
      variable=="gcs" ~ "GCS",
      variable=="treatment" ~ "Treatment arm",
      variable=="neut" ~ "Neutrophils (x10^9/L)",
      variable=="haem" ~ "Haemoglobin (g/L)",
      variable == 'csf_op' ~ 'CSF Opening Pressure (cmH20)',
      variable == 'csf_qcult' ~ 'Log10 CSF Quantitative Culture (cfu/mL)',
      .default = NA_character_
    ),
    x = case_when(
      xlabel!="ECOG" ~ x,
      x=="Normal" ~"0",
      x=="Restricted activity" ~"1",
      x=="Ambulatory" ~"2",
      x=="Limited self-care" ~"3",
      x=="Bedbound" ~"4")
  ) %>% 
  mutate(
    x = case_when(
      xlabel!="Treatment arm" ~ x,
      x=="1wk AmBd+5FC/Ambisome" ~"1",
      x=="1wk AmBd+FLU" ~"2",
      x=="2wks AmBd+5FC" ~"3",
      x=="2wks AmBd+FLU" ~"4",
      x=="FLU+5FC PO" ~"5")
  )


get_marg_plot_logistic <- function(dataset){
  
  dataset1 <- dataset
  dataset1$pred <- dataset1$yhat
  
  numerics <- dataset1 %>%  
    filter(x_var_type=="numeric") %>%
    ggplot() +
    geom_line(aes(x=as.numeric(x), y=pred)) +
    geom_ribbon(aes(x=as.numeric(x), ymin=lower, ymax=upper), linetype=2, alpha=0.1) +
    facet_wrap(~xlabel, scales = "free_x", ncol = 3) +
    theme_pubr() +
    xlab("") +
    coord_cartesian(ylim=c(0,0.35)) +
    ylab("")
  
  factors <- dataset1 %>%  
    filter(x_var_type=="character") %>%
    ggplot() +
    geom_point(aes(x=x, y=pred), size=1) +
    ggplot2::geom_errorbar(aes(x=x, ymin=lower, ymax=upper), width=0.01, color="black") +
    facet_wrap(~xlabel, scales = "free_x", ncol=3) +
    theme_pubr() +
    xlab("") +
    coord_cartesian(ylim=c(0,0.35)) +
    ylab("") +
    theme(axis.text.x = element_text(angle = 0, 
                                     vjust = 0.7,
                                     # size = 6
    ))
  
  plot_marginals_research_2w <- ggarrange(factors, numerics, 
                                          ncol=1, heights = c(1.2,2.2))
  return(plot_marginals_research_2w)
}

research_plot_assoc_mi <- get_marg_plot_logistic(dataset=research_predictions_df)

research_plot_assoc_mi <- annotate_figure(
  research_plot_assoc_mi,
  bottom = 'Predictor value',
  left = 'Risk of mortality (%)'
)

research_plot_assoc_mi



#### Table of regression parameters

research_regress_table <- research_mi_model %>% 
  tidy_and_attach(exponentiate = F, conf.int = T) %>% 
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  mutate(across(where(is.numeric), ~ round(. , 2))) %>% 
  mutate(ci = paste0(conf.low, ' to ', conf.high),
         p.value = case_when(
           p.value >= 0.1 ~ round(p.value, 1),
           p.value < 0.1 & p.value >= 0.01 ~ round(p.value, 2),
           p.value < 0.01 & p.value >= 0.001 ~ round(p.value, 3),
           p.value <0.001 ~ p.value
         )) %>% 
  mutate(
    p.value = ifelse(p.value <0.001, "<0.001", p.value)
  ) %>% 
  select(term, estimate, ci#, 
         # p.value
  ) %>% 
  rename(Variable = term, Estimate = estimate, `95% Confidence Int.` = ci#, 
         # `p value` = p.value
  ) %>% 
  flextable() %>% 
  add_footer_lines(values = c("Restricted cubic spline knot positions are:",
                              paste0("- Neutrophils = ",round(research_mi_model$Design$parms$neut,1)[1], ", ",
                                     round(research_mi_model$Design$parms$neut,1)[2],", ",
                                     round(research_mi_model$Design$parms$neut,1)[3],";"),
                              paste0("- Haemoglobin = ",round(research_mi_model$Design$parms$haem_gl,1)[1], ", ",
                                     round(research_mi_model$Design$parms$haem_gl,1)[2],", ",
                                     round(research_mi_model$Design$parms$haem_gl,1)[3],";"),
                              paste0("- CSF Opening Pressure = ",round(research_mi_model$Design$parms$csf_OP,1)[1], ", ",
                                     round(research_mi_model$Design$parms$csf_OP,1)[2],", ",
                                     round(research_mi_model$Design$parms$csf_OP,1)[3],";"),
                              paste0("- CSF QCulture (log) = ",round(research_mi_model$Design$parms$csf_qculture_log,1)[1], ", ",
                                     round(research_mi_model$Design$parms$csf_qculture_log,1)[2],", ",
                                     round(research_mi_model$Design$parms$csf_qculture_log,1)[3],".")))

# Rename variables
research_regress_table[["body"]][["content"]][["data"]][[2]][['txt']] <- 'Glasgow Coma Score: 11-14'
research_regress_table[["body"]][["content"]][["data"]][[3]][['txt']] <- 'Glasgow Coma Score: ≤10'
research_regress_table[["body"]][["content"]][["data"]][[4]][['txt']] <- 'ECOG: Restricted activity'
research_regress_table[["body"]][["content"]][["data"]][[5]][['txt']] <- 'ECOG: Ambulatory'
research_regress_table[["body"]][["content"]][["data"]][[6]][['txt']] <- 'ECOG: Limited self-care'
research_regress_table[["body"]][["content"]][["data"]][[7]][['txt']] <- 'ECOG: Bedbound'
research_regress_table[["body"]][["content"]][["data"]][[8]][['txt']] <- 'Treatment: 1wk AmBd+FLU'
research_regress_table[["body"]][["content"]][["data"]][[9]][['txt']] <- 'Treatment: 2wk AmBd+5FC'
research_regress_table[["body"]][["content"]][["data"]][[10]][['txt']] <- 'Treatment: 2wk AmBd+FLU'
research_regress_table[["body"]][["content"]][["data"]][[11]][['txt']] <- 'Treatment: 5FC+FLU'
research_regress_table[["body"]][["content"]][["data"]][[12]][['txt']] <- 'Neutrophils'
research_regress_table[["body"]][["content"]][["data"]][[13]][['txt']] <- 'Neutrophils (spline 1)'
research_regress_table[["body"]][["content"]][["data"]][[14]][['txt']] <- 'Haemoglobin'
research_regress_table[["body"]][["content"]][["data"]][[15]][['txt']] <- 'Haemoglobin (spline 1)'
research_regress_table[["body"]][["content"]][["data"]][[16]][['txt']] <- 'CSF Opening Pressure'
research_regress_table[["body"]][["content"]][["data"]][[17]][['txt']] <- 'CSF Opening Pressure (spline 1)'
research_regress_table[["body"]][["content"]][["data"]][[18]][['txt']] <- 'CSF QCulture (log)'
research_regress_table[["body"]][["content"]][["data"]][[19]][['txt']] <- 'CSF QCulture (log) (spline 1)'

research_regress_table


# Save the research model
saveRDS(research_mi_model, file = "research_mi_model_251023.rds")


### Combine the variable selection tables:

## Combine variable selection tables for basic and research model:
tbl_var_sel_b_d2 <- var_sel_basic_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("Number of MI datasets\nselected (total n=10)"))
tbl_var_sel_r_d2 <- var_sel_r_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("Number of MI datasets\nselected (total n=10)"))

table_var_sel_complete <- tbl_merge(
  tbls = list(tbl_var_sel_b_d2, tbl_var_sel_r_d2),
  tab_spanner = c("**Baseline Model**", "**Research Model**")
) %>%
  as_flex_table()

table_var_sel_complete


### Produce density plots of the predictions made by both models

all_mi_aregi <- data.frame(all_mi_aregi)

all_mi_aregi$basic_predictions2 <- predict(basic_mi_model, newdata = all_mi_aregi, type = 'fitted')
all_mi_aregi$research_predictions2 <- predict(research_mi_model, newdata = all_mi_aregi, type = 'fitted')

all_mi_aregi_imputed <- all_mi_aregi %>% filter(.imp>0)

predictions_density <- ggplot(all_mi_aregi_imputed, aes(x = basic_predictions2, fill = 'basic_predictions2')) +
  geom_density(alpha = 0.4) +
  geom_density(aes(x = research_predictions2, fill = 'research_predictions2'), alpha = 0.4) +
  scale_fill_manual(values = c('red', 'blue'), name = "Model",
                    labels = c('Basic', 'Research')) +
  labs(y = "Density", x = "Predicted Mortality Risk") 
predictions_density

median_research_prediction <-  as.numeric(summary(all_mi_aregi_imputed$research_predictions2)[3])
median_basic_prediction <- as.numeric(summary(all_mi_aregi_imputed$basic_predictions2)[3])



#### Fit an ECOG only model to use in DCA and paired DeLong testing

ecog_model <- fit.mult.impute(
  death_2week ~ ecog,
  lrm,
  dev_mi_aregi_mids
)

#Save
saveRDS(ecog_model, 'ecogonly_mi_model.rds')


### Fit the no-treatment basic and research model for the treatment effect analysis

basic_model_notreat <- fit.mult.impute(
  death_2week ~ ecog + gcs_bin +
    rcs(neut, 3) + rcs(haem_gl, 3),
  lrm,
  dev_mi_aregi_mids
)

research_model_notreat <- fit.mult.impute(death_2week ~
                                            gcs_bin + ecog +
                                            rcs(neut, 3) + rcs(haem_gl, 3) + 
                                            rcs(csf_OP, 3) + 
                                            rcs(csf_qculture_log, 3),
                                          lrm,
                                          dev_mi_aregi_mids)

# Save both

saveRDS(basic_model_notreat, 'basicmodel_2w_notreat.rds')
saveRDS(research_model_notreat, 'researchmodel_2w_notreat.rds')

