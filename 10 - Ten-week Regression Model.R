# Cryptococcal meningitis prognostic model
# 10: 10-week mortality model train and validate
# Description: Variable selection, IECV and validation of secondary models 
# Started: 16/11/23
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

label(dev_mi_aregi$death_2week) <- "Two-week Mortality"
label(dev_mi_aregi$treatment) <- "Treatment received"
label(dev_mi_aregi$gcs_bin) <- "GCS score"
label(dev_mi_aregi$ecog) <- "ECOG"
label(dev_mi_aregi$neut) <- "Neutrophil count"
label(dev_mi_aregi$haem_gl) <- "Haemoglobin"
label(dev_mi_aregi$csf_OP) <- "CSF opening pressure"
label(dev_mi_aregi$age) <- "Age"
label(dev_mi_aregi$weight) <- "Weight"
label(dev_mi_aregi$seizure) <- "Seizures on admission"
label(dev_mi_aregi$csf_cellcount) <- "CSF White cell count"
label(dev_mi_aregi$csf_qculture_log) <- "Log(CSF quantitative culture)"

###################### BASIC MODEL ################################


# Set datadist (basic model)

attach(dev_mi_aregi)
ddist <- datadist(age, sex, weight, seizure, gcs_bin, ecog,
                  neut, haem_gl, 
                  treatment, 
                  death_10week,
                  q.display = c(0.025, 0.975))
options(datadist = 'ddist')
detach(dev_mi_aregi)
attached <- search()
attached [!grepl("package", attached)]



### Variable selection in each MI dataset using logistic regression

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

# Define the formula to input into the variable selector function (same as complete case basic formula)
f <- death_10week ~ 
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

tbl_var_sel_b_detail_10 <- var_sel_basic_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("Number of MI datasets\nselected (total n=10)")) %>%
  as_flex_table()
tbl_var_sel_b_detail_10

dev_mi_aregi_mids <- dev_mi_aregi %>%
  select(.imp, study_id,
         age, sex, weight, seizure, gcs_bin, ecog,
         neut, haem_gl,
         treatment,
         death_10week) %>%
  as.mids(.id = "study_id")

# Run the fit.mult.impute with the selected variables to create the model
basic_mi_model_10 <- fit.mult.impute(death_10week ~ rcs(age, 3) + rcs(weight, 3) +
                                       seizure + ecog + treatment +
                                       rcs(neut, 3) + rcs(haem_gl, 3),
                                     lrm,
                                     dev_mi_aregi_mids)

#Examine associations

basic_predictions_10_df <- Predict(basic_mi_model_10, fun = plogis) %>% 
  rownames_to_column(var = 'variable') %>%
  as.data.frame() %>% 
  mutate(variable = case_when(
    substr(variable, 1, 3) == 'age' ~ 'age',
    substr(variable, 1, 6) == 'weight' ~ 'weight',
    substr(variable, 1, 7) == 'seizure' ~ 'seizure',
    substr(variable, 1, 4) == 'ecog' ~ 'ecog',
    substr(variable, 1, 9) == 'treatment' ~ 'treatment',
    substr(variable, 1, 4) == 'neut' ~ 'neut',
    substr(variable, 1, 4) == 'haem' ~ 'haem',
    .default = variable
  )) %>% 
  mutate(
    x_var_type = ifelse(variable == 'neut' | variable == 'haem' |
                          variable == 'age' | variable == 'weight',
                        'numeric', 'character'),
    x = case_when(
      variable == 'age' ~ as.character(age),
      variable == 'weight' ~ as.character(weight),
      variable == 'seizure' ~ seizure,
      variable == 'ecog' ~ ecog,
      variable == 'treatment' ~ treatment,
      variable == 'neut' ~ as.character(neut),
      variable == 'haem' ~ as.character(haem_gl),
      .default = 'Failed'
    ),
    xlabel = case_when(
      variable=="age" ~ "Age (years)",
      variable=="weight" ~ "Weight (kg)",
      variable=="ecog" ~ "ECOG",
      variable=="seizure" ~ "Seizures",
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
    coord_cartesian(ylim=c(0,0.5)) +
    ylab("")
  
  factors <- dataset1 %>%  
    filter(x_var_type=="character") %>%
    ggplot() +
    geom_point(aes(x=x, y=pred), size=1) +
    ggplot2::geom_errorbar(aes(x=x, ymin=lower, ymax=upper), width=0.01, color="black") +
    facet_wrap(~xlabel, scales = "free_x", ncol=2) +
    theme_pubr() +
    xlab("") +
    coord_cartesian(ylim=c(0,0.5)) +
    ylab("") +
    theme(axis.text.x = element_text(angle = 30, 
                                     vjust = 0.7,
                                     size = 6
    ))
  
  plot_marginals_research_2w <- ggarrange(numerics, factors, 
                                          ncol=2, widths = c(2,2))
  return(plot_marginals_research_2w)
}

basic_plot_assoc_mi_10 <- get_marg_plot_logistic(dataset = basic_predictions_10_df)

basic_plot_assoc_mi_10 <- annotate_figure(
  basic_plot_assoc_mi_10,
  bottom = 'Predictor value',
  left = 'Risk of mortality (%)'
)

basic_plot_assoc_mi_10


## Create a table of model coefficients

basic_regress_table_10w <- basic_mi_model_10 %>% 
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
                              paste0("- Age = ",round(basic_mi_model_10$Design$parms$age,1)[1], ", ",
                                     round(basic_mi_model_10$Design$parms$age,1)[2],", ",
                                     round(basic_mi_model_10$Design$parms$age,1)[3],";"),
                              paste0("- Weight = ",round(basic_mi_model_10$Design$parms$weight,1)[1], ", ",
                                     round(basic_mi_model_10$Design$parms$weight,1)[2],", ",
                                     round(basic_mi_model_10$Design$parms$weight,1)[3],";"),
                              paste0("- Neutrophils = ",round(basic_mi_model_10$Design$parms$neut,1)[1], ", ",
                                     round(basic_mi_model_10$Design$parms$neut,1)[2],", ",
                                     round(basic_mi_model_10$Design$parms$neut,1)[3],";"),
                              paste0("- Haemoglobin = ",round(basic_mi_model_10$Design$parms$haem_gl,1)[1], ", ",
                                     round(basic_mi_model_10$Design$parms$haem_gl,1)[2],", ",
                                     round(basic_mi_model_10$Design$parms$haem_gl,1)[3],".")))
# Rename variables
basic_regress_table_10w[["body"]][["content"]][["data"]][[2]][['txt']] <- 'Age'
basic_regress_table_10w[["body"]][["content"]][["data"]][[3]][['txt']] <- 'Age (spline 1)'
basic_regress_table_10w[["body"]][["content"]][["data"]][[4]][['txt']] <- 'Weight'
basic_regress_table_10w[["body"]][["content"]][["data"]][[5]][['txt']] <- 'Weight (spline 1)'
basic_regress_table_10w[["body"]][["content"]][["data"]][[6]][['txt']] <- 'Seizure: Yes'
basic_regress_table_10w[["body"]][["content"]][["data"]][[7]][['txt']] <- 'ECOG: Restricted activity'
basic_regress_table_10w[["body"]][["content"]][["data"]][[8]][['txt']] <- 'ECOG: Ambulatory'
basic_regress_table_10w[["body"]][["content"]][["data"]][[9]][['txt']] <- 'ECOG: Limited self-care'
basic_regress_table_10w[["body"]][["content"]][["data"]][[10]][['txt']] <- 'ECOG: Bedbound'
basic_regress_table_10w[["body"]][["content"]][["data"]][[11]][['txt']] <- 'Treatment: 1wk AmBd+FLU'
basic_regress_table_10w[["body"]][["content"]][["data"]][[12]][['txt']] <- 'Treatment: 2wk AmBd+5FC'
basic_regress_table_10w[["body"]][["content"]][["data"]][[13]][['txt']] <- 'Treatment: 2wk AmBd+FLU'
basic_regress_table_10w[["body"]][["content"]][["data"]][[14]][['txt']] <- 'Treatment: 5FC+FLU'
basic_regress_table_10w[["body"]][["content"]][["data"]][[15]][['txt']] <- 'Neutrophils'
basic_regress_table_10w[["body"]][["content"]][["data"]][[16]][['txt']] <- 'Neutrophils (spline 1)'
basic_regress_table_10w[["body"]][["content"]][["data"]][[17]][['txt']] <- 'Haemoglobin'
basic_regress_table_10w[["body"]][["content"]][["data"]][[18]][['txt']] <- 'Haemoglobin (spline 1)'

basic_regress_table_10w

########################## RESEARCH MODEL ###################################

attach(dev_mi_aregi)
ddist <- datadist(age, sex, weight, seizure, gcs_bin, ecog,
                  neut, haem_gl, 
                  cd4, csf_OP, csf_cellcount, csf_qculture_log,
                  treatment, 
                  death_10week,
                  q.display = c(0.025, 0.975))
options(datadist = 'ddist')
detach(dev_mi_aregi)
attached <- search()
attached [!grepl("package", attached)]

### Variable selection in each MI dataset using logistic regression

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

f_r <- death_10week ~ 
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

tbl_var_sel_r_detail_10 <- var_sel_r_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("Number of MI datasets\nselected (total n=10)")) %>%
  as_flex_table()
tbl_var_sel_r_detail_10


dev_mi_aregi_mids <- dev_mi_aregi %>%
  select(.imp, study_id,
         age, sex, weight, seizure, gcs_bin, ecog,
         neut, haem_gl,
         treatment,
         cd4, csf_OP, csf_cellcount, csf_qculture_log,
         death_10week) %>%
  as.mids(.id = "study_id")

research_mi_model_10 <- fit.mult.impute(death_10week ~ rcs(age, 3) +
                                          gcs_bin + ecog + treatment +
                                          rcs(neut, 3) + rcs(haem_gl, 3) +
                                          rcs(csf_cellcount, 3) +
                                          rcs(csf_qculture_log, 3),
                                        lrm,
                                        dev_mi_aregi_mids)

#Examine associations

research_predictions_df_10 <- Predict(research_mi_model_10, fun = plogis) %>% 
  rownames_to_column(var = 'variable') %>%
  as.data.frame() %>% 
  mutate(variable = case_when(
    substr(variable, 1, 3) == 'age' ~ 'age',
    substr(variable, 1, 3) == 'gcs' ~ 'gcs',
    substr(variable, 1, 4) == 'ecog' ~ 'ecog',
    substr(variable, 1, 9) == 'treatment' ~ 'treatment',
    substr(variable, 1, 4) == 'neut' ~ 'neut',
    substr(variable, 1, 4) == 'haem' ~ 'haem',
    substr(variable, 1, 8) == 'csf_cell' ~ 'csf_cell',
    substr(variable, 1, 9) == 'csf_qcult' ~ 'csf_qcult',
    .default = variable
  )) %>% 
  mutate(
    x_var_type = ifelse(variable == 'neut' | variable == 'haem' |
                          variable == 'age' |
                          variable == 'csf_cell' | variable == 'csf_qcult',
                        'numeric', 'character'),
    x = case_when(
      variable == 'age' ~ as.character(age),
      variable == 'gcs' ~ gcs_bin,
      variable == 'ecog' ~ ecog,
      variable == 'treatment' ~ treatment,
      variable == 'neut' ~ as.character(neut),
      variable == 'haem' ~ as.character(haem_gl),
      variable == 'csf_cell' ~ as.character(csf_cellcount),
      variable == 'csf_qcult' ~ as.character(csf_qculture_log),
      .default = 'Failed'
    ),
    xlabel = case_when(
      variable=='age' ~ 'Age (years)',
      variable=="ecog" ~ "ECOG",
      variable=="gcs" ~ "GCS",
      variable=="treatment" ~ "Treatment arm",
      variable=="neut" ~ "Neutrophils (x10^9/L)",
      variable=="haem" ~ "Haemoglobin (g/L)",
      variable == 'csf_cell' ~ 'CSF cell count (WBC per mm3)',
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
    coord_cartesian(ylim=c(0,0.5)) +
    ylab("")
  
  factors <- dataset1 %>%  
    filter(x_var_type=="character") %>%
    ggplot() +
    geom_point(aes(x=x, y=pred), size=1) +
    ggplot2::geom_errorbar(aes(x=x, ymin=lower, ymax=upper), width=0.01, color="black") +
    facet_wrap(~xlabel, scales = "free_x", ncol=3) +
    theme_pubr() +
    xlab("") +
    coord_cartesian(ylim=c(0,0.5)) +
    ylab("") +
    theme(axis.text.x = element_text(angle = 0, 
                                     vjust = 0.7,
                                     # size = 6
    ))
  
  plot_marginals_research_2w <- ggarrange(factors, numerics, 
                                          ncol=1, heights = c(1.2,2.2))
  return(plot_marginals_research_2w)
}

research_plot_assoc_mi_10 <- get_marg_plot_logistic(dataset=research_predictions_df_10)

research_plot_assoc_mi_10 <- annotate_figure(
  research_plot_assoc_mi_10,
  bottom = 'Predictor value',
  left = 'Risk of mortality (%)'
)

research_plot_assoc_mi_10



## Create a table of model coefficients

research_regress_table_10w <- research_mi_model_10 %>% 
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
                              paste0("- Age = ",round(research_mi_model_10$Design$parms$age,1)[1], ", ",
                                     round(research_mi_model_10$Design$parms$age,1)[2],", ",
                                     round(research_mi_model_10$Design$parms$age,1)[3],";"),
                              paste0("- Neutrophils = ",round(research_mi_model_10$Design$parms$neut,1)[1], ", ",
                                     round(research_mi_model_10$Design$parms$neut,1)[2],", ",
                                     round(research_mi_model_10$Design$parms$neut,1)[3],";"),
                              paste0("- Haemoglobin = ",round(research_mi_model_10$Design$parms$haem_gl,1)[1], ", ",
                                     round(research_mi_model_10$Design$parms$haem_gl,1)[2],", ",
                                     round(research_mi_model_10$Design$parms$haem_gl,1)[3],"."),
                              paste0("- CSF WCC = ",round(research_mi_model_10$Design$parms$csf_cellcount,1)[1], ", ",
                                     round(research_mi_model_10$Design$parms$csf_cellcount,1)[2],", ",
                                     round(research_mi_model_10$Design$parms$csf_cellcount,1)[3],";"),
                              paste0("- log CSF QCC = ",round(research_mi_model_10$Design$parms$csf_qculture_log,1)[1], ", ",
                                     round(research_mi_model_10$Design$parms$csf_qculture_log,1)[2],", ",
                                     round(research_mi_model_10$Design$parms$csf_qculture_log,1)[3],";")))
# Rename variables
research_regress_table_10w[["body"]][["content"]][["data"]][[2]][['txt']] <- 'Age'
research_regress_table_10w[["body"]][["content"]][["data"]][[3]][['txt']] <- 'Age (spline 1)'
research_regress_table_10w[["body"]][["content"]][["data"]][[4]][['txt']] <- 'Glasgow Coma Score: 11-14'
research_regress_table_10w[["body"]][["content"]][["data"]][[5]][['txt']] <- 'Glasgow Coma Score: ≤10'
research_regress_table_10w[["body"]][["content"]][["data"]][[6]][['txt']] <- 'ECOG: Restricted activity'
research_regress_table_10w[["body"]][["content"]][["data"]][[7]][['txt']] <- 'ECOG: Ambulatory'
research_regress_table_10w[["body"]][["content"]][["data"]][[8]][['txt']] <- 'ECOG: Limited self-care'
research_regress_table_10w[["body"]][["content"]][["data"]][[9]][['txt']] <- 'ECOG: Bedbound'
research_regress_table_10w[["body"]][["content"]][["data"]][[10]][['txt']] <- 'Treatment: 1wk AmBd+FLU'
research_regress_table_10w[["body"]][["content"]][["data"]][[11]][['txt']] <- 'Treatment: 2wk AmBd+5FC'
research_regress_table_10w[["body"]][["content"]][["data"]][[12]][['txt']] <- 'Treatment: 2wk AmBd+FLU'
research_regress_table_10w[["body"]][["content"]][["data"]][[13]][['txt']] <- 'Treatment: 5FC+FLU'
research_regress_table_10w[["body"]][["content"]][["data"]][[14]][['txt']] <- 'Neutrophils'
research_regress_table_10w[["body"]][["content"]][["data"]][[15]][['txt']] <- 'Neutrophils (spline 1)'
research_regress_table_10w[["body"]][["content"]][["data"]][[16]][['txt']] <- 'Haemoglobin'
research_regress_table_10w[["body"]][["content"]][["data"]][[17]][['txt']] <- 'Haemoglobin (spline 1)'
research_regress_table_10w[["body"]][["content"]][["data"]][[18]][['txt']] <- 'CSF White Cell Count'
research_regress_table_10w[["body"]][["content"]][["data"]][[19]][['txt']] <- 'CSF White Cell Count (spline 1)'
research_regress_table_10w[["body"]][["content"]][["data"]][[20]][['txt']] <- 'CSF QCulture (log)'
research_regress_table_10w[["body"]][["content"]][["data"]][[21]][['txt']] <- 'CSF QCulture (log) (spline 1)'

research_regress_table_10w


# Save models
saveRDS(basic_mi_model_10, file = "basic_mi_model_10week_161123.rds")
saveRDS(research_mi_model_10, file = "research_mi_model_10week_161123.rds")


###### Combine variable selection tables:

tbl_var_sel_b_d10 <- var_sel_basic_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("Number of MI datasets\nselected (total n=10)"))
tbl_var_sel_r_d10 <- var_sel_r_master %>%
  filter(imp > 0) %>%
  select(-n, -imp) %>%
  tbl_summary() %>%
  modify_header(label ~ md("**Variable**"), 
                stat_0 ~ md("Number of MI datasets\nselected (total n=10)"))

table_var_sel_complete_10week <- tbl_merge(
  tbls = list(tbl_var_sel_b_d10, tbl_var_sel_r_d10),
  tab_spanner = c("**Baseline Model**", "**Research Model**")
) %>%
  as_flex_table()
table_var_sel_complete_10week



# Chosen Models for IECV

f_b_10 <- death_10week ~ rcs(age, 3) + rcs(weight, 3) + seizure + ecog + treatment +
  rcs(neut, 3) + rcs(haem_gl, 3)

# Chosen Research Model for IECV

f_r_10 <- death_10week ~ rcs(age, 3) + gcs_bin + ecog + treatment +
  rcs(neut, 3) + rcs(haem_gl, 3) + rcs(csf_cellcount, 3) + rcs(csf_qculture_log, 3)

# Remove labelling from df
dev_mi_aregi <- data.frame(dev_mi_aregi)


get_iecv_10 <- function(dataset, impno, location, f.chosen){
  
  #Split into dev and val
  ipd.dev <- dataset %>% filter(.imp == impno & country_iecv2 != location)
  ipd.val <- dataset %>% filter(.imp == impno & country_iecv2 == location)
  
  fit.dev <- lrm(f.chosen, data = ipd.dev)
  
  ipd.val$predy.val <- predict(fit.dev, type = 'fitted', newdata = ipd.val)
  ipd.val$predlp.val <- predict(fit.dev, type = 'lp', newdata = ipd.val)
  
  fit.val.slope <- glm(death_10week ~ predlp.val, data = ipd.val, family = 'binomial')
  fit.val.citl <- glm(death_10week ~ 1, offset = predlp.val, data = ipd.val, family = 'binomial')
  
  cstat.val <- roc(ipd.val$death_10week, ipd.val$predy.val, ci = T)$ci[1:3]
  
  fit.recal <- fit.dev
  fit.recal$coefficients[1] <- fit.recal$coefficients[1] + fit.val.citl$coef[1]
  ipd.val$recal.predy.val <- predict(fit.recal, type = 'fitted', newdata = ipd.val)
  ipd.val$recal.predlp.val <- predict(fit.recal, type = 'lp', newdata = ipd.val)
  
  ipd.val <- ipd.val %>% filter(!is.na(predy.val) & !is.na(death_10week)) 
  ipd.val$country_iecv <- location
  ipd.val$cstat.lb <- cstat.val[1]
  ipd.val$cstat <- cstat.val[2]
  ipd.val$cstat.ub <- cstat.val[3]
  ipd.val$cstat.se <- sqrt(var(roc(ipd.val$death_10week ~ ipd.val$predy.val)))
  ipd.val$citl <- fit.val.citl$coef[1]
  ipd.val$citl.se <- sqrt(diag(vcov(fit.val.citl)))[1]
  ipd.val$calslope <- fit.val.slope$coef[2]
  ipd.val$calslope.se <- sqrt(diag(vcov(fit.val.slope)))[2]
  ipd.val$O <- sum(ipd.val$death_10week=='Died')
  ipd.val$E <- sum(ipd.val$predy.val)
  ipd.val$N <- dim(ipd.val)[1]
  ipd.val$recal.int <- fit.recal$coefficients[1]
  ipd.val$death_10week <- as.numeric(ipd.val$death_10week) - 1
  ipd.val$crypto_model <- ipd.val$predy.val
  
  ipd.val <- ipd.val %>% select(.imp, study_id,
                                country_iecv2,
                                cstat.lb, cstat, cstat.ub, cstat.se,
                                citl, citl.se,
                                calslope, calslope.se,
                                recal.int,
                                O, E, N,
                                predy.val, predlp.val, recal.predy.val,
                                crypto_model, 
                                death_10week)
  return(ipd.val)
}

### Pool MI estimates function (agnostic of basic v research)

pool_iecv_values <- function(dataset, location){
  dataset1 <- dataset %>% filter(country_iecv2 == location)
  q1 <- dataset1 %>% select(cstat, calslope, citl)
  se1 <- dataset1 %>% select(cstat.se, calslope.se, citl.se)
  pooled <- mi.meld(q1, se1, byrow = T)
  pooled.q <- as.data.frame(pooled$q.mi)
  pooled.se <- as.data.frame(pooled$se.mi)
  iecv_master_pooled1 <- cbind(pooled.q, pooled.se)
  iecv_master_pooled1$location <- location
  iecv_master_pooled1$O <- round(mean(dataset1$O), 0)
  iecv_master_pooled1$E <- mean(dataset1$E)
  iecv_master_pooled1$N <- dataset1$N[1]
  iecv_master_pooled1 <- iecv_master_pooled1 %>% select(location, O, E, N, cstat, cstat.se,
                                                        calslope, calslope.se, citl, citl.se)
  return(iecv_master_pooled1)
}


#################### BASIC MODEL IECV ###################

iecv_master_b10 <- data.frame()

for (i in c(1:10)) {
  for (j in c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda")) {
    iecv10b <- get_iecv_10(dataset = dev_mi_aregi, f.chosen = f_b_10, impno = i, location = j)
    iecv_master_b10 <- bind_rows(iecv_master_b10, iecv10b)
  }
}

#Concise summary of IECV metrics
iecv_master_summary_b10 <- iecv_master_b10 %>%
  filter(.imp > 0) %>%
  distinct(country_iecv2, .imp, .keep_all = T) %>%
  select(country_iecv2, cstat, cstat.se, calslope, calslope.se, citl, citl.se, recal.int, O, E, N)


#Run pooled MI estimates function

iecv_master_pooled_b10 <- data.frame()

for (i in rev(c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda"))) {
  iecv_master_pooled10b <- pool_iecv_values(dataset = iecv_master_summary_b10, location = i)
  iecv_master_pooled_b10 <- rbind(iecv_master_pooled_b10, iecv_master_pooled10b)
}


################ RESEARCH MODEL IECV #########################


iecv_master_r10 <- data.frame()

for (i in c(1:10)) {
  for (j in c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda")) {
    iecv10r <- get_iecv_10(dataset = dev_mi_aregi, f.chosen = f_r_10, impno = i, location = j)
    iecv_master_r10 <- bind_rows(iecv_master_r10, iecv10r)
  }
}

#Concise summary of IECV metrics
iecv_master_summary_r10 <- iecv_master_r10 %>%
  filter(.imp > 0) %>%
  distinct(country_iecv2, .imp, .keep_all = T) %>%
  select(country_iecv2, cstat, cstat.se, calslope, calslope.se, citl, citl.se, recal.int, O, E, N)


#Run pooled MI estimates function

iecv_master_pooled_r10 <- data.frame()

for (i in rev(c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda"))) {
  iecv_master_pooled10r <- pool_iecv_values(dataset = iecv_master_summary_r10, location = i)
  iecv_master_pooled_r10 <- rbind(iecv_master_pooled_r10, iecv_master_pooled10r)
}



################## RANDOM EFFECTS META-ANALYSIS PLOTS ##################

## PLOT - BASIC MODEL

pdf("Graphs and Tables/iecv_metrics_basicmodel_10wk_v1.pdf", height=6, width=16)

#CSTAT
par(fig=c(0,0.4,0,1))
model.cstat.b.10 <- rma(yi = cstat, sei = cstat.se,  method = "REML", test = "knha", data = iecv_master_pooled_b10)
summary(model.cstat.b.10)
metafor::forest(model.cstat.b.10,
                slab=iecv_master_pooled_b10$location, xlab = "C-statistic", 
                alim=c(0.5, 1.0), cex=1, cex.lab=1, top=0,mlab = "Meta-analysis",
                xlim=c(0.3, 1.2))

## CITL

par(fig=c(0.4,0.7, 0,1), new=T)
model.CITL.b.10 <- rma(yi = citl, sei = citl.se, method = "REML",test = "knha", data = iecv_master_pooled_b10)
summary(model.CITL.b.10)
metafor::forest(model.CITL.b.10, slab=NA, xlab = "Calibration in the large", alim=c(-1.0, 1.0), 
                refline=0, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(-1.0, 2.0))
## Slope

par(fig=c(0.7,1, 0,1), new=T)
model.CS.b.10 <- rma(yi = calslope, sei = calslope.se, method = "REML", test = "knha", data = iecv_master_pooled_b10)
summary(model.CS.b.10)
metafor::forest(model.CS.b.10, slab=NA, xlab = "Calibration slope", alim=c(0.0, 2.0), 
                refline=1, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(0.0, 3.0))

title("Meta-Analysis of IECV Metrics - Basic 10wk Model", line = -1.3, outer = T, cex.main = 1.2)

dev.off()

# Clear plot matrix
plot.new()
par(mfrow = c(1, 1))

# Save heterogeneity figures
basic10w_cstat_i2 <- model.cstat.b.10$I2 %>% round(1)
basic10w_citl_i2 <- model.CITL.b.10$I2 %>% round(1)
basic10w_slope_i2 <- model.CS.b.10$I2 %>% round(1)

#### PLOT - RESEARCH MODEL

pdf("Graphs and Tables/iecv_metrics_researchmodel_10wk_v1.pdf", height=6, width=16)

#CSTAT
par(fig=c(0,0.4,0,1))
model.cstat.r.10 <- rma(yi = cstat, sei = cstat.se,  method = "REML", test = "knha", data = iecv_master_pooled_r10)
summary(model.cstat.r.10)
metafor::forest(model.cstat.r.10,
                slab=iecv_master_pooled_r10$location, xlab = "C-statistic", 
                alim=c(0.5, 1.0), cex=1, cex.lab=1, top=0,mlab = "Meta-analysis",
                xlim=c(0.3, 1.2))

## CITL

par(fig=c(0.4,0.7, 0,1), new=T)
model.CITL.r.10 <- rma(yi = citl, sei = citl.se, method = "REML",test = "knha", data = iecv_master_pooled_r10)
summary(model.CITL.r.10)
metafor::forest(model.CITL.r.10, slab=NA, xlab = "Calibration in the large", alim=c(-1.0, 1.0), 
                refline=0, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(-1.0, 2.1), at=c(-1, -0.5, 0, 0.5, 1))
## Slope

par(fig=c(0.7,1, 0,1), new=T)
model.CS.r.10 <- rma(yi = calslope, sei = calslope.se, method = "REML", test = "knha", data = iecv_master_pooled_r10)
summary(model.CS.r.10)
metafor::forest(model.CS.r.10, slab=NA, xlab = "Calibration slope", alim=c(0.0, 2.0), 
                refline=1, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(0.0, 3.0))

title("Meta-Analysis of IECV Metrics - Research 10wk Model", line = -1.3, outer = T, cex.main = 1.2)

dev.off()

plot.new()
par(mfrow = c(1, 1))

# Measures of heterogeneity
research10w_cstat_i2 <- model.cstat.r.10$I2 %>% round(1)
research10w_citl_i2 <- model.CITL.r.10$I2 %>% round(1)
research10w_slope_i2 <- model.CS.r.10$I2 %>% round(1)

# Calibration plots for IECV - Research Model
calibration_iecv_r10 <- ggplot(NULL) +
  geom_smooth(data=iecv_master_r10[iecv_master_r10$.imp>0,], 
              aes(predy.val, death_10week),
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_r10[iecv_master_r10$.imp==1,], aes(predy.val), sides="b", linewidth =0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_r10

calibration_iecv_b10 <- ggplot(NULL) +
  geom_smooth(data=iecv_master_b10[iecv_master_b10$.imp>0,], 
              aes(predy.val, death_10week),
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_b10[iecv_master_b10$.imp==1,], aes(predy.val), sides="b", linewidth =0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_b10

#Calibration plots with recalibration
calibration_iecv_recal_r10 <- ggplot(NULL) +
  geom_smooth(data=iecv_master_r10[iecv_master_r10$.imp>0,], 
              aes(recal.predy.val, death_10week),
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_r10[iecv_master_r10$.imp==1,], aes(recal.predy.val), sides="b", linewidth=0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_recal_r10




################# 10-WEEK MODEL HELD-OUT VALIDATION ######################



# Predictions
val_mi_aregi <- data.frame(val_mi_aregi)

val_mi_aregi$r_predict_10model <- predict(research_mi_model_10, newdata = val_mi_aregi, type = 'fitted')
val_mi_aregi$b_predict_10model <- predict(basic_mi_model_10, newdata = val_mi_aregi, type = 'fitted')

val_mi_aregi$death_10week <- as.numeric(val_mi_aregi$death_10week)-1

get_validation_summary <- function(dataset, impno, chosen_model){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_10week
  
  if (chosen_model == "basic") {
    dataset1$score = dataset1$b_predict_10model
  }
  else if (chosen_model == 'research') {
    dataset1$score = dataset1$r_predict_10model
  }
  else{
    stop("Invalid model choice. Please use 'basic' or 'research'.")
  }
  
  ci.auc <- roc(dataset1$outcome, dataset1$score, ci = T)$ci[1:3]
  dataset1$auc <- ci.auc[2]
  dataset1$auc.se <- sqrt(var(roc(dataset1$outcome ~ dataset1$score)))
  dataset1$auc.lower <- ci.auc[1]
  dataset1$auc.upper <- ci.auc[3]
  dataset1$auc.ci <- paste0(round(dataset1$auc,2), " (", round(dataset1$auc.lower,2), " - ", round(dataset1$auc.upper, 2), ")")
  dataset1$imp <- paste(impno)
  dataset1$O <- sum(dataset1$outcome == 1)
  dataset1$E <- sum(dataset1$score)
  dataset1$N <- dim(dataset1)[1]
  
  dataset1$logit <- log(dataset1$score/(1 - dataset1$score))
  
  fit.val.slope <- glm(outcome ~ logit, data=dataset1, family=binomial)
  dataset1$slope<- fit.val.slope$coef[2]
  dataset1$slope.se <- sqrt(diag(vcov(fit.val.slope)))[2]
  
  fit.val.citl <- glm(outcome ~ 1, offset=logit, data=dataset1, family="binomial")
  dataset1$citl <- fit.val.citl$coef[1]
  dataset1$citl.se <- sqrt(diag(vcov(fit.val.citl)))[1]
  
  return(dataset1) 
}

validation_master_r_10weekmodel <- data.frame()
validation_master_b_10weekmodel <- data.frame()

# Reseach model run validation summary function and create df summary
for (i in c(1:10)) {
  val1r <- get_validation_summary(dataset = val_mi_aregi, chosen_model = "research", impno = i)
  validation_master_r_10weekmodel <- bind_rows(validation_master_r_10weekmodel, val1r)
}

val_master_summary_r_10weekmodel <- validation_master_r_10weekmodel %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)



# Basic model run validation summary function and create df summary
for (i in c(1:10)) {
  val1b <- get_validation_summary(dataset = val_mi_aregi, chosen_model = "basic", impno = i)
  validation_master_b_10weekmodel <- bind_rows(validation_master_b_10weekmodel, val1b)
}

val_master_summary_b_10weekmodel <- validation_master_b_10weekmodel %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)


# Pool validation function

pool_validation <- function(dataset){
  dataset1 <- dataset
  q1 <- dataset1 %>% select(auc, slope, citl)
  se1 <- dataset1 %>% select(auc.se, slope.se, citl.se)
  pooled <- mi.meld(q1, se1, byrow = T)
  pooled.q <- as.data.frame(pooled$q.mi)
  pooled.se <- as.data.frame(pooled$se.mi)
  val_master_pooled1 <- cbind(pooled.q, pooled.se)
  val_master_pooled1$O <- round(mean(dataset1$O), 0)
  val_master_pooled1$E <- mean(dataset1$E)
  val_master_pooled1$N <- dataset1$N[1]
  val_master_pooled1 <- val_master_pooled1 %>% select(auc, auc.se, slope, slope.se, 
                                                      citl, citl.se, O, E, N)
  return(val_master_pooled1)
}

# Run the pool validation function for both models
validation_summary_pooled_r_10weekmodel <- pool_validation(val_master_summary_r_10weekmodel)
validation_summary_pooled_b_10weekmodel <- pool_validation(val_master_summary_b_10weekmodel)


# Function to create 95%CI

get_confidence <- function(dataset){
  dataset1<- dataset
  dataset1$auc.lower <- dataset1$auc - 1.96*(dataset1$auc.se)
  dataset1$auc.upper <- dataset1$auc + 1.96*(dataset1$auc.se)
  dataset1$slope.lower <- dataset1$slope - 1.96*(dataset1$slope.se)
  dataset1$slope.upper <- dataset1$slope + 1.96*(dataset1$slope.se)
  dataset1$citl.lower <- dataset1$citl - 1.96*(dataset1$citl.se)
  dataset1$citl.upper <- dataset1$citl + 1.96*(dataset1$citl.se)
  dataset1$auc.ci <- paste0(round(dataset1$auc,2), " (", round(dataset1$auc.lower,2), " - ",
                            round(dataset1$auc.upper, 2), ")")
  dataset1$slope.ci <- paste0(round(dataset1$slope,2), " (", round(dataset1$slope.lower,2), " - ",
                              round(dataset1$slope.upper, 2), ")")
  dataset1$citl.ci <- paste0(round(dataset1$citl,2), " (", round(dataset1$citl.lower,2), " - ",
                             round(dataset1$citl.upper, 2), ")")
  
  return(dataset1)
}

validation_summary_pooled_r_10weekmodel <- get_confidence(validation_summary_pooled_r_10weekmodel)
validation_summary_pooled_b_10weekmodel <- get_confidence(validation_summary_pooled_b_10weekmodel)

# Create a summary table

combined_valdata_10weekmodel <- bind_rows(
  validation_summary_pooled_r_10weekmodel %>% mutate(model = "Research"),
  validation_summary_pooled_b_10weekmodel %>% mutate(model = "Basic")
)

validation_results_10weekmodel <- combined_valdata_10weekmodel %>%
  select(model, auc.ci, slope.ci, citl.ci) %>%
  rename(Model=model, AUROC=auc.ci, CITL=citl.ci, Slope=slope.ci)

### Calibration curves in validation data

calibration_val_r_10week <- ggplot(NULL) +
  geom_smooth(data=validation_master_r_10weekmodel[validation_master_r_10weekmodel$.imp>0,], 
              aes(score, outcome),
              color="#E6AB02",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_r_10weekmodel[validation_master_r_10weekmodel$imp==1,], 
           aes(score), sides="b", linewidth=0.2, alpha=0.4) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 
calibration_val_r_10week

calibration_val_b_10week <- ggplot(NULL) +
  geom_smooth(data=validation_master_b_10weekmodel[validation_master_b_10weekmodel$.imp>0,], 
              aes(score, outcome),
              color="#66A61E",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_b_10weekmodel[validation_master_b_10weekmodel$imp==1,], 
           aes(score), sides="b", linewidth=0.2, alpha=0.4) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 
calibration_val_b_10week
