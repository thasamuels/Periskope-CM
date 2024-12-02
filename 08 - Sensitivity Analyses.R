# Cryptococcal meningitis prognostic model - Periskope-CM
# 8: Sensitivity Analyses
# Description: 10 week model performance, single factor analysis
# Started: 18/12/23
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

# Read in data, split into dev and val, and load the two model objects from 06
all_mi_aregi <- read_rds("complete_mi_aregimpute_131123.rds")

dev_mi_aregi <- all_mi_aregi %>% filter(trial == "ACTA" | (trial == "Ambition" & country != "Malawi"))
val_mi_aregi <- all_mi_aregi %>% filter(country == "Malawi" & trial == "Ambition")

basic_model <- readRDS("basic_mi_model_251023.rds")
research_model <- readRDS("research_mi_model_251023.rds")

# Clear formatting
dev_mi_aregi <- data.frame(dev_mi_aregi)
val_mi_aregi <- data.frame(val_mi_aregi)

val_mi_aregi$death_10week <- as.numeric(val_mi_aregi$death_10week)-1
val_mi_aregi$death_2week <- as.numeric(val_mi_aregi$death_2week)-1

#Predictions
val_mi_aregi$crypto_model_r <- predict(research_model, newdata = val_mi_aregi, type = 'fitted')
val_mi_aregi$crypto_model_b <- predict(basic_model, newdata = val_mi_aregi, type = 'fitted')


########### MAIN MODEL PERFORMANCE FOR 10 WEEK MORTALITY ENDPOINT ###############

get_validation_summary_10w <- function(dataset, impno, chosen_model){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_10week
  
  if (chosen_model == "basic") {
    dataset1$score = dataset1$crypto_model_b
  }
  else if (chosen_model == 'research') {
    dataset1$score = dataset1$crypto_model_r
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

get_iecv <- function(dataset, impno, location, f.chosen){
  
  #Split into dev and val
  ipd.dev <- dataset %>% filter(.imp == impno & country_iecv2 != location)
  ipd.val <- dataset %>% filter(.imp == impno & country_iecv2 == location)
  
  #Train model in dev
  fit.dev <- lrm(f.chosen, data = ipd.dev)
  
  #Predictions
  ipd.val$predy.val <- predict(fit.dev, type = 'fitted', newdata = ipd.val)
  ipd.val$predlp.val <- predict(fit.dev, type = 'lp', newdata = ipd.val)
  
  #Determine calibration slope in validation data
  fit.val.slope <- glm(death_2week ~ predlp.val, data = ipd.val, family = 'binomial')
  
  #Determine calibration-in-the-large in the validation data
  fit.val.citl <- glm(death_2week ~ 1, offset = predlp.val, data = ipd.val, family = 'binomial')
  
  #Determine discriminative performance in validation data (c-statistic)
  cstat.val <- roc(ipd.val$death_2week, ipd.val$predy.val, ci = T)$ci[1:3]
  
  #Recalibrate to validation data
  fit.recal <- fit.dev
  fit.recal$coefficients[1] <- fit.recal$coefficients[1] + fit.val.citl$coef[1]
  ipd.val$recal.predy.val <- predict(fit.recal, type = 'fitted', newdata = ipd.val)
  ipd.val$recal.predlp.val <- predict(fit.recal, type = 'lp', newdata = ipd.val)
  
  #Extract statistics
  ipd.val <- ipd.val %>% filter(!is.na(predy.val) & !is.na(death_2week)) 
  ipd.val$country_iecv <- location
  ipd.val$cstat.lb <- cstat.val[1]
  ipd.val$cstat <- cstat.val[2]
  ipd.val$cstat.ub <- cstat.val[3]
  ipd.val$cstat.se <- sqrt(var(roc(ipd.val$death_2week ~ ipd.val$predy.val)))
  ipd.val$citl <- fit.val.citl$coef[1]
  ipd.val$citl.se <- sqrt(diag(vcov(fit.val.citl)))[1]
  ipd.val$calslope <- fit.val.slope$coef[2]
  ipd.val$calslope.se <- sqrt(diag(vcov(fit.val.slope)))[2]
  ipd.val$O <- sum(ipd.val$death_2week=='Died')
  ipd.val$E <- sum(ipd.val$predy.val)
  ipd.val$N <- dim(ipd.val)[1]
  ipd.val$recal.int <- fit.recal$coefficients[1]
  ipd.val$death_2week <- as.numeric(ipd.val$death_2week) - 1
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
                                death_2week)
  return(ipd.val)
}

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


validation_master_r_10w <- data.frame()
validation_master_b_10w <- data.frame()

# Reseach model run validation summary function and create df summary
for (i in c(1:10)) {
  val1r <- get_validation_summary_10w(dataset = val_mi_aregi, chosen_model = "research", impno = i)
  validation_master_r_10w <- bind_rows(validation_master_r_10w, val1r)
}

val_master_summary_r_10w <- validation_master_r_10w %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)



# Basic model run validation summary function and create df summary
for (i in c(1:10)) {
  val1b <- get_validation_summary_10w(dataset = val_mi_aregi, chosen_model = "basic", impno = i)
  validation_master_b_10w <- bind_rows(validation_master_b_10w, val1b)
}

val_master_summary_b_10w <- validation_master_b_10w %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

validation_summary_pooled_r_10w <- pool_validation(val_master_summary_r_10w)
validation_summary_pooled_b_10w <- pool_validation(val_master_summary_b_10w)

validation_summary_pooled_r_10w <- get_confidence(validation_summary_pooled_r_10w)
validation_summary_pooled_b_10w <- get_confidence(validation_summary_pooled_b_10w)

combined_valdata_10w <- bind_rows(
  validation_summary_pooled_r_10w %>% mutate(model = "Research"),
  validation_summary_pooled_b_10w %>% mutate(model = "Basic")
)

validation_results_mainmodel_10w <- combined_valdata_10w %>%
  select(model, auc.ci, slope.ci, citl.ci) %>%
  rename(Model=model, AUROC=auc.ci, CITL=citl.ci, Slope=slope.ci)

validation_table_10w <- combined_valdata_10w %>%
  select(model, auc.ci, slope.ci, citl.ci) %>%
  rename(Model=model, AUROC=auc.ci, CITL=citl.ci, Slope=slope.ci) %>%
  gt() %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_footnote(
    footnote = "Brackets represent 95% Confidence Intervals",
    locations = cells_column_labels(columns = c('AUROC', 'Slope', 'CITL'))
  ) %>%
  tab_footnote(
    footnote = "AUROC - Area Under Receiver Operator Characteristic curve",
    locations = cells_column_labels(columns = 'AUROC')
  ) %>%
  tab_footnote(
    footnote = "CITL - Calibration-in-the-large",
    locations = cells_column_labels(columns = 'CITL')
  )
validation_table_10w


calibration_val_r_10w <- ggplot(NULL) +
  geom_smooth(data=validation_master_r_10w[validation_master_r_10w$.imp>0,], 
              aes(score, outcome),
              color="#E6AB02",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_r_10w[validation_master_r_10w$imp==1,], 
           aes(score), sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_r_10w


calibration_val_b_10w <- ggplot(NULL) +
  geom_smooth(data=validation_master_b_10w[validation_master_b_10w$.imp>0,], 
              aes(score, outcome),
              color="#66A61E",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_b_10w[validation_master_b_10w$imp==1,], 
           aes(score), sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_b_10w



#################### SINGLE FACTOR IECV ANALYSIS #############################

formula_list <- list(
  ecog = death_2week ~ ecog, 
  gcs = death_2week ~ gcs_bin,
  haem = death_2week ~ haem_gl,
  neut = death_2week ~ neut,
  csfop = death_2week ~ csf_OP,
  csfqc = death_2week ~ csf_qculture_log,
  treatment = death_2week ~ treatment
)

iecv_master_list <- list()

dev_mi_aregi <- data.frame(dev_mi_aregi)

for (var_name in names(formula_list)) 
{
  f <- formula_list[[var_name]]
  
  iecv_master_single <- data.frame()
  
  for (i in c(1:10)) 
  {
    for (j in c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda")) 
    {
      iecv <- get_iecv(dataset = dev_mi_aregi, f.chosen = f, impno = i, location = j)
      iecv_master_single <- bind_rows(iecv_master_single, iecv)
    }
  }
  
  
  iecv_master_summary_single <- iecv_master_single %>%
    filter(.imp > 0) %>%
    distinct(country_iecv2, .imp, .keep_all = TRUE) %>%
    select(country_iecv2, cstat, cstat.se, calslope, calslope.se, citl, citl.se, recal.int, O, E, N)
  
  iecv_master_pooled_single <- data.frame()
  
  for (i in rev(c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda"))) 
  {
    iecv_master_pooled1 <- pool_iecv_values(dataset = iecv_master_summary_single, location = i)
    iecv_master_pooled_single <- rbind(iecv_master_pooled_single, iecv_master_pooled1)
  }
  
  iecv_result <- iecv_master_pooled_single %>%
    select(location, cstat, cstat.se) %>%
    rename(!!paste0("location.", var_name) := location,
           !!paste0("cstat.", var_name) := cstat,
           !!paste0("cstat.se.", var_name) := cstat.se)
  
  iecv_master_list[[var_name]] <- iecv_result
}

# Combine the results into a master dataframe
iecv_singlefactor_master <- bind_cols(iecv_master_list)

iecv_singlefactor_master <- iecv_singlefactor_master %>%
  select(-location.gcs, -location.haem, -location.neut, -location.csfop, -location.csfqc, -location.treatment) %>%
  rename(Location=location.ecog)

columns_to_process <- c("ecog", "gcs", "haem", "neut", "csfop", "csfqc", "treatment")

for (col in columns_to_process) {
  se_col <- paste0("cstat.se.", col)
  lower_col <- paste0("cstat.lb.", col)  
  upper_col <- paste0("cstat.ub.", col)
  ci_col <- paste0("cstat.ci.", col)
  
  iecv_singlefactor_master[[lower_col]] <- iecv_singlefactor_master[[paste0("cstat.", col)]] - 1.96 * iecv_singlefactor_master[[se_col]]
  iecv_singlefactor_master[[upper_col]] <- iecv_singlefactor_master[[paste0("cstat.", col)]] + 1.96 * iecv_singlefactor_master[[se_col]]
  
  iecv_singlefactor_master[[paste0("cstat.", col)]] <- round(iecv_singlefactor_master[[paste0("cstat.", col)]], 2)
  
  # Create a column combining cstat and CI
  iecv_singlefactor_master[[ci_col]] <- paste0(
    " (",
    round(iecv_singlefactor_master[[lower_col]], 2),
    " - ",
    round(iecv_singlefactor_master[[upper_col]], 2),
    ")"
  )
}


iecv_singlefactor_master_tableready <- iecv_singlefactor_master %>%
  select(Location, cstat.ecog, cstat.ci.ecog, cstat.gcs, cstat.ci.gcs, cstat.haem, cstat.ci.haem, 
         cstat.neut, cstat.ci.neut, cstat.csfop, cstat.ci.csfop, cstat.csfqc, cstat.ci.csfqc, 
         cstat.treatment, cstat.ci.treatment) %>%
  rename(ECOG=cstat.ecog, ECOG.ci = cstat.ci.ecog, 
         GCS=cstat.gcs, GCS.ci = cstat.ci.gcs, 
         Haemoglobin=cstat.haem, Haem.ci = cstat.ci.haem,
         Neutrophils=cstat.neut, Neut.ci = cstat.ci.neut,
         CSF_Open_pressure= cstat.csfop, CSF_OP.ci = cstat.ci.csfop, 
         CSF_Quant_Culture= cstat.csfqc, CSF_QC.ci = cstat.ci.csfqc,
         Treatment = cstat.treatment, Treatment.ci = cstat.ci.treatment)

iecv_singlefactor_master_tableready <- iecv_singlefactor_master_tableready %>%
  select(Location, ECOG, CSF_Quant_Culture, Neutrophils, GCS, CSF_Open_pressure,Haemoglobin, Treatment)

iecv_singlefactor_auc_table <- iecv_singlefactor_master_tableready %>%
  gt() %>%
  cols_label(CSF_Open_pressure = "CSF Opening Pressure",
             CSF_Quant_Culture = "CSF Quantitative Culture") %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels(c(1:8))) %>%
  tab_spanner(columns = c(2:8), label = "C-statistic") %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_spanners()) %>%
  tab_footnote(footnote = "C-statistic results coloured from a lighter to darker blue as C-statistic increases",
               locations = cells_column_spanners(spanners = "C-statistic")) %>%
  data_color(columns = c(2:8), 
             method = "numeric",
             palette = "Blues",
             domain = c(0.46,0.78)
  )
iecv_singlefactor_auc_table


############################ SINGLE PREDICTOR VALIDATION ANALYSIS #########################

# Set up variables

# Continuous
val_mi_aregi$Haemoglobin <- as.numeric(val_mi_aregi$haem_gl)
val_mi_aregi$Neutrophils <- as.numeric(val_mi_aregi$neut)
val_mi_aregi$CSF_OpenPressure <- as.numeric(val_mi_aregi$csf_OP)
val_mi_aregi$CSF_Qculture_log <- as.numeric(val_mi_aregi$csf_qculture_log)

# Factorial
val_mi_aregi$GCS <- factor(val_mi_aregi$gcs_bin, ordered = T)
val_mi_aregi$ECOG <- factor(val_mi_aregi$ecog, ordered = T)
# Cant use Treatment because val_mi_aregi has only one level

single_predictors_cat <- c("GCS", "ECOG") 
single_predictors_cont <- c("Haemoglobin", "Neutrophils", "CSF_OpenPressure", "CSF_Qculture_log")
single_predictors <- c("GCS", "ECOG", "Haemoglobin", "Neutrophils", "CSF_OpenPressure", "CSF_Qculture_log")

# Create receiving df
validation_single_master <- data.frame()

# Function to get AUCs for individual predictors by imp

get_validation_summary_single <- function(dataset, variable, impno){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_2week
  dataset1$score <- dataset1[[variable]]
  
  dataset1 <- dataset1 %>% select(study_id, score, outcome)
  
  
  ci.auc <- roc(dataset1$outcome, dataset1$score, ci = T)$ci[1:3]
  
  
  dataset1$auc <- ci.auc[2]
  dataset1$auc.se <- sqrt(var(roc(dataset1$outcome ~ dataset1$score)))
  dataset1$auc.lower <- ci.auc[1]
  dataset1$auc.upper <- ci.auc[3]
  dataset1$auc.ci <- paste0(round(dataset1$auc,2), " (", round(dataset1$auc.lower,2), " - ", round(dataset1$auc.upper, 2), ")")
  
  dataset1$variable <- paste0(variable)
  dataset1$imp <- paste(impno)
  
  return(dataset1) 
}

for (i in c(1:10)) {
  for (j in single_predictors_cont) {
    val_sp1 <- get_validation_summary_single(dataset = val_mi_aregi, variable = j, impno = i)
    validation_single_master <- bind_rows(validation_single_master, val_sp1)
  }
}

val_sp2 <- data.frame()
for (i in c(1:10)) {
  val_sp2.1 <- get_validation_summary_single(dataset = val_mi_aregi, variable = 'GCS', impno = i)
  val_sp2 <- bind_rows(val_sp2, val_sp2.1)
}

val_sp3 <- data.frame()
for (i in c(1:10)) {
  val_sp3.1 <- get_validation_summary_single(dataset = val_mi_aregi, variable = 'ECOG', impno = i)
  val_sp3 <- bind_rows(val_sp3, val_sp3.1)
}

# Change levels of categorical variables to numeric to allow combining of data
# GCS - 15 = 1, 11-14 = 2 ≤10 = 3
val_sp2 <- val_sp2 %>%
  mutate(
    score = as.numeric(score)
  )

# ECOG - 0 = Normal and so on
val_sp3 <- val_sp3 %>%
  mutate(
    score = (as.numeric(score) - 1)
  )

validation_single_master <- validation_single_master %>%
  bind_rows(validation_single_master, val_sp2, val_sp3)

validation_summary_imp <- validation_single_master %>% 
  filter(imp>0) %>%
  distinct(variable, imp, .keep_all = T) %>%
  select(variable, imp, auc, auc.se)


pool_validation_single <- function(dataset, variable_name){
  dataset1 <- dataset %>% filter(variable==variable_name)
  q1 <- dataset1 %>% select(auc)
  se1 <- dataset1 %>% select(auc.se)
  pooled <- mi.meld(q=q1, se=se1, byrow = TRUE)
  pooled.q <- as.data.frame(pooled$q.mi)
  pooled.ses <- as.data.frame(pooled$se.mi)
  pooled1 <- cbind(pooled.q, pooled.ses)
  pooled1$variable_name <- paste(variable_name)
  
  pooled1 <- pooled1 %>% select(variable_name, auc, auc.se)
  return(pooled1)
}

validation_summary_single_pooled <- data.frame()

for (i in c(all_of(single_predictors))) {
  pooled1 <- pool_validation_single(dataset = validation_summary_imp, variable_name = i)
  validation_summary_single_pooled <- bind_rows(validation_summary_single_pooled, pooled1)
}

validation_summary_single_pooled$auc.lower <- validation_summary_single_pooled$auc - 
  1.96*(validation_summary_single_pooled$auc.se)
validation_summary_single_pooled$auc.upper <- validation_summary_single_pooled$auc + 
  1.96*(validation_summary_single_pooled$auc.se)
validation_summary_single_pooled$auc.ci <- paste0(round(validation_summary_single_pooled$auc,2), 
                                                  " (", round(validation_summary_single_pooled$auc.lower,2), 
                                                  " - ", round(validation_summary_single_pooled$auc.upper, 2), ")")

validation_single_predictors <- validation_summary_single_pooled


validation_single_predictors <- validation_single_predictors %>%
  select(variable_name, auc.ci) %>%
  rename(Predictor = variable_name, AUROC = auc.ci) %>%
  mutate(Predictor = case_when(
    Predictor == "GCS" ~ "Glasgow Coma Score",
    Predictor == "ECOG" ~ "ECOG Performance Status",
    Predictor == "Haemoglobin" ~ "Haemoglobin (g/L)",
    Predictor == "Neutrophils" ~ "Neutrophils (x10^9/L)",
    Predictor == "CSF_OpenPressure" ~ "CSF Opening Pressure (cmH20)",
    Predictor == "CSF_Qculture_log" ~ "CSF Quantitative Culture (log cfu/ml)"
  )) %>%
  arrange(desc(AUROC))

## Tabulate

single_predictor_table <- validation_single_predictors %>%
  gt() %>%
  cols_label(AUROC = 'C-statistic') %>% 
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
  # tab_footnote(footnote = "AUROC = Area under the receiver operator characteristic curve",
  #              locations = cells_column_labels(columns = c(2))) %>%
  tab_footnote(footnote = "Brackets show 95% confidence intervals",
               locations = cells_column_labels(columns = c(2)))
single_predictor_table
