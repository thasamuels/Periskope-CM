# Cryptococcal meningitis prognostic model - Periskope-CM
# 06: Held-out Validation
# Description: Validation of the model on held out validation cohort
# Started: 27/10/23
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

# Read in data, split into dev and val, and load the model objects

all_mi_aregi <- read_rds("complete_mi_aregimpute_131123.rds")

dev_mi_aregi <- all_mi_aregi %>% filter(trial == "ACTA" | (trial == "Ambition" & country != "Malawi"))
val_mi_aregi <- all_mi_aregi %>% filter(country == "Malawi" & trial == "Ambition")

basic_model <- readRDS("basic_mi_model_251023.rds")
research_model <- readRDS("research_mi_model_251023.rds")
ecog_model <- readRDS('ecogonly_mi_model.rds') # for benchmarking

## Generate predictions
val_mi_aregi <- data.frame(val_mi_aregi)
samplesize_val <- as.numeric(length(val_mi_aregi$study_id[val_mi_aregi$.imp==1])) #sample size object for validation

#Research
val_mi_aregi$crypto_model_r <- predict(research_model, newdata = val_mi_aregi, type = 'fitted')
val_mi_aregi$death_2week <- as.numeric(val_mi_aregi$death_2week)-1

#Basic and ECOG
val_mi_aregi$crypto_model_b <- predict(basic_model, newdata = val_mi_aregi, type = 'fitted')
val_mi_aregi$crypto_model_ecog <- predict(ecog_model, newdata = val_mi_aregi, type = 'fitted')



################### EXTERNAL VALIDATION ######################

### Validation function

get_validation_summary <- function(dataset, impno, chosen_model){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_2week
  
  if (chosen_model == "basic") {
    dataset1$score = dataset1$crypto_model_b
  }
  else if (chosen_model == 'research') {
    dataset1$score = dataset1$crypto_model_r
  }
  else if (chosen_model == 'ecog') {
    dataset1$score = dataset1$crypto_model_ecog
  }
  else{
    stop("Invalid model choice. Please use 'basic', 'research' or 'ecog'.")
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

validation_master_r <- data.frame()
validation_master_b <- data.frame()
validation_master_e <- data.frame()



# Research model run validation summary function and create df summary
for (i in c(1:10)) {
  val1r <- get_validation_summary(dataset = val_mi_aregi, chosen_model = "research", impno = i)
  validation_master_r <- bind_rows(validation_master_r, val1r)
}

saveRDS(validation_master_r, 'validation_master_dataframe_researchmodel.rds')

val_master_summary_r <- validation_master_r %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)


# Basic model run validation summary function and create df summary
for (i in c(1:10)) {
  val1b <- get_validation_summary(dataset = val_mi_aregi, chosen_model = "basic", impno = i)
  validation_master_b <- bind_rows(validation_master_b, val1b)
}

saveRDS(validation_master_b, 'validation_master_dataframe_basicmodel.rds')

val_master_summary_b <- validation_master_b %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

# ecog model run summary for dca analysis

for (i in c(1:10)) {
  val1e <- get_validation_summary(dataset = val_mi_aregi, chosen_model = "ecog", impno = i)
  validation_master_e <- bind_rows(validation_master_e, val1e)
}

val_master_summary_e <- validation_master_e %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

saveRDS(validation_master_e, 'validation_master_dataframe_ecogmodel.rds')


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
validation_summary_pooled_r <- pool_validation(val_master_summary_r)
validation_summary_pooled_b <- pool_validation(val_master_summary_b)
validation_summary_pooled_e <- pool_validation(val_master_summary_e)

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

validation_summary_pooled_r <- get_confidence(validation_summary_pooled_r)
validation_summary_pooled_b <- get_confidence(validation_summary_pooled_b)
validation_summary_pooled_e <- get_confidence(validation_summary_pooled_e)

# Create a summary table

combined_valdata <- bind_rows(
  validation_summary_pooled_r %>% mutate(model = "Research"),
  validation_summary_pooled_b %>% mutate(model = "Basic"),
  validation_summary_pooled_e %>% mutate(model = 'ECOG')
)

validation_table <- combined_valdata %>%
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
validation_table


################### BENCHMARKING TO ECOG ########################

# Compare the AUROCs using paired DeLongs
# Extract original DeLong values, by MI dataset

get_delong_summary <- function(dataset, impno){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_2week
  dataset1$imp <- paste(impno)
  
  roc_research <- roc(dataset1$outcome, dataset1$crypto_model_r)
  roc_basic <- roc(dataset1$outcome, dataset1$crypto_model_b)
  roc_ecog <- roc(dataset1$outcome, dataset1$crypto_model_ecog)
  
  delong_research <- roc.test(roc_research, roc_ecog, method = 'delong')
  delong_basic <- roc.test(roc_basic, roc_ecog, method = 'delong')
  
  dataset1$delta_auroc_research <- delong_research$estimate[[1]] - delong_research$estimate[[2]]
  dataset1$delta_auroc_basic <- delong_basic$estimate[[1]] - delong_basic$estimate[[2]]
  dataset1$z_delta_research <- (delong_research$statistic)
  dataset1$z_delta_basic <- (delong_basic$statistic)
  dataset1$se_delta_research <- (dataset1$delta_auroc_research) / (delong_research$statistic)
  dataset1$se_delta_basic <- (dataset1$delta_auroc_basic) / (delong_basic$statistic)
  
  return(dataset1)
}

delong_master <- data.frame()

for (i in c(1:10)) {
  val_del <- get_delong_summary(dataset = val_mi_aregi, impno = i)
  delong_master <- bind_rows(delong_master, val_del)
}
delong_master_summary <- delong_master %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, delta_auroc_research, delta_auroc_basic,
         se_delta_research, se_delta_basic)

# Pool results
pool_delong <- function(dataset){
  dataset1 <- dataset
  q1 <- dataset1 %>% select(delta_auroc_research, delta_auroc_basic)
  se1 <- dataset1 %>% select(se_delta_research, se_delta_basic)
  pooled <- mi.meld(q1, se1, byrow = T)
  pooled.q <- as.data.frame(pooled$q.mi)
  pooled.se <- as.data.frame(pooled$se.mi)
  val_master_pooled1 <- cbind(pooled.q, pooled.se)

  return(val_master_pooled1)
}

pooled_delong <- pool_delong(delong_master_summary)

# add CIs and p values
delong_r <- pooled_delong %>% 
  select(delta_auroc_research, se_delta_research) %>% 
  rename(delta = delta_auroc_research,
         se = se_delta_research) %>% 
  mutate(model = 'Research')

delong_b <- pooled_delong %>% 
  select(delta_auroc_basic, se_delta_basic) %>% 
  rename(delta = delta_auroc_basic,
         se = se_delta_basic) %>% 
  mutate(model = 'Basic')

delong <- delong_r %>% add_row(delong_b)

delong <- delong %>%
  mutate(
    z_score = delta / se,
    upper_ci = delta + (1.96 * se),
    lower_ci = delta - (1.96 * se)
  ) %>% 
  mutate(
    p_value = 2 * (1 - pnorm(abs(z_score)))
  )

delong_table <- delong %>%
  mutate(delta = round(delta,2)) %>% 
  mutate(
    cstat_diff_ci = 
      paste0(delta, " (", round(lower_ci, 2), " to ", round(upper_ci, 2), ")"),
    p_value = case_when(
      p_value >= 0.1 ~ round(p_value, 1),
      p_value < 0.1 & p_value >= 0.01 ~ round(p_value, 3),
      p_value < 0.01 & p_value >= 0.001 ~ round(p_value, 3),
      p_value <0.001 ~ p_value)
  ) %>% 
  select(model, cstat_diff_ci, p_value)
delong_table



# For table combine with other model validation results
validation_results_main <- combined_valdata %>%
  filter(model != 'ECOG') %>% 
  select(model, auc.ci, slope.ci, citl.ci) %>%
  rename(Model=model, AUROC=auc.ci, CITL=citl.ci, Slope=slope.ci)



######################## CALIBRATION PLOTS ########################

# Pooled calibration plot - research

calibration_val_r <- ggplot(NULL) +
  geom_smooth(data=validation_master_r[validation_master_r$.imp>0,], 
              aes(score, outcome),
              color="#E6AB02",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_r[validation_master_r$imp==1,], 
           aes(score), sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_r

# By MI dataset
validation_master_r$imp <- as.numeric(as.character(validation_master_r$imp))

calibration_val_each_mi <- ggplot(validation_master_r, 
                                  aes(score, outcome)) +
  geom_smooth(method = "loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  facet_wrap(~imp) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_each_mi


# Pooled calibration plot - basic

calibration_val_b <- ggplot(NULL) +
  geom_smooth(data=validation_master_b[validation_master_b$.imp>0,], 
              aes(score, outcome),
              color="#66A61E",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_b[validation_master_b$imp==1,], 
           aes(score), sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_b


# By MI dataset - basic

validation_master_b$imp <- as.numeric(as.character(validation_master_b$imp))

calibration_val_each_mi_b <- ggplot(validation_master_b, 
                                    aes(score, outcome)) +
  geom_smooth(method = "loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  facet_wrap(~imp) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_each_mi_b



###################### OPTIONAL RE-CALIBRATION BY INTERCEPT ##########################


# Optional re-calibration of intercept for Malawi validation data

val_recal <- function(dataset, impno, chosen_model){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_2week
  
  #Generate predictions depending on input model, and define relevant model object. 
  if (chosen_model == "basic") {
    dataset1$predy = predict(basic_model, type = "fitted", newdata = dataset1)
    dataset1$predlp = predict(basic_model, type = "lp", newdata = dataset1)
    fit.recal <- basic_model
    fit.citl = glm(outcome ~ 1, offset = predlp, data = dataset1, family = "binomial")
    fit.recal$coefficients[1] <- fit.recal$coefficients[1] + fit.citl$coef[1]
    dataset1$recal.predy <- predict(fit.recal, type = 'fitted', newdata = dataset1)
  }
  else if (chosen_model == 'research') {
    dataset1$predy = predict(research_model, type = "fitted", newdata = dataset1)
    dataset1$predlp = predict(research_model, type = "lp", newdata = dataset1)
    fit.recal <- research_model
    fit.citl = glm(outcome ~ 1, offset = predlp, data = dataset1, family = "binomial")
    fit.recal$coefficients[1] <- fit.recal$coefficients[1] + fit.citl$coef[1]
    dataset1$recal.predy <- predict(fit.recal, type = 'fitted', newdata = dataset1)
  }
  else{
    stop("Invalid model choice. Please use 'basic' or 'research'.")
  }
  
  return(dataset1)
}

# Set the receiving data frames and delabel to avoid the categorical variable level recognition problem, then recode outcome
recal_dataset_r <- data.frame()
recal_dataset_b <- data.frame()
val_mi_aregi <- data.frame(val_mi_aregi)

# Research
for (i in c(1:10)) {
  recal_dataset1r <- val_recal(dataset = val_mi_aregi, impno = i, chosen_model = "research")
  recal_dataset_r <- bind_rows(recal_dataset_r, recal_dataset1r)
}

# Basic
for (i in c(1:10)) {
  recal_dataset1b <- val_recal(dataset = val_mi_aregi, impno = i, chosen_model = "basic")
  recal_dataset_b <- bind_rows(recal_dataset_b, recal_dataset1b)
}


# Re-calibrated calibration plots for validation data

#Research

calibration_r_recal <- ggplot(NULL) +
  geom_smooth(data = recal_dataset_r[recal_dataset_r$.imp>0,],
              aes(recal.predy, death_2week),
              method = "loess", se = T,
              method.args = loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim = c(0,1), xlim = c(0,1)) +
  geom_rug(data = recal_dataset_r[recal_dataset_r$.imp == 1,], 
           aes(recal.predy), sides = 'b', linewidth = 0.1, alpha = 0.3) +
  geom_abline(color = "grey50") +
  theme_pubr() +
  xlab("Predicted Risk") +
  ylab("Observed Risk")

calibration_r_recal


#Basic
calibration_b_recal <- ggplot(NULL) +
  geom_smooth(data = recal_dataset_b[recal_dataset_b$.imp>0,],
              aes(recal.predy, death_2week),
              method = "loess", se = T,
              method.args = loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim = c(0,1), xlim = c(0,1)) +
  geom_rug(data = recal_dataset_b[recal_dataset_b$.imp == 1,], 
           aes(recal.predy), sides = 'b', linewidth = 0.1, alpha = 0.3) +
  geom_abline(color = "grey50") +
  theme_pubr() +
  xlab("Predicted Risk") +
  ylab("Observed Risk")

calibration_b_recal


######################## PREDICTION DENSITIES ########################

#### Density plot of predictions

val_mi_aregi <- data.frame(val_mi_aregi)
val_mi_aregi$research_predictions <- predict(research_model, type = 'fitted', newdata = val_mi_aregi)
val_mi_aregi$basic_predictions <- predict(basic_model, type = 'fitted', newdata = val_mi_aregi)

density_predict_val <- val_mi_aregi %>%
  filter(.imp >0)

iecv_prediction_density_val <- ggplot(data = density_predict_val, aes(x = research_predictions, fill = as.factor(.imp))) +
  geom_density(alpha = 0.4) + 
  labs(title = "Prediction distribution in Validation data by MI dataset",
       x = "Predicted probability",
       y = "Density") +
  scale_fill_discrete(name = "MI dataset")
iecv_prediction_density_val

val_mi_aregi <- val_mi_aregi %>%
  mutate(death_2week = ifelse(death_2week == 1, "Died", "Alive")) %>%
  mutate(death_2week = factor(death_2week, levels = c("Died", "Alive")))

density_predict_val <- val_mi_aregi %>%
  filter(.imp >0)

prediction_density_research_val <- ggplot(data = density_predict_val, 
                                          aes(x = research_predictions, fill = death_2week)) +
  geom_density(alpha = 0.6, size = 0.5) +
  labs(x = "Predicted Probability of Mortality",
       y = "Density") +
  scale_fill_manual(values = c('#E6AB02', '#FFD94A'), name = "2-week mortality") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 9)) +
  theme_pubr()
prediction_density_research_val

prediction_density_basic_val <- ggplot(data = density_predict_val, 
                                       aes(x = basic_predictions, fill = death_2week)) +
  geom_density(alpha = 0.6, size = 0.5) +
  labs(x = "Predicted Probability of Mortality",
       y = "Density") +
  scale_fill_manual(values = c('#508517', '#A5F184'), name = "2-week mortality") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 9)) +
  theme_pubr()
prediction_density_basic_val

### Median risk predictions by outcome

median_r_pred_death <- val_mi_aregi %>% filter(!is.na(crypto_model_r) & death_2week == 1) %>% 
  summarise(median_val = median(crypto_model_r),
            q1 = quantile(crypto_model_r, 0.25),
            q3 = quantile(crypto_model_r, 0.75))
median_r_pred_alive <- val_mi_aregi %>% filter(!is.na(crypto_model_r) & death_2week == 0) %>% 
  summarise(median_val = median(crypto_model_r),
            q1 = quantile(crypto_model_r, 0.25),
            q3 = quantile(crypto_model_r, 0.75))

median_b_pred_death <- val_mi_aregi %>% filter(!is.na(crypto_model_b) & death_2week == 1) %>% 
  summarise(median_val = median(crypto_model_b),
            q1 = quantile(crypto_model_r, 0.25),
            q3 = quantile(crypto_model_b, 0.75))
median_b_pred_alive <- val_mi_aregi %>% filter(!is.na(crypto_model_b) & death_2week == 0) %>% 
  summarise(median_val = median(crypto_model_b),
            q1 = quantile(crypto_model_r, 0.25),
            q3 = quantile(crypto_model_b, 0.75))



##################### External Validation stratified by gender, age ####################


# by gender

val_mi_female <- val_mi_aregi %>% 
  filter(sex == 'Female')

validation_master_fem_r <- data.frame()
validation_master_fem_b <- data.frame()

for (i in c(1:10)) {
  val1fr <- get_validation_summary(dataset = val_mi_female, chosen_model = "research", impno = i)
  validation_master_fem_r <- bind_rows(validation_master_fem_r, val1fr)
}

validation_summary_fem_r <- validation_master_fem_r %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

for (i in c(1:10)) {
  val1fb <- get_validation_summary(dataset = val_mi_female, chosen_model = "basic", impno = i)
  validation_master_fem_b <- bind_rows(validation_master_fem_b, val1fb)
}

validation_summary_fem_b <- validation_master_fem_b %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)


validation_summary_pooled_fem_r <- pool_validation(validation_summary_fem_r)
validation_summary_pooled_fem_b <- pool_validation(validation_summary_fem_b)

validation_summary_pooled_fem_r <- get_confidence(validation_summary_pooled_fem_r)
validation_summary_pooled_fem_b <- get_confidence(validation_summary_pooled_fem_b)

validation_female_r <- validation_summary_pooled_fem_r %>% 
  mutate(model = 'Research',
         strat = 'Female',
         variable = 'Sex') %>% 
  select(variable, strat, model, auc.ci)

validation_female_b <- validation_summary_pooled_fem_b %>% 
  mutate(model = 'Basic',
         strat = 'Female',
         variable = 'Sex') %>% 
  select(variable, strat, model, auc.ci)

validation_strat <- validation_female_r %>% 
  add_row(validation_female_b)


# Male

val_mi_male <- val_mi_aregi %>% 
  filter(sex == 'Male')

validation_master_male_r <- data.frame()
validation_master_male_b <- data.frame()

for (i in c(1:10)) {
  val1mr <- get_validation_summary(dataset = val_mi_male, chosen_model = "research", impno = i)
  validation_master_male_r <- bind_rows(validation_master_male_r, val1mr)
}

validation_summary_male_r <- validation_master_male_r %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

for (i in c(1:10)) {
  val1mb <- get_validation_summary(dataset = val_mi_male, chosen_model = "basic", impno = i)
  validation_master_male_b <- bind_rows(validation_master_male_b, val1mb)
}

validation_summary_male_b <- validation_master_male_b %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)


validation_summary_pooled_male_r <- pool_validation(validation_summary_male_r)
validation_summary_pooled_male_b <- pool_validation(validation_summary_male_b)

validation_summary_pooled_male_r <- get_confidence(validation_summary_pooled_male_r)
validation_summary_pooled_male_b <- get_confidence(validation_summary_pooled_male_b)

validation_male_r <- validation_summary_pooled_male_r %>% 
  mutate(model = 'Research',
         strat = 'Male',
         variable = 'Sex') %>% 
  select(variable, strat, model, auc.ci)

validation_male_b <- validation_summary_pooled_male_b %>% 
  mutate(model = 'Basic',
         strat = 'Male',
         variable = 'Sex') %>% 
  select(variable, strat, model, auc.ci)

validation_strat <- validation_strat %>% 
  add_row(validation_male_r) %>% 
  add_row(validation_male_b)


### by Age

val_mi_aregi <- val_mi_aregi %>% 
  mutate(
    age_quart = case_when(
      ntile(age, 2) == 1 ~ "<38", # splits in various centile groups, defined by user
      ntile(age, 2) == 2 ~ "38+"
    )
  )

age_groups <- c("<38", "38+")

get_validation_summary_agestrat <- function(dataset, impno, chosen_model, age){
  dataset1 <- dataset %>% filter(.imp==impno)
  age_select <- paste0(age)
  dataset1 <- dataset1 %>% filter(age_quart == age_select)
  
  dataset1$outcome <- dataset1$death_2week
  
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
  dataset1$agecat <- paste(age)
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

validation_master_age_r <- data.frame()
validation_master_age_b <- data.frame()

for (i in c(1:10)) {
  for (j in age_groups) {
    val1ar <- get_validation_summary_agestrat(dataset = val_mi_aregi, age = j, 
                                              chosen_model = "research", impno = i)
    validation_master_age_r <- bind_rows(validation_master_age_r, val1ar)
  }
}

validation_summary_age_r <- validation_master_age_r %>% 
  filter(imp>0) %>%
  distinct(age_quart, imp, .keep_all = T) %>%
  select(age_quart, imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

for (i in age_groups) {
  for (j in c(1:10)) {
    val1ab <- get_validation_summary_agestrat(dataset = val_mi_aregi, age = i, 
                                              chosen_model = "basic", impno = j)
    validation_master_age_b <- bind_rows(validation_master_age_b, val1ab)
  }
}

validation_summary_age_b <- validation_master_age_b %>% 
  filter(imp>0) %>%
  distinct(age_quart, imp, .keep_all = T) %>%
  select(age_quart, imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

pool_validation_agestrat <- function(dataset, age){
  agecat <- paste(age)
  dataset1 <- dataset %>% filter(age_quart==agecat)
  q1 <- dataset1 %>% select(auc, slope, citl)
  se1 <- dataset1 %>% select(auc.se, slope.se, citl.se)
  pooled <- mi.meld(q=q1, se=se1, byrow = TRUE)
  pooled.q <- as.data.frame(pooled$q.mi)
  pooled.ses <- as.data.frame(pooled$se.mi)
  pooled1 <- cbind(pooled.q, pooled.ses)
  pooled1$age_cat <- paste(age)
  pooled1 <- pooled1 %>% select(age_cat, auc, auc.se, slope, slope.se, citl, citl.se)
  return(pooled1)
}

validation_summary_pooled_age_r <- data.frame()
validation_summary_pooled_age_b <- data.frame()

for(i in age_groups){
  pooled1 <- pool_validation_agestrat(dataset=validation_summary_age_r, age = i)
  validation_summary_pooled_age_r <- bind_rows(pooled1, validation_summary_pooled_age_r)
}

for(i in age_groups){
  pooled1b <- pool_validation_agestrat(dataset=validation_summary_age_b, age = i)
  validation_summary_pooled_age_b <- bind_rows(pooled1b, validation_summary_pooled_age_b)
}

validation_summary_pooled_age_r <- get_confidence(validation_summary_pooled_age_r)
validation_summary_pooled_age_b <- get_confidence(validation_summary_pooled_age_b)


validation_age_r <- validation_summary_pooled_age_r %>% 
  mutate(model = 'Research',
         strat = age_cat,
         variable = 'Age') %>% 
  select(variable, strat, model, auc.ci)

validation_age_b <- validation_summary_pooled_age_b %>% 
  mutate(model = 'Basic',
         strat = age_cat,
         variable = 'Age') %>% 
  select(variable, strat, model, auc.ci)

validation_strat <- validation_strat %>% 
  add_row(validation_age_r) %>% 
  add_row(validation_age_b)

add_pool_strat <- combined_valdata %>%
  filter(model != 'ECOG') %>%
  mutate(strat = '-',
         variable = 'Pooled') %>% 
  select(variable, strat, model, auc.ci)

validation_strat <- validation_strat %>% 
  add_row(add_pool_strat) %>%
  mutate(variable = factor(variable, levels = c('Pooled', 'Sex', 'Age'))) %>% 
  arrange(variable) %>% 
  pivot_wider(
    id_cols = c(variable, strat),
    names_from = model,
    values_from = auc.ci
  )

# Sample sizes

val_ss <- dim(val_mi_aregi[val_mi_aregi$.imp == 1,])[1]

fem_ss <- dim(val_mi_female[val_mi_female$.imp == 1,])[1]
male_ss <- dim(val_mi_male[val_mi_male$.imp == 1,])[1]

age1_ss <- val_mi_aregi %>% filter(age_quart == '<38' & .imp == 1) %>% dim()
age1_ss <- age1_ss[1]
age2_ss <- val_mi_aregi %>% filter(age_quart == '38+' & .imp == 1) %>% dim()
age2_ss <- age2_ss[1]

sample_size <- data.frame(sample_size = c(val_ss, fem_ss, male_ss, age1_ss, age2_ss))

validation_strat <- validation_strat %>% 
  add_column(sample_size) %>%
  select(variable, strat, sample_size, Basic, Research)

