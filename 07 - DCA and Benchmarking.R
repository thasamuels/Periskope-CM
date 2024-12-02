# Cryptococcal meningitis prognostic model - Periskope-CM
# 7: Comparative Analysis with Existing Models and Decision Curve Analysis
# Description: Validating Zhao et al model on validation dataset, DCA for all models 
# Started: 11/03/2024
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
library(gt)
library(Hmisc)
library(readxl)
library(gtools)

#Read in data
all_mi <- readRDS("complete_mi_aregimpute_131123.rds")

# Clear labels
all_mi <- data.frame(all_mi)

# Create variable categories and Zhao score

all_mi <- all_mi %>%
  mutate(
    age_zhao = case_when(
      age <= 49 ~ '≤49',
      age >= 50 & age <= 59 ~ '50-59',
      age >= 60 & age <= 69 ~ '60-69',
      age >= 70 ~ '≥70',
      T ~ NA_character_
    ) %>% factor(levels = c('≤49', '50-59', '60-69', '≥70')),
    gcs_zhao = case_when(
      gcs_bin == '15' ~ 'No',
      gcs_bin != '15' ~ 'Yes',
      T ~ NA_character_
    ) %>% factor(levels = c('No', 'Yes')),
    csf_op_zhao_cat = case_when(
      csf_OP < 20 ~ '<200',
      csf_OP >= 20 & csf_OP <= 25 ~ '200-250',
      csf_OP > 25 & csf_OP <= 35 ~ '251-350',
      csf_OP > 35 ~ '>350',
      T ~ NA_character_
    ) %>% factor(levels = c('<200', '200-250', '251-350', '>350')),
    csf_op_zhao_cont = csf_OP * 10,
    urea_zhao = case_when(
      urea <= 7.1 ~ '≤7.1',
      urea > 7.1 ~ '>7.1',
      T ~ NA_character_
    ) %>% factor(levels = c('≤7.1', '>7.1')),
    cd4_zhao = case_when(
      cd4 > 50 ~ '>50',
      cd4 <= 50 ~ '≤50',
      T ~ NA_character_
    ) %>% factor(levels = c('>50', '≤50'))
  ) %>%
  mutate(
    age_score_z = case_when(
      age_zhao == '≤49' ~ 0.0,
      age_zhao == '50-59' ~ 1.5,
      age_zhao == '60-69' ~ 3.0,
      age_zhao == '≥70' ~ 4.5
    ),
    gcs_score_z = case_when(
      gcs_zhao == 'No' ~ 0.0,
      gcs_zhao == 'Yes' ~ 2.5
    ),
    neck_stiff_score_z = case_when(
      neck_stiff == 'No' ~ 0,
      neck_stiff == 'Yes' ~ 2.0
    ),
    csf_op_score_z = case_when(
      csf_op_zhao_cat == '<200' ~ 0.0,
      csf_op_zhao_cat == '200-250' ~ 1.0,
      csf_op_zhao_cat == '251-350' ~ 2.0,
      csf_op_zhao_cat == '>350' ~ 3.0
    ),
    urea_score_z = case_when(
      urea_zhao == '≤7.1' ~ 0,
      urea_zhao == '>7.1' ~ 1.0
    ),
    cd4_score_z = case_when(
      cd4_zhao == '>50' ~ 0.0,
      cd4_zhao == '≤50' ~ 2.0
    )
  ) %>%
  mutate(
    zhao_score = case_when(
      !is.na(age_zhao) & !is.na(gcs_zhao) & !is.na(csf_op_zhao_cat) & !is.na(neck_stiff) & !is.na(urea_zhao) & !is.na(cd4_zhao) ~
        (age_score_z + gcs_score_z + neck_stiff_score_z + csf_op_score_z + urea_score_z + cd4_score_z),
      T ~ NA_real_
    )
  ) 

# Save coefficients
age_coef_z <- 1.11
gcs_coef_z <- 1.841
neck_stiff_coef_z <- 1.502
csf_op_coef_z <- 0.788
urea_coef_z <- 0.904
cd4_coef_z <- 1.505


# multiply these by score in Table 3 Zhao et al to get a offset value

all_mi <- all_mi %>%
  mutate(
    age_coef_score = case_when(
      age <= 49 ~ 0,
      age >= 50 & age <= 59 ~ 1,
      age >= 60 & age <= 69 ~ 2,
      age >= 70 ~ 3,
      T ~ NA_real_
    ),
    gcs_coef_score = case_when(
      gcs_bin == '15' ~ 0,
      gcs_bin != '15' ~ 1,
      T ~ NA_real_
    ),
    neck_stiff_coef_score = case_when(
      neck_stiff == 'No' ~ 0,
      neck_stiff == 'Yes' ~ 1,
      T ~ NA_real_
    ),
    csf_op_coef_score = case_when(
      csf_OP < 20 ~ 0,
      csf_OP >= 20 & csf_OP <= 25 ~ 1,
      csf_OP > 25 & csf_OP <= 35 ~ 2,
      csf_OP > 35 ~ 3,
      T ~ NA_real_
    ),
    urea_coef_score = case_when(
      urea_zhao == '≤7.1' ~ 0,
      urea_zhao == '>7.1' ~ 1,
      T ~ NA_real_
    ),
    cd4_coef_score = case_when(
      cd4_zhao == '>50' ~ 0,
      cd4_zhao == '≤50' ~ 1,
      T ~ NA_real_
    )
  ) %>%
  mutate(
    total_zhao_coef_score = ((age_coef_z * age_coef_score) + 
                               (gcs_coef_z * gcs_coef_score) + 
                               (neck_stiff_coef_z * neck_stiff_coef_score) +
                               (csf_op_coef_z * csf_op_coef_score) +
                               (urea_coef_z * urea_coef_score) +
                               (cd4_coef_z * cd4_coef_score))
  )





######################## External validation ##########################


# Set up validation cohort
val_mi <- all_mi %>% filter(country == "Malawi" & trial == "Ambition")
val_mi$death_28d <- as.numeric(val_mi$death_28d) - 1
val_mi$death_2week <- as.numeric(val_mi$death_2week) - 1

## Generate a coefficient for Zhao model

zhao_intercept_model <- glm(death_28d ~ 1, offset = total_zhao_coef_score, data = val_mi, family = 'binomial')
zhao_intercept <- zhao_intercept_model$coefficients[1]

#Create model predictions

val_mi <- val_mi %>%
  mutate(
    zhao_probability = inv.logit(
      zhao_intercept +
        (age_coef_z * age_coef_score) +
        (gcs_coef_z * gcs_coef_score) + 
        (neck_stiff_coef_z * neck_stiff_coef_score) +
        (csf_op_coef_z * csf_op_coef_score) +
        (urea_coef_z * urea_coef_score) +
        (cd4_coef_z * cd4_coef_score)
    )
  )

# Validation summary function, adapted from 08

get_validation_summary_point <- function(dataset, impno, outcome_name, method){
  dataset1 <- dataset %>% filter(.imp==impno)
  
  if (outcome_name == '28-day mortality'){
    dataset1$outcome <- dataset1$death_28d
    dataset1$outcome_name <- paste0('28-day mortality')
  }
  else if (outcome_name == '2-week mortality'){
    dataset1$outcome <- dataset1$death_2week
    dataset1$outcome_name <- paste0('2-week mortality')
  }
  else{
    stop('Invalid outcome choice. Please enter one of "28-day mortality" or "2-week mortality"')
  }
  
  if (method == 'point score'){
    dataset1$score <- dataset1$zhao_score
    dataset1$method <- paste0('point_score')
  }
  else if (method == 'probability'){
    dataset1$score <- dataset1$zhao_probability
    dataset1$method <- paste0('probability')
  }
  
  
  
  dataset1 <- dataset1 %>% select(study_id, score, outcome, outcome_name, method)
  dataset1 <- dataset1 %>% filter(!is.na(score) & !is.na(outcome))
  
  ci.auc <- roc(dataset1$outcome, dataset1$score, ci = T)$ci[1:3]
  
  dataset1$auc <- ci.auc[2]
  dataset1$auc.se <- sqrt(var(roc(dataset1$outcome ~ dataset1$score)))
  dataset1$auc.lower <- ci.auc[1]
  dataset1$auc.upper <- ci.auc[3]
  dataset1$auc.ci <- paste0(round(dataset1$auc,2), " (", round(dataset1$auc.lower,2), " - ", round(dataset1$auc.upper, 2), ")")
  
  dataset1$score_name <- paste0('Zhao')
  dataset1$imp <- paste(impno)
  dataset1$O <- sum(dataset1$outcome == 1)
  dataset1$N <- dim(dataset1)[1]
  
  if (method == 'probability'){
    dataset1$E <- sum(dataset1$score)
    dataset1$logit <- log(dataset1$score/(1 - dataset1$score))
    
    fit.val.slope <- glm(outcome ~ logit, data=dataset1, family=binomial)
    dataset1$slope<- fit.val.slope$coef[2]
    dataset1$slope.se <- sqrt(diag(vcov(fit.val.slope)))[2]
    
    fit.val.citl <- glm(outcome ~ 1, offset=logit, data=dataset1, family="binomial")
    dataset1$citl <- fit.val.citl$coef[1]
    dataset1$citl.se <- sqrt(diag(vcov(fit.val.citl)))[1]
    
    return(dataset1) 
  }
  else if (method == 'point score'){
    return(dataset1)
  }
  
}

validation_master_zhao_point <- data.frame()
validation_master_zhao_prob <- data.frame()
validation_master_zhao_2w_point <- data.frame()
validation_master_zhao_2w_prob <- data.frame()


######### Validation summaries

# 28d mortality, point score

for (i in c(1:10)) {
  val1z_point <- get_validation_summary_point(dataset = val_mi, impno = i, outcome = '28-day mortality',
                                              method = 'point score')
  validation_master_zhao_point <- bind_rows(validation_master_zhao_point, val1z_point)
}

val_master_summary_z_point <- validation_master_zhao_point %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, score_name, outcome_name, O, N)


# 28d mortality, probability score

for (i in c(1:10)) {
  val1z_prob <- get_validation_summary_point(dataset = val_mi, impno = i, outcome = '28-day mortality',
                                             method = 'probability')
  validation_master_zhao_prob <- bind_rows(validation_master_zhao_prob, val1z_prob)
}

val_master_summary_z_prob <- validation_master_zhao_prob %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, score_name, outcome_name, O, E, N)


# 2-week mortality, point score

for (i in c(1:10)) {
  val2z_point <- get_validation_summary_point(dataset = val_mi, impno = i, outcome = '2-week mortality',
                                              method = 'point score')
  validation_master_zhao_2w_point <- bind_rows(validation_master_zhao_2w_point, val2z_point)
}

val_master_summary_z_2w_point <- validation_master_zhao_2w_point %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, score_name, outcome_name, O, N)

# 2-week mortality, Probability score
for (i in c(1:10)) {
  val2z_prob <- get_validation_summary_point(dataset = val_mi, impno = i, outcome = '2-week mortality',
                                             method = 'probability')
  validation_master_zhao_2w_prob <- bind_rows(validation_master_zhao_2w_prob, val2z_prob)
}

val_master_summary_z_2w_prob <- validation_master_zhao_2w_prob %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, score_name, outcome_name, O, E, N)




# Validation pooling for points score
pool_validation_point <- function(dataset, method){
  dataset1 <- dataset
  
  if (method == 'point score'){
    q1 <- dataset1 %>% select(auc)
    se1 <- dataset1 %>% select(auc.se) 
  }
  else if (method == 'probability'){
    q1 <- dataset1 %>% select(auc, slope, citl)
    se1 <- dataset1 %>% select(auc.se, slope.se, citl.se)
  }
  pooled <- mi.meld(q1, se1, byrow = T)
  pooled.q <- as.data.frame(pooled$q.mi)
  pooled.se <- as.data.frame(pooled$se.mi)
  val_master_pooled1 <- cbind(pooled.q, pooled.se)
  
  val_master_pooled1$O <- round(mean(dataset1$O), 0)
  val_master_pooled1$N <- dataset1$N[1]
  val_master_pooled1$score_name <- paste0('Zhao')
  val_master_pooled1$outcome_name <- unique(dataset1$outcome_name)
  
  if (method == 'point score'){
    val_master_pooled1$method <- paste0('Point Score')
    val_master_pooled1 <- val_master_pooled1 %>% select(score_name, outcome_name, method, auc, auc.se, O, N)
  }
  else if (method == 'probability'){
    val_master_pooled1$E <- mean(dataset1$E) 
    val_master_pooled1$method <- paste0('Probability')
    val_master_pooled1 <- val_master_pooled1 %>% select(score_name, outcome_name, method,
                                                        auc, auc.se, slope, slope.se, 
                                                        citl, citl.se, O, E, N)
  }  
  
  return(val_master_pooled1)
}

# Run the pool validation function for all outcomes and methods
validation_summary_pooled_z <- pool_validation_point(val_master_summary_z_point, method = 'point score')
validation_summary_pooled_z_2w <- pool_validation_point(val_master_summary_z_2w_point, method = 'point score')
validation_summary_pooled_z_prob <- pool_validation_point(val_master_summary_z_prob, method = 'probability')
validation_summary_pooled_z_2w_prob <- pool_validation_point(val_master_summary_z_2w_prob, method = 'probability')


# Function to create 95%CI

get_confidence <- function(dataset, method){
  dataset1<- dataset
  dataset1$auc.lower <- dataset1$auc - 1.96*(dataset1$auc.se)
  dataset1$auc.upper <- dataset1$auc + 1.96*(dataset1$auc.se)
  
  dataset1$auc.ci <- paste0(round(dataset1$auc,2), " (", round(dataset1$auc.lower,2), " - ",
                            round(dataset1$auc.upper, 2), ")")
  
  if (method == 'point score'){
    return(dataset1)
  }
  else if (method == 'probability'){
    dataset1$slope.lower <- dataset1$slope - 1.96*(dataset1$slope.se)
    dataset1$slope.upper <- dataset1$slope + 1.96*(dataset1$slope.se)
    dataset1$citl.lower <- dataset1$citl - 1.96*(dataset1$citl.se)
    dataset1$citl.upper <- dataset1$citl + 1.96*(dataset1$citl.se)
    
    dataset1$slope.ci <- paste0(round(dataset1$slope,2), " (", round(dataset1$slope.lower,2), " - ",
                                round(dataset1$slope.upper, 2), ")")
    dataset1$citl.ci <- paste0(round(dataset1$citl,2), " (", round(dataset1$citl.lower,2), " - ",
                               round(dataset1$citl.upper, 2), ")")
    return(dataset1)
  }
}

# Generate CIs
validation_summary_pooled_z_point <- get_confidence(validation_summary_pooled_z, method = 'point score')
validation_summary_pooled_z_2w_point <- get_confidence(validation_summary_pooled_z_2w,  method = 'point score')
validation_summary_pooled_z_prob <- get_confidence(validation_summary_pooled_z_prob, method = 'probability')
validation_summary_pooled_z_2w_prob <- get_confidence(validation_summary_pooled_z_2w_prob, 'probability')


# Format to add 28d mortality results to other validation results
zhao_auc <- validation_summary_pooled_z_point$auc.ci

validation_results_zhao <- validation_summary_pooled_z_prob %>%
  select(score_name, auc.ci, slope.ci, citl.ci) %>%
  rename(Model = score_name, 
         AUROC = auc.ci,
         Slope = slope.ci,
         CITL = citl.ci) %>%
  mutate(
    AUROC = zhao_auc,
    CITL = '-'
  )
validation_results_zhao


######################### Calibration  ####################################

#Pooled calibration plot - Zhao

calibration_val_z <- ggplot(NULL) +
  geom_smooth(data=validation_master_zhao_prob[validation_master_zhao_prob$imp>0,], 
              aes(score, outcome),
              color="#4D8FA6",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_zhao_prob[validation_master_zhao_prob$imp==1,], 
           aes(score), sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_z


################## Decision Curve Analysis ###############################

# All for 2-week mortality

# Read in Basic and Research Validation Master dataframes

validation_master_r <- readRDS('validation_master_dataframe_researchmodel.rds')
validation_master_b <- readRDS('validation_master_dataframe_basicmodel.rds')
validation_master_e <- readRDS('validation_master_dataframe_ecogmodel.rds')
validation_master_z_dca <- validation_master_zhao_2w_prob

validation_master_r_dca <- validation_master_r %>%
  mutate(score_name = 'Research',
         outcome_name = '2-week mortality') 

validation_master_b_dca <- validation_master_b %>%
  mutate(score_name = 'Basic',
         outcome_name = '2-week mortality') 

validation_master_e_dca <- validation_master_e %>%
  mutate(score_name = 'ECOG',
         outcome_name = '2-week mortality') 

### DCA function 

get_dca_val_combined <- function(dataset, imp, chosen_model){
  dataset1 <- dataset %>% filter(imp==imp)
  dataset1 <- dataset1 %>% select(score, outcome)
  
  nb.opt.out.score <- decision_curve(outcome ~ score,
                                     data = dataset1, 
                                     policy = 'opt-in',
                                     thresholds = seq(0, .5, by = 0.01),
                                     bootstraps = 1, 
                                     fitted.risk = T)
  
  nb <- nb.opt.out.score$derived.data
  
  nb$score_name <- ifelse(nb$model != "All" & nb$model != "None", chosen_model, nb$model)
  
  nb$imp <- imp
  return(nb)
}

dca_master_val_combined <- data.frame()

for (i in c(1:10)) {
  nb1 <- get_dca_val_combined(dataset = validation_master_r_dca, imp = i, chosen_model = 'Research')
  dca_master_val_combined <- bind_rows(nb1, dca_master_val_combined)
}
for (i in c(1:10)) {
  nb1 <- get_dca_val_combined(dataset = validation_master_b_dca, imp = i, chosen_model = 'Basic')
  dca_master_val_combined <- bind_rows(nb1, dca_master_val_combined)
}
for (i in c(1:10)) {
  nb1 <- get_dca_val_combined(dataset = validation_master_z_dca, imp = i, chosen_model = 'Zhao')
  dca_master_val_combined <- bind_rows(nb1, dca_master_val_combined)
}
for (i in c(1:10)) {
  nb1 <- get_dca_val_combined(dataset = validation_master_e_dca, imp = i, chosen_model = 'ECOG')
  dca_master_val_combined <- bind_rows(nb1, dca_master_val_combined)
}

dca_master_val_combined$score_name <- factor(dca_master_val_combined$score_name, 
                                             levels = c(
                                               'All', 'None', 'Basic', 'Research', 'Zhao', 'ECOG'
                                             ))

linetype_mapping <- c("None" = "solid", "All" = "solid", "Basic" = "dotdash", "Research" = "longdash", "Zhao" = "dashed", "ECOG" = "twodash")
color_mapping <- c("None" = "red", "All" = "black", "Basic" = "#66A61E", "Research" = "#E6AB02", "Zhao" = "#4D8FA6", "ECOG" = '#B399D4')
linewidth_mapping <- c("None" = 0.6, "All" = 0.6, "Basic" = 0.7, "Research" = 0.7, "Zhao" = 0.7, "ECOG" = 0.7)

dca_plot_val_combined <- ggplot(dca_master_val_combined, 
                                aes(x=thresholds, y=NB,
                                    color = score_name, 
                                    linetype = score_name,
                                    linewidth = score_name
                                )) + 
  geom_smooth(aes(color = score_name),
              method="loess", se=F) +
  xlab("Threshold Probability") + ylab("Net Benefit") + 
  coord_cartesian(ylim=c(0, 0.1), xlim = c(0,0.4)) +
  theme_pubr(legend = "right") +
  theme(legend.title = element_blank()) +
  scale_linetype_manual(values = linetype_mapping) +
  scale_color_manual(values = color_mapping) +
  scale_linewidth_manual(values = linewidth_mapping)

dca_plot_val_combined

