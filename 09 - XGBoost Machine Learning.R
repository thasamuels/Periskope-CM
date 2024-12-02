# Cryptococcal meningitis prognostic model - Periskope-CM
# 09: Machine learning using XGBoost
# Description: Compares performance of XGBoost model to (primary) logistic regression approach
# Started: 12/12/23
# Author: Rishi K Gupta, University College London

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
library(Hmisc)
library(corrplot)
library(xgboost)
library(caret)
library(gtools)
library(grid)
#library(recipes)

# Read in MI data, split into dev and val----

all_mi_aregi <- read_rds("complete_mi_aregimpute_131123.rds")
#all_mi_aregi$death_2week <- as.numeric(all_mi_aregi$death_2week)-1

dev_mi_aregi <- all_mi_aregi %>% 
  filter(trial == "ACTA" | (trial == "Ambition" & country != "Malawi")) %>%
  filter(.imp>0)

val_mi_aregi <- all_mi_aregi %>% 
  filter(country == "Malawi" & trial == "Ambition") %>%
  filter(.imp>0)

# First use default parameters----

## Default grid----

grid_default <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

## Default control----

control_default <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = F # FALSE for reproducible results 
)

## Train model----

xgb_base <- caret::train(
  death_2week ~ 
    gcs_bin + ecog+
    neut+ haem_gl+
    csf_OP+ csf_qculture_log+
    treatment,
  trControl = control_default,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE, 
  data=dev_mi_aregi %>%
    filter(.imp>0)
)

## Quick look at validation performance in stacked MI sets using default parameters----

val_mi_aregi$xgb_caret <- predict(xgb_base, val_mi_aregi, type="prob") %>%
  as_tibble %>%
  pull(2)

roc(val_mi_aregi$death_2week ~ val_mi_aregi$xgb_caret)

# Parameter tuning----

## Define grid search----

tune_grid <- expand.grid(
  nrounds = seq(from = 50, to = 150, by = 50),
  eta = c(0.1, 0.3, 0.5),
  max_depth = c(4, 6, 8),
  gamma = c(0,1),
  colsample_bytree = c(0.8, 1),
  min_child_weight = c(1, 2, 3),
  subsample = c(0.8, 1)
)

## Define control----

tune_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 10, # with n folds 
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE, # FALSE for reproducible results 
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

## Train model in 1 MI set as a test run----

# xgb_tune <- caret::train(
#   factor(death_2week) ~ 
#     gcs_bin + ecog+
#     neut+ haem_gl+
#     csf_OP+ csf_qculture_log+
#     treatment,
#   trControl = tune_control,
#   method = "xgbTree",
#   data=dev_mi_aregi %>%
#     filter(.imp==1),
#   tuneGrid = tune_grid,
#   verbose = TRUE, 
#   metric="ROC"
# )

# Function for tuning parameters in each MI set----

## Define the function----

tuning_mi <- function(data, imp){
  
  dataset1 <- data %>% filter(.imp == imp)
  
  set.seed(12345)
  xgb_tune1 <- caret::train(
    factor(death_2week) ~ 
      gcs_bin + ecog+
      neut+ haem_gl+
      csf_OP+ csf_qculture_log+
      treatment,
    trControl = tune_control,
    method = "xgbTree",
    data=dataset1,
    tuneGrid = tune_grid,
    verbose = F, 
    metric="ROC"
  )
  
  bestTune1 <-  xgb_tune1$bestTune
  bestTune1$n <- nrow(dataset1)
  bestTune1$imp <- imp
  return(bestTune1)
}

## Run the function (takes a while)---- hash out if want to avoid time retraining

tuning_mi_master <- data.frame()

 for (i in c(1:10)) {
  bestTune1 <- tuning_mi(imp = i, data = dev_mi_aregi)
  tuning_mi_master <- bind_rows(tuning_mi_master, bestTune1)
 }

sys.time()
tuning_mi_master
saveRDS(tuning_mi_master, "XGBoost data/tuning_mi_master.rds")
# Hash out to here

tuning_mi_master <- readRDS("XGBoost data/tuning_mi_master.rds")

## Find best tuning options (defined as the mode for each)----

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

grid_final <- tuning_mi_master %>%
  select(-imp, -n) %>%
  summarise(across(everything(), ~getmode(.)))

# saveRDS(grid_final, "XGBoost data/grid_final.rds")
grid_final <- readRDS("XGBoost data/grid_final.rds")

# Now train the model across the development MI sets, using the chosen parameters----

xgb_final <- caret::train(
  
  death_2week ~ 
    gcs_bin + ecog+
    neut+ haem_gl+
    csf_OP+ csf_qculture_log+
    treatment,
  
  trControl = control_default,
  tuneGrid = grid_final,
  method = "xgbTree",
  verbose = TRUE, 
  data=dev_mi_aregi)

 saveRDS(xgb_final, "XGBoost data/xgb_final.rds") # hash out if want to avoid time retraining
xgb_final <- readRDS("XGBoost data/xgb_final.rds")

### Extract the variable importance and create df

# Create baseline rows for categorical variables
baseline_categories <- data.frame(
  Variable = c('ECOG performance status', 'Treatment', 'GCS score'),
  Overall = c(rep('-', 3)),
  Level = c('Normal', '1wk AmBd+5FC', '15') 
)

xgb_variableimport <- varImp(xgb_final)$importance %>% 
  rownames_to_column(var = 'Variable') %>%
  mutate(Overall = round(Overall, 1)) %>% 
  mutate(
    Level = case_when(
      grepl("^ecog", Variable) ~ sub("^ecog", "", Variable),
      grepl("^gcs_bin", Variable) ~ sub("^gcs_bin", "", Variable),
      grepl("^treatment", Variable) ~ sub("^treatment", "", Variable),
      .default = ''
    )
  ) %>% 
  mutate(
    Variable = case_when(
      grepl("^ecog", Variable) ~ sub("ecog.*", "ECOG performance status", Variable),
      grepl("^gcs_bin", Variable) ~ sub("gcs_bin.*", "GCS score", Variable),
      grepl("^treatment", Variable) ~ sub("treatment.*", "Treatment", Variable),
      Variable == 'csf_qculture_log' ~ 'CSF QCC (log)',
      Variable == 'neut' ~ 'Neutrophils',
      Variable == 'csf_OP' ~ 'CSF Opening Pressure',
      Variable == 'haem_gl' ~ 'Haemoglobin',
      .default = Variable
    ),
    Level = ifelse(Level == 'FLU+5FC PO', 'FLU+5FC', Level)
  ) %>%
  arrange(desc(Overall)) %>%
  mutate(Overall = as.character(Overall)) %>% 
  add_row(baseline_categories) %>% 
  select(Variable, Level, Overall) %>% 
  rename('Variable Importance' = Overall) %>% 
  arrange(Variable)

xgb_variableimport

# IECV to assess generalisability of XGB approach----

## Define ICEV calculation function (adapted from 07)----

get_iecv_xgb <- function(dataset, impno, location){
  
  #Split into dev and val
  ipd.dev <- dataset %>% filter(.imp == impno & country_iecv2 != location)
  ipd.val <- dataset %>% filter(.imp == impno & country_iecv2 == location)
  
  #Train model in dev
  xgb_fit_dev <- caret::train(
    death_2week ~ 
      gcs_bin + ecog+
      neut+ haem_gl+
      csf_OP+ csf_qculture_log+
      treatment,
    trControl = control_default,
    tuneGrid = grid_final,
    method = "xgbTree",
    verbose = TRUE, 
    data=ipd.dev)
  
  #Predictions
  ipd.val$predy.val <- predict(xgb_fit_dev, ipd.val, type="prob") %>%
    as_tibble %>%
    pull(2)
  ipd.val$predlp.val <- logit(ipd.val$predy.val)
  
  #Determine calibration slope in validation data
  fit.val.slope <- glm(death_2week ~ predlp.val, data = ipd.val, family = 'binomial')
  
  #Determine calibration-in-the-large in the validation data
  fit.val.citl <- glm(death_2week ~ 1, offset = predlp.val, data = ipd.val, family = 'binomial')
  
  #Determine discriminative performance in validation data (c-statistic)
  cstat.val <- roc(ipd.val$death_2week, ipd.val$predy.val, ci = T)$ci[1:3]
  
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
  ipd.val$death_2week <- as.numeric(ipd.val$death_2week) - 1
  ipd.val$xgb <- ipd.val$predy.val
  
  ipd.val <- ipd.val %>% select(.imp, study_id,
                                country_iecv2,
                                cstat.lb, cstat, cstat.ub, cstat.se,
                                citl, citl.se,
                                calslope, calslope.se,
                                O, E, N,
                                predy.val, predlp.val,
                                xgb, 
                                death_2week)
  return(ipd.val)
}

## Run the IECV function for the XGB model----

iecv_master_x <- data.frame()

for (i in c(1:10)) {
  for (j in c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda")) {
    iecv1x <- get_iecv_xgb(dataset = dev_mi_aregi, impno = i, location = j)
    iecv_master_x <- bind_rows(iecv_master_x, iecv1x)
  }
}

## Concise summary of IECV metrics----

iecv_master_summary_x <- iecv_master_x %>%
  filter(.imp > 0) %>%
  distinct(country_iecv2, .imp, .keep_all = T) %>%
  select(country_iecv2, cstat, cstat.se, calslope, calslope.se, citl, citl.se, O, E, N)

## Pool MI estimates (function from 07)----

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

iecv_master_pooled_x <- data.frame()

for (i in rev(c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda"))) {
  iecv_master_pooled1x <- pool_iecv_values(dataset = iecv_master_summary_x, location = i)
  iecv_master_pooled_x <- rbind(iecv_master_pooled_x, iecv_master_pooled1x)
}

saveRDS(iecv_master_pooled_x, "XGBoost data/iecv_master_pooled_x.rds")

## IECV forest plots----

pdf("XGBoost data/iecv_metrics_xgb.pdf", height=6, width=16)

### CSTAT
par(fig=c(0,0.4,0,1))
model.cstat.x <- rma(yi = cstat, sei = cstat.se,  method = "REML", test = "knha", data = iecv_master_pooled_x)
summary(model.cstat.x)
metafor::forest(model.cstat.x,
                slab=iecv_master_pooled_x$location, xlab = "C-statistic", 
                alim=c(0.5, 1.0), cex=1, cex.lab=1, top=0,mlab = "Meta-analysis",
                xlim=c(0.3, 1.22), at=c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

### CITL

par(fig=c(0.4,0.7, 0,1), new=T)
model.CITL.x <- rma(yi = citl, sei = citl.se, method = "REML",test = "knha", data = iecv_master_pooled_x)
summary(model.CITL.x)
metafor::forest(model.CITL.x, slab=NA, xlab = "Calibration in the large", alim=c(-1.0, 1.0), 
                refline=0, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(-1.0, 2.1), at=c(-1, -0.5, 0, 0.5, 1))
### Slope

par(fig=c(0.7,1, 0,1), new=T)
model.CS.x <- rma(yi = calslope, sei = calslope.se, method = "REML", test = "knha", data = iecv_master_pooled_x)
summary(model.CS.x)
metafor::forest(model.CS.x, slab=NA, xlab = "Calibration slope", alim=c(0.0, 2.0), 
                refline=1, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(0.0, 3.0))

title("Meta-Analysis of IECV Metrics - XGBoost Model", line = -1.3, outer = T, cex.main = 1.2)

dev.off()

# Save heterogeneity
x_cstat_i2 <- model.cstat.x$I2 %>% round(1)
x_citl_i2 <- model.CITL.x$I2 %>% round(1)
x_slope_i2 <- model.CS.x$I2 %>% round(1)

## Calibration plots for IECV XGB----

calibration_iecv_x <- ggplot(NULL) +
  geom_smooth(data=iecv_master_x[iecv_master_x$.imp>0,], 
              aes(predy.val, death_2week),
              color = 'red',
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_x[iecv_master_x$.imp==1,], aes(predy.val), sides="b", linewidth =0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) +
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_x

# Validation set----

## Generate XGB and LRM predictions in validation set----

val_mi_aregi$xgb <- predict(xgb_final, val_mi_aregi, type="prob") %>%
  as_tibble %>%
  pull(2)

research_model <- readRDS("research_mi_model_251023.rds")
val_mi_aregi <- data.frame(val_mi_aregi)
val_mi_aregi$crypto_model_r <- predict(research_model, val_mi_aregi, type="fitted")

ggplot(val_mi_aregi) + 
  geom_point(aes(x=crypto_model_r, y = xgb)) +
  theme_pubr()

## Density plot----

prediction_density_x_val <- 
  val_mi_aregi %>%
  mutate(death_2week = factor(death_2week, levels = c('Died', 'Alive'))) %>% 
  filter(.imp>0) %>%
  ggplot(aes(x = xgb, fill = death_2week)) +
  geom_density(alpha = 0.6, linewidth = 0.5) +
  labs(x = "Predicted Probability of Mortality",
       y = "Density") +
  scale_fill_manual(values = c('#CC0000', '#FF6666'), name = "2-week mortality") + 
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 9)) +
  theme_pubr()

prediction_density_x_val
saveRDS(prediction_density_x_val, "XGBoost data/prediction_density_xgb.rds")

## Validate performance----

### Function (adapted from 08)----

get_validation_summary <- function(dataset, impno, chosen_model){
  dataset1 <- dataset %>% filter(.imp==impno)
  dataset1$outcome <- dataset1$death_2week
  
  if (chosen_model == "basic") {
    dataset1$score = dataset1$crypto_model_b
  }
  else if (chosen_model == 'research') {
    dataset1$score = dataset1$crypto_model_r
  }
  else if (chosen_model == 'xgb') {
    dataset1$score = dataset1$xgb
  }
  else{
    stop("Invalid model choice. Please use basic, research or xgb")
  }
  
  ci.auc <- roc(dataset1$outcome, dataset1$score, ci = T)$ci[1:3]
  dataset1$auc <- ci.auc[2]
  dataset1$auc.se <- sqrt(var(roc(dataset1$outcome ~ dataset1$score)))
  dataset1$auc.lower <- ci.auc[1]
  dataset1$auc.upper <- ci.auc[3]
  dataset1$auc.ci <- paste0(round(dataset1$auc,2), " (", round(dataset1$auc.lower,2), " - ", round(dataset1$auc.upper, 2), ")")
  dataset1$imp <- paste(impno)
  dataset1$O <- sum(dataset1$outcome == "Died")
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

### Run function----

validation_master_x <- data.frame()

for (i in c(1:10)) {
  val1x <- get_validation_summary(dataset = val_mi_aregi, chosen_model = "xgb", impno = i)
  validation_master_x <- bind_rows(validation_master_x, val1x)
}

val_master_summary_x <- validation_master_x %>%
  filter(imp>0) %>%
  distinct(imp, .keep_all = T) %>%
  select(imp, auc, auc.se, slope, slope.se, citl, citl.se, O, E, N)

val_master_summary_x

## Pool validation (function from 08)----

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

validation_summary_pooled_x <- pool_validation(val_master_summary_x)
validation_summary_pooled_x

## Function to create 95%CI (function from 08)----

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

validation_summary_pooled_x <- get_confidence(validation_summary_pooled_x)
validation_summary_pooled_x %>%
  select(N, O, E, auc.ci, slope.ci, citl.ci)
saveRDS(validation_summary_pooled_x, "XGBoost data/validation_summary_pooled_x.rds")

## Pooled calibration plot----

calibration_val_x <- ggplot(NULL) +
  geom_smooth(data=validation_master_x[validation_master_x$.imp>0,], 
              aes(score, as.numeric(outcome)-1),
              color="red",
              method="loess", se=T,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) + 
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=validation_master_x[validation_master_x$imp==1,], 
           aes(score), sides="b", linewidth=0.15, alpha=0.25) +
  geom_abline(color="grey50") + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk") 

calibration_val_x
saveRDS(calibration_val_x, "XGBoost data/calibration_val_x.rds")

# Examine associations in XGB model----

## Create a data frame of 100 records with 'base case' average values----

marginals_data <- dev_mi_aregi %>%
  select(gcs_bin, ecog,
         neut, haem_gl,
         csf_OP, csf_qculture_log,
         treatment)

## Define function to vary each variable across a specified range, including 2 x 2 interactions----

get_marginals <- function(dataset_marginals, xvar, yvar){
  dataset1 <- dataset_marginals
  
  if (is.numeric(dataset1[[xvar]])) {
    x_q5 <-  dataset_marginals[[xvar]] %>% quantile(0.05) # Distribution of continuous vars will be shown betweem 5th and 95th percentiles
    x_q95 <-  dataset_marginals[[xvar]] %>% quantile(0.95)
    xvalues <- seq(from = x_q5, to = x_q95, length.out = 100) # 100 values per continuous variable
  }
  
  if (is.factor(dataset1[[xvar]])) {
    xvalues <-  dataset_marginals[[xvar]] %>% levels # Each level will be represented for factor variables
  }
  
  dataset1 <- dataset1 %>%
    mutate(
      across(where(is.numeric), ~median(.x)), # use median for continuous vars as base case
      across(where(is.factor), ~getmode(.x))) %>%  # use mode for factor vars as base case
    dplyr::slice(1:length(xvalues))
  
  dataset1[[xvar]] <- xvalues
  
  if(xvar != yvar) {
    if (is.factor(dataset1[[yvar]])) { # y represents the potential interaction variable 
      yvalues <-  dataset_marginals[[yvar]] %>% levels() # For factor vars, y will show each level of the factor variable
    }
    
    if (is.numeric(dataset1[[yvar]])) {
      q12.5 <- dataset_marginals[[yvar]] %>% quantile(0.125) %>% round(., 1) # For continuous y values, plots will be produced
      q37.5 <- dataset_marginals[[yvar]] %>% quantile(0.375) %>% round(., 1) # for the midpoint of each quarter of the distribution
      q62.5 <- dataset_marginals[[yvar]] %>% quantile(0.625) %>% round(., 1)
      q87.5 <- dataset_marginals[[yvar]] %>% quantile(0.875) %>% round(., 1)
      yvalues <-  c(q12.5, q37.5, q62.5, q87.5)
    }
    
    dataset1  <- uncount(dataset1, length(yvalues), .id="id")
    dataset1[[yvar]] <- rep(yvalues, nrow(dataset1)/length(yvalues))
    
  }
  
  if(xvar == yvar){
    dataset1$id <- 1
  }
  
  dataset1$xvar <- xvar
  dataset1$x_var_type <- class(dataset1[[xvar]])
  dataset1$x <- as.character(dataset1[[xvar]])
  
  dataset1$yvar <- yvar
  dataset1$y_var_type <- class(dataset1[[yvar]])
  dataset1$y <- as.character(dataset1[[yvar]])
  
  marginals_master <<- rbind(dataset1, marginals_master)
}

## Run function----

marginals_master <- data.frame()

vars <- c("csf_qculture_log", "haem_gl", "neut","csf_OP", "ecog",  "gcs_bin", "treatment")

for (i in vars) {
  for (j in vars) {
    get_marginals(dataset_marginals = marginals_data, xvar=i, yvar=j)
  }
}

## Add XGB predictions----

marginals_master$xgb <- predict(xgb_final, marginals_master, type="prob") %>%
  as_tibble %>%
  pull(2)

## Add LRM predictions----

marginals_master <- data.frame(marginals_master)
marginals_master$lrm <- predict(research_model, marginals_master, type="fitted")

## Add treatment labels back in to data frame (lost somewhere)----

marginals_master <- marginals_master %>%
  mutate(treatment = 
           #factor(treatment, levels = c("Ambition_arms", "AmBd_Flu_1w", "AmBd_5FC_2w", "AmBd_Flu_2w", "Flu_5FC_PO")),
           as.numeric(factor(treatment, levels = c("1wk AmBd+5FC/Ambisome", "1wk AmBd+FLU", "2wks AmBd+5FC", "2wks AmBd+FLU", # coded as numbered factor to tidy up plot
                                                   "FLU+5FC PO"))),
         
         treatment = factor(treatment),
         
         ecog = 
           #factor(treatment, levels = c("Ambition_arms", "AmBd_Flu_1w", "AmBd_5FC_2w", "AmBd_Flu_2w", "Flu_5FC_PO")),
           factor(as.numeric(ecog)-1), 
         
         xlabel = case_when(
           xvar=="ecog" ~ "ECOG",
           xvar=="gcs_bin" ~ "GCS",
           xvar=="treatment" ~ "Treatment arm",
           xvar=="neut" ~ "Neutrophils (x10^9/L)",
           xvar=="haem_gl" ~ "Haemoglobin (g/L)",
           xvar=="csf_qculture_log" ~ "Quantitative culture (log cfu/mL)",
           xvar=="csf_OP" ~ "Opening pressure (cmH20)"), 
         
         ylabel = case_when(
           yvar=="ecog" ~ "ECOG",
           yvar=="gcs_bin" ~ "GCS",
           yvar=="treatment" ~ "Treatment arm",
           yvar=="neut" ~ "Neutrophils (x10^9/L)",
           yvar=="haem_gl" ~ "Haemoglobin (g/L)",
           yvar=="csf_qculture_log" ~ "Quantitative culture (log cfu/mL)",
           yvar=="csf_OP" ~ "Opening pressure (cmH20)"))

## Function to create marginal effect plots----

get_marg_plot <- function(xvariable, yvariable, dataset){
  
  dataset1 <- dataset %>%
    filter(xvar==xvariable, yvar==yvariable)
  
  dataset1$x <- dataset1[[xvariable]]
  dataset1$y <- dataset1[[yvariable]]
  
  if (xvariable == yvariable & is.numeric(dataset1$x)) {
    plot1 <- ggplot(dataset1) +
      aes(x = x, y = xgb) +
      geom_smooth(linewidth = 1, se=F, color="grey50") +
      theme_pubr(base_size = 10, legend="right") +
      xlab(ifelse(dataset1$yvar[1]=="treatment", dataset1$xlabel[1], "")) +
      theme(axis.title.y = element_blank())+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), breaks = scales::breaks_extended(n = 4))
    
  }
  
  else if (xvariable == yvariable & is.factor(dataset1$x)) {
    plot1 <- ggplot(dataset1) +
      aes(x = x, y = xgb) +
      geom_point(size = 2, color="grey50") +
      geom_line(aes(group=1), linetype="dashed",  color="grey50") +
      theme_pubr(base_size = 10, legend="right") +
      xlab(ifelse(dataset1$yvar[1]=="treatment", dataset1$xlabel[1], "")) +
      theme(axis.title.y = element_blank())+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), breaks = scales::breaks_extended(n = 4))
    
  }
  
  else  if(is.numeric(dataset1$x)) {
    plot1 <- ggplot(dataset1) +
      aes(x = x, y = xgb, color=factor(y)) +
      geom_smooth(linewidth = 1, se=F) +
      theme_pubr(base_size = 10, legend="right") +
      xlab(ifelse(dataset1$yvar[1]=="treatment", dataset1$xlabel[1], "")) +
      theme(axis.title.y = element_blank()) +
      labs(colour=dataset1$ylabel[1])+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), breaks = scales::breaks_extended(n = 4)) +
      scale_color_manual(values = wesanderson::wes_palette(name = "Darjeeling1"))
  }
  
  else if(is.factor(dataset1$x)) {
    plot1 <- ggplot(dataset1) +
      aes(x = x, y = xgb, color=factor(y)) +
      geom_line(aes(group=factor(y)), linetype="dashed") +
      geom_point(size = 2) +
      theme_pubr(base_size = 10, legend="right") +
      xlab(ifelse(dataset1$yvar[1]=="treatment", dataset1$xlabel[1], "")) +
      theme(axis.title.y = element_blank()) +
      labs(colour=dataset1$ylabel[1]) +
      #theme(legend.title.align=0.5, 
      #      legend.justification = "centre")  +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01), breaks = scales::breaks_extended(n = 4)) +
      scale_color_manual(values = wesanderson::wes_palette(name = "Darjeeling1"))
    #theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    #theme(axis.text.x=element_blank())
  }
  
  return(plot1)
  
}


## Create empty list and then run function to add all plots----

myplots <- list()

for (i in vars) {
  for (j in vars) {
    plot1 <- get_marg_plot(dataset=marginals_master, xvariable=i, yvariable=j)
    name <- paste0(i, "_", j)
    myplots[[name]] <- plot1
  }
}

## Combine plots into "rows" with common legends----

csf_qculture_plots <- ggarrange(plotlist  = myplots[seq(1, 49, by=7)], ncol=7, common.legend = T, legend="top")
haem_gl_plots <- ggarrange(plotlist  = myplots[seq(2, 49, by=7)], ncol=7, common.legend = T, legend="top")
neut_plots <- ggarrange(plotlist  = myplots[seq(3, 49, by=7)], ncol=7, common.legend = T, legend="top")
csf_op_plots <- ggarrange(plotlist  = myplots[seq(4, 49, by=7)], ncol=7, common.legend = T, legend="top")
ecog_plots <- ggarrange(plotlist  = myplots[seq(5, 49, by=7)], ncol=7, common.legend = T, legend="top")
gcs_plots <- ggarrange(plotlist  = myplots[seq(6, 49, by=7)], ncol=7, common.legend = T, legend="top")
treatment_plots <- ggarrange(plotlist  = myplots[seq(7, 49, by=7)], ncol=7, common.legend = T, legend="top")

## Combine rows into plot matrix----

marginals_plot_x <- ggarrange(csf_qculture_plots, haem_gl_plots, neut_plots, csf_op_plots, ecog_plots, gcs_plots,treatment_plots, ncol=1, align="hv")
marginals_plot_x <- annotate_figure(marginals_plot_x, left = textGrob("Predicted 2-week mortality", rot = 90, vjust = 1, gp = gpar(cex = 0.9)))
marginals_plot_x

#Save plot
saveRDS(marginals_plot_x, "XGBoost data/marginals_plot.rds")
