# Cryptococcal meningitis prognostic model - Periskope-CM
# 05: Internal-external cross validation
# Description: Performs IECV using MI datasets by country in development cohort
# Started: 26/10/23
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

#Read in data
all_mi_aregi <- read_rds("complete_mi_aregimpute_131123.rds")

# Split into dev and val
dev_mi_aregi <- all_mi_aregi %>% filter(trial == "ACTA" | (trial == "Ambition" & country != "Malawi"))
val_mi_aregi <- all_mi_aregi %>% filter(country == "Malawi" & trial == "Ambition")


# Chosen Basic Model for IECV
f_basic_chosen <- death_2week ~
  gcs_bin + ecog + treatment +
  rcs(neut, 3) + rcs(haem_gl, 3)

# Chosen Research Model for IECV
f_research_chosen <- death_2week ~
  gcs_bin + ecog + treatment +
  rcs(neut, 3) + rcs(haem_gl, 3) + 
  rcs(csf_OP, 3) + 
  rcs(csf_qculture_log, 3)


# Remove labelled from df to allow functions to work and save samplesize

dev_mi_aregi <- data.frame(dev_mi_aregi)
sampsize_dev <- as.numeric(length(dev_mi_aregi$subjid[dev_mi_aregi$.imp==1])) 


########################## IECV Function ##############################


### Define ICEV calculation function (agnostic of basic v research)

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


################################ BASIC MODEL ########################


# Run the IECV calculation function for the basic model

iecv_master_b <- data.frame()

for (i in c(1:10)) {
  for (j in c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda")) {
    iecv1b <- get_iecv(dataset = dev_mi_aregi, f.chosen = f_basic_chosen, impno = i, location = j)
    iecv_master_b <- bind_rows(iecv_master_b, iecv1b)
  }
}

#Concise summary of IECV metrics
iecv_master_summary_b <- iecv_master_b %>%
  filter(.imp > 0) %>%
  distinct(country_iecv2, .imp, .keep_all = T) %>%
  select(country_iecv2, cstat, cstat.se, calslope, calslope.se, citl, citl.se, recal.int, O, E, N)


#Run pooled MI estimates function

iecv_master_pooled_b <- data.frame()

for (i in rev(c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda"))) {
  iecv_master_pooled1b <- pool_iecv_values(dataset = iecv_master_summary_b, location = i)
  iecv_master_pooled_b <- rbind(iecv_master_pooled_b, iecv_master_pooled1b)
}



########################## RESEARCH MODEL ############################

# Run the IECV calculation function for the research model

iecv_master_r <- data.frame()

for (i in c(1:10)) {
  for (j in c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda")) {
    iecv1r <- get_iecv(dataset = dev_mi_aregi, f.chosen = f_research_chosen, impno = i, location = j)
    iecv_master_r <- bind_rows(iecv_master_r, iecv1r)
  }
}

#Concise summary of IECV metrics
iecv_master_summary_r <- iecv_master_r %>%
  filter(.imp > 0) %>%
  distinct(country_iecv2, .imp, .keep_all = T) %>%
  select(country_iecv2, cstat, cstat.se, calslope, calslope.se, citl, citl.se, recal.int, O, E, N)


#Run pooled MI estimates function

iecv_master_pooled_r <- data.frame()

for (i in rev(c("Bots/SAfr", "Cameroon", "Malawi", "Tan/Zam/Zim", "Uganda"))) {
  iecv_master_pooled1r <- pool_iecv_values(dataset = iecv_master_summary_r, location = i)
  iecv_master_pooled_r <- rbind(iecv_master_pooled_r, iecv_master_pooled1r)
}


################## RANDOM EFFECTS META-ANALYSIS ##################



#################### IECV FOREST PLOT - BASIC MODEL #####################

pdf("Graphs and Tables/iecv_metrics_basicmodel_v3.pdf", height=6, width=16)

#CSTAT
par(fig=c(0,0.4,0,1))
model.cstat.b <- rma(yi = cstat, sei = cstat.se,  method = "REML", test = "knha", data = iecv_master_pooled_b)
summary(model.cstat.b)
metafor::forest(model.cstat.b,
                slab=iecv_master_pooled_b$location, xlab = "C-statistic", 
                alim=c(0.4, 1.0), cex=1, cex.lab=1, top=0,mlab = "Meta-analysis",
                xlim=c(0.23, 1.22))

## CITL

par(fig=c(0.4,0.7, 0,1), new=T)
model.CITL.b <- rma(yi = citl, sei = citl.se, method = "REML",test = "knha", data = iecv_master_pooled_b)
summary(model.CITL.b)
metafor::forest(model.CITL.b, slab=NA, xlab = "Calibration in the large", alim=c(-1.0, 1.0), 
                refline=0, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(-1.0, 2.1), at=c(-1, -0.5, 0, 0.5, 1))
## Slope

par(fig=c(0.7,1, 0,1), new=T)
model.CS.b <- rma(yi = calslope, sei = calslope.se, method = "REML", test = "knha", data = iecv_master_pooled_b)
summary(model.CS.b)
metafor::forest(model.CS.b, slab=NA, xlab = "Calibration slope", alim=c(0.0, 2.0), 
                refline=1, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(0.0, 3.0))

title("Meta-Analysis of IECV Metrics - Basic Model", line = -1.3, outer = T, cex.main = 1.2)

dev.off()

# Save the heterogeneity figures for Rmd

basic_cstat_i2 <- model.cstat.b$I2 %>% round(1)
basic_citl_i2 <- model.CITL.b$I2 %>% round(1)
basic_slope_i2 <- model.CS.b$I2 %>% round(1)

# Clear plot matrix
plot.new()
par(mfrow = c(1, 1))



####################### IECV FOREST PLOT - RESEARCH MODEL ##################

pdf("Graphs and Tables/iecv_metrics_researchmodel_v3.pdf", height=6, width=16)

#CSTAT
par(fig=c(0,0.4,0,1))
model.cstat.r <- rma(yi = cstat, sei = cstat.se,  method = "REML", test = "knha", data = iecv_master_pooled_r)
summary(model.cstat.r)
metafor::forest(model.cstat.r,
                slab=iecv_master_pooled_r$location, xlab = "C-statistic", 
                alim=c(0.5, 1.0), cex=1, cex.lab=1, top=0,mlab = "Meta-analysis",
                xlim=c(0.3, 1.22), at=c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

## CITL

par(fig=c(0.4,0.7, 0,1), new=T)
model.CITL.r <- rma(yi = citl, sei = citl.se, method = "REML",test = "knha", data = iecv_master_pooled_r)
summary(model.CITL.r)
metafor::forest(model.CITL.r, slab=NA, xlab = "Calibration in the large", alim=c(-1.0, 1.0), 
                refline=0, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(-1.0, 2.1), at=c(-1, -0.5, 0, 0.5, 1))
## Slope

par(fig=c(0.7,1, 0,1), new=T)
model.CS.r <- rma(yi = calslope, sei = calslope.se, method = "REML", test = "knha", data = iecv_master_pooled_r)
summary(model.CS.r)
metafor::forest(model.CS.r, slab=NA, xlab = "Calibration slope", alim=c(0.0, 2.0), 
                refline=1, cex=1, cex.lab=1, top=0, mlab = "",
                xlim=c(0.0, 3.0))

title("Meta-Analysis of IECV Metrics - Research Model", line = -1.3, outer = T, cex.main = 1.2)

dev.off()

plot.new()
par(mfrow = c(1, 1))

# Save measures of heterogeneity
research_cstat_i2 <- model.cstat.r$I2 %>% round(1)
research_citl_i2 <- model.CITL.r$I2 %>% round(1)
research_slope_i2 <- model.CS.r$I2 %>% round(1)



############################# Calibration plots for IECV ##########################

#Research
calibration_iecv_r <- ggplot(NULL) +
  geom_smooth(data=iecv_master_r[iecv_master_r$.imp>0,], 
              aes(predy.val, death_2week),
              color = 'red',
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_r[iecv_master_r$.imp==1,], aes(predy.val), sides="b", linewidth =0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) +
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_r


# Basic
calibration_iecv_b <- ggplot(NULL) +
  geom_smooth(data=iecv_master_b[iecv_master_b$.imp>0,], 
              aes(predy.val, death_2week),
              color = 'red',
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_b[iecv_master_b$.imp==1,], aes(predy.val), sides="b", linewidth =0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_b


###################### Calibration plots with recalibration of intercept ##################

# Research
calibration_iecv_recal_r <- ggplot(NULL) +
  geom_smooth(data=iecv_master_r[iecv_master_r$.imp>0,], 
              aes(recal.predy.val, death_2week),
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_r[iecv_master_r$.imp==1,], aes(recal.predy.val), sides="b", linewidth=0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_recal_r


# Basic
calibration_iecv_recal_b <- ggplot(NULL) +
  geom_smooth(data=iecv_master_b[iecv_master_b$.imp>0,], 
              aes(recal.predy.val, death_2week),
              method="loess", se=F,
              method.args=loess.control(statistics = "approximate", trace.hat = "approximate")) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
  geom_rug(data=iecv_master_b[iecv_master_b$.imp==1,], aes(recal.predy.val), sides="b", linewidth=0.1, alpha=0.1) +
  geom_abline(color="grey50") + 
  facet_wrap(~country_iecv2) + 
  theme_pubr() +
  xlab("Predicted risk") +
  ylab("Observed risk")

calibration_iecv_recal_b



## Density plot for distribution of predictions to illustrate.
research_model_1 <- lrm(f_research_chosen, data = dev_mi_aregi)

dev_mi_aregi <- data.frame(dev_mi_aregi)
dev_mi_aregi$research_predictions <- predict(research_model_1, type = 'fitted', newdata = dev_mi_aregi)

density_predict_dev <- dev_mi_aregi %>%
  filter(.imp >0)

iecv_prediction_density <- ggplot(data = density_predict_dev, aes(x = research_predictions, fill = as.factor(.imp))) +
  geom_density(alpha = 0.4) + 
  labs(title = "Prediction distribution in Development data by MI dataset",
       x = "Predicted probability",
       y = "Density") +
  scale_fill_discrete(name = "MI dataset")
iecv_prediction_density


