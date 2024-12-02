# Cryptococcal meningitis prognostic model - Periskope-CM
# 03: Multiple Imputation for CM model using aRegImpute
# Description: Uses aRegImpute to create 10 MI datasets
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

##Read in data
cm_main <- read_rds("cm_main_2023_10_02.rds")

# Prep
mi_input_all_mri <- cm_main %>%
  select(study_id, 
         age, sex, weight, seizure, gcs_bin, ecog, 
         novision, headache_dur, temp, resp_rate, 
         tb_hist,
         wcc, neut, haem_gl, cd4,
         csf_OP, csf_cellcount, csf_gluc, csf_prot, csf_qculture_log, 
         death_2week, death_10week, death_28d,
         treatment, trial_arm,
         papilloedema, on_arvs,
         urea, creatinine, neck_stiff,
         country, country_iecv, country_iecv2, trial,
         time_days_numeric
  )



# Impute all data combined

mi_input_complete <- mi_input_all_mri %>% select(-study_id, -country, -country_iecv, -country_iecv2, -trial, -trial_arm, -time_days_numeric)
mi_input_complete <- data.frame(mi_input_complete)

set.seed(19790)
mi_output_complete <- aregImpute(~ age + sex + weight + seizure + gcs_bin + ecog +
                                   novision + headache_dur + temp + resp_rate + 
                                   tb_hist + 
                                   wcc + neut + haem_gl + cd4 + 
                                   papilloedema + on_arvs +
                                   treatment +
                                   urea + creatinine + neck_stiff +
                                   csf_OP + csf_cellcount + csf_gluc + csf_prot + csf_qculture_log +
                                   death_2week + death_10week + death_28d,
                                 nk = 3, data = mi_input_complete, x = T, n.impute = 10)


####################### Stack MI datasets #####################


### Function to combine

fill_data_complete <- function(impute, data, im) {
  imp <- cbind.data.frame(impute.transcan(x = impute, imputation = im, data = data, list.out = T, pr = F))
  imp$imp <- im
  id <- data %>% select(study_id, country, country_iecv, country_iecv2, trial, trial_arm, time_days_numeric)
  imp <- cbind(id, imp)
  return(imp)
}

complete_mi <- mi_input_all_mri
complete_mi$imp <- 0

for (i in c(1:10)) {
  fulldat1 <- fill_data_complete(im = i, impute = mi_output_complete, data = mi_input_all_mri)
  complete_mi <- rbind(complete_mi, fulldat1)
}


complete_mi <- complete_mi %>% rename(.imp = imp)

# Save the output

saveRDS(complete_mi, file = "complete_mi_aregimpute_131123.rds")
