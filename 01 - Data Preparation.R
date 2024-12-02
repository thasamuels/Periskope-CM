# Cryptococcal meningitis prognostic model
# 01: Data Prep
# Description: Cleaning and Labeling of variables. Interrogation of pattern of missingness. 
# Started: 27/09/2023
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


#load master dataset

cm_main <- read_csv("master_acta_ambition_complete.csv")


####Convert dates to correct structures
cm_main <- cm_main %>%
  mutate(
    hosadmdate_scr = as.Date(hosadmdate_scr, format = "%d-%b-%y"),
    inclusiondate_scr = as.Date(inclusiondate_scr, format = "%d-%b-%y"),
    datedeath_dth = as.Date(datedeath_dth, format = "%d-%b-%y"),
    date_term_end = as.Date(date_term_end, format = "%d-%b-%y")
  ) 


##########Creating the primary outcome variable and secondary mortality endpoints:

cm_main <- cm_main %>%
  mutate(
    death_2week = case_when(died == 1 & study_day_death <= 15 ~ "Died", 
                            died == 1 & study_day_death >15 | died == 0 ~ "Alive",
                            T ~ NA_character_),
    death_10week = case_when(died == 1 & study_day_death <= 71 ~ "Died", 
                             died == 1 & study_day_death >71 | died == 0 ~ "Alive",
                             T ~ NA_character_)
  )

###### Code trial, country and trial arm.

cm_main <- cm_main %>%
  mutate(
    trial = factor(trial),
    country = factor(country,
                     levels = c("Botswana", "Cameroon", "Malawi", "South Africa", "Tanzania", "Uganda", "Zambia", "Zimbabwe")),
    trial_arm = ifelse(randassgmt_scr == "Oral - FLU_5FC", "FLU+5FC PO", randassgmt_scr)
  ) %>%
  mutate(
    trial_arm = factor(trial_arm,
                       levels = c("Ambisome", "1wk AmBd+5FC", "1wk AmBd+FLU", "2wks AmBd+5FC", "2wks AmBd+FLU", "FLU+5FC PO"))
  )

####### Clean and rename the predictor variables of interest

cm_main <- cm_main %>%
  mutate(
    age = case_when( #checks age, if >110 sets to NA, if not, remains age
      age_scr >110 ~ NA_real_,
      T ~ age_scr),
    sex= factor(sex_scr, levels = c("Female", "Male")),
    weight = ifelse(is.na(weight_pth), NA_real_, round(weight_pth)),
    seizure = factor(seizls72curr_pth),
    novision = factor(vsnlosscurr_pth),
    headache_dur = case_when( #there is a headache of duration 500 here - I have deleted. Needs a sense check
      hdachedur_pth >200 ~ NA_real_,
      T ~ hdachedur_pth),
    gcscore_adm = case_when(
      !is.na(gcscore_pth) & gcscore_pth >= 3 & gcscore_pth <= 15 & gcscore_pth %%1 == 0 ~ gcscore_pth,
      T ~ NA_integer_),
    ecog = case_when(
      !is.na(ecogscore_pth) & ecogscore_pth >= 0 & ecogscore_pth <= 4 & ecogscore_pth %%1 == 0 ~ as.character(ecogscore_pth),
      T ~ NA_character_) %>%
      recode_factor("0" = "Normal",
                    "1" = "Restricted activity", 
                    "2" = "Ambulatory",
                    "3" = "Limited self-care",
                    "4" = "Bedbound"),
    temp = ifelse(is.na(temp_pth), NA_real_, round(temp_pth, 1)),
    papill_pth = as.factor(papill_pth), #included while trying to resolve factor problem
    papilloedema = ifelse(papill_pth == "Unchecked", "No", "Yes") %>%
      factor(levels = c("No", "Yes")),
    resp_rate = case_when(
      rr_pth > 60 ~ NA_real_,
      T ~ rr_pth),
    on_arvs = case_when(
      !is.na(onarvs_pth) & onarvs_pth != "Blank" ~ onarvs_pth,
      T ~ NA_character_),
    tb_hist = case_when(!is.na(pmh_tb_pth) ~ pmh_tb_pth, T ~ NA_character_),
    #currtb doesnt need editing, but please note below - high degree of missingness
    wcc = ifelse(is.na(wbc_brf), NA_real_, round(wbc_brf, 2)),
    neut = ifelse(is.na(absneutrophil_brf), NA_real_, round(absneutrophil_brf, 2)),
    haem_gl = case_when(!is.na(hemoglobin_brf) ~ hemoglobin_brf * 10,
                        T ~ NA_real_),
    cd4 = ifelse(is.na(abscd4_brf), NA_real_, round(abscd4_brf)),
    csf_OP = case_when(!is.na(openpressure_csf) ~ openpressure_csf, T ~ NA_integer_),
    csf_cellcount = case_when(!is.na(totcellcount_csf) ~ totcellcount_csf, T ~ NA_integer_),
    csf_gluc = ifelse(is.na(glucose_mmol_csf), NA_real_, round(glucose_mmol_csf, 2)),
    csf_prot = ifelse(is.na(protein_csf), NA_real_, round(protein_csf, 2)),
    csf_qculture = case_when(!is.na(qculture_csf_0) ~ qculture_csf_0, T ~ NA_real_),
    hosadmdate_scr = case_when( #make sure missing dates marked as NA
      is.finite(hosadmdate_scr) ~ hosadmdate_scr, 
      T ~ as.Date(NA_character_)), 
    inclusiondate_scr = case_when(
      is.finite(inclusiondate_scr) ~ inclusiondate_scr, 
      T ~ as.Date(NA_character_)),
    datedeath_dth = case_when(
      is.finite(datedeath_dth) ~ datedeath_dth,
      T ~ as.Date(NA_character_)),
    date_term_end = case_when(
      is.finite(date_term_end) ~ date_term_end,
      T ~ as.Date(NA_character_))
    
  )

##### Clean and rename predictor variables for trial outcome modelling 
cm_main <- cm_main %>%
  mutate(
    haem_d7 = case_when(!is.na(haemoglobin_d7) ~ haemoglobin_d7 * 10, T ~ NA_real_),
    haem_d14 = case_when(!is.na(haemoglobin_d14) ~ haemoglobin_d14 * 10, T ~ NA_real_),
    creat = ifelse(is.na(creat_umol_brf), NA_real_, round(creat_umol_brf)),
    creat_d7 = ifelse(is.na(creatinine_d7), NA_real_, round(creatinine_d7)),
    creat_d14 = ifelse(is.na(creatinine_d14), NA_real_, round(creatinine_d14)),
    csf_qculture_d7 = case_when(!is.na(quantitative_culture_d7) ~ quantitative_culture_d7, T ~ NA_real_),
    csf_qculture_d14 = case_when(!is.na(quantitative_culture_d14) ~ quantitative_culture_d14, T ~ NA_real_),
    csf_op_d7 = case_when(!is.na(opening_pressure_d7) ~ opening_pressure_d7, T ~ NA_real_),
    csf_op_d14 = case_when(!is.na(opening_pressure_d14) ~ opening_pressure_d14, T ~ NA_real_),
  )

#Log1 transform qcultures. First make all 0 = 1 to stop -inf values, then log10 transform
cm_main <- cm_main %>%
  mutate(
    cssfqculture_nozero = ifelse(csf_qculture == 0, 1, csf_qculture),
    cssfqculture_nozero_d7 = ifelse(csf_qculture_d7 == 0, 1, csf_qculture_d7),
    cssfqculture_nozero_d14 = ifelse(csf_qculture_d14 == 0, 1, csf_qculture_d14)
  )

cm_main <- cm_main %>%
  mutate(
    csf_qculture_log = case_when(!is.na(cssfqculture_nozero) ~ log10(as.numeric(cssfqculture_nozero)), T ~ NA_real_),
    csf_qculture_d7_log = case_when(!is.na(cssfqculture_nozero_d7) ~ log10(as.numeric(cssfqculture_nozero_d7)), T ~ NA_real_),
    csf_qculture_d14_log = case_when(!is.na(cssfqculture_nozero_d14) ~ log10(as.numeric(cssfqculture_nozero_d14)), T ~ NA_real_),
  )

# Create a proxy categorical 3-bin GCS variable (decision to treat this way instead of continuous)

cm_main <- cm_main %>%
  mutate(
    gcs_bin = case_when(
      gcscore_adm >= 3 & gcscore_adm <= 10 ~ "<=10",
      gcscore_adm > 10 & gcscore_adm < 15 ~ "11-14",
      gcscore_adm == 15 ~ "15", 
      T ~ NA_character_
    ) %>%
      factor(levels = c("<=10", "11-14", "15"))
  )

# Create a dichotomous ECOG variable if need to reduce sample size footprint of ecog

cm_main <- cm_main %>%
  mutate(
    ecog_bin = case_when(
      ecogscore_pth >= 0 & ecogscore_pth <= 1 ~ "0-1",
      ecogscore_pth > 1 & ecogscore_pth <= 4 ~ "2+",
      T ~ NA_character_
    ) %>%
      factor(levels = c("0-1", "2+"))
  )

# Creating trial arm and IECV country factor variables.

cm_main <- cm_main %>%
  mutate(
    flucytosine = case_when(
      grepl('5FC', trial_arm) | trial_arm == "Ambisome" ~ "5FC",
      TRUE ~ 'No 5FC'
    ) %>%
      factor(levels = c("No 5FC", "5FC")),
    treatment = case_when( 
      trial_arm == "FLU+5FC PO" ~ "FLU+5FC PO",
      trial_arm == "2wks AmBd+FLU" ~ "2wks AmBd+FLU",
      trial_arm == "2wks AmBd+5FC" ~ "2wks AmBd+5FC",
      trial_arm == "1wk AmBd+FLU" ~ "1wk AmBd+FLU",
      trial_arm == "Ambisome" | trial_arm == "1wk AmBd+5FC" ~ "1wk AmBd+5FC/Ambisome",
    ) %>%
      factor(levels = c("1wk AmBd+5FC/Ambisome", "1wk AmBd+FLU", "2wks AmBd+5FC", "2wks AmBd+FLU", 
                        "FLU+5FC PO")),
    country_iecv = case_when(
      country == 'Botswana' | country == 'South Africa' | country == 'Tanzania' |
        country == 'Zambia' | country == 'Zimbabwe' ~ "Bots/SAf/Tanz/Zam/Zim",
      country == 'Cameroon' ~ "Cameroon",
      country == 'Malawi' ~ "Malawi",
      country == 'Uganda' ~ "Uganda"
    ),
    country_iecv2 = case_when(
      country == 'Botswana' | country == 'South Africa' ~ "Bots/SAfr",
      country == 'Zambia' | country == 'Zimbabwe' | country == 'Tanzania' ~ "Tan/Zam/Zim",
      country == 'Cameroon' ~ "Cameroon",
      country == 'Malawi' ~ "Malawi",
      country == 'Uganda' ~ "Uganda"
    )
  )

########## Clean extreme/incorrect values ############

# ID the relevant extreme values for CSF

csf_prot_filter <- cm_main %>% filter(csf_prot > 20) %>%
  select(study_id, csf_prot, trial) %>% print()

csf_gluc_filter <- cm_main %>% filter(csf_gluc > 10) %>%
  select(study_id, trial, csf_gluc) %>% print()

# Remove extreme values from relevant variables

cm_main <- cm_main %>%
  mutate(
    csf_gluc = round(ifelse(csf_gluc > 10, csf_gluc / 10, csf_gluc),2),
    csf_prot = ifelse(csf_prot > 20, NA_real_, csf_prot)
  )

##### MAKE SURE ALL FACTOR VARIABLES ARE TREATED AS SUCH


#To analyse which factor variables need correcting
filter4factor <- cm_main %>% 
  select(sex, seizure, novision, gcs_bin, ecog, ecog_bin, papilloedema, on_arvs, tb_hist, death_2week, death_10week,
         trial, trial_arm, country, country_iecv, country_iecv2, flucytosine, treatment) %>%
  map(~ levels(.x)) %>%
  print()

# If labelled levels are NULL, needs following command
cm_main <- cm_main %>%
  mutate_at(vars(papilloedema, on_arvs, tb_hist, death_2week, death_10week, country_iecv, country_iecv2), as.factor)

# Check which reference levels are set for factor variables - if they are listed first
# under the filter4factor then they are the reference for the regression.

#### If you need to relevel any variables for the regression

cm_main$ecog <- relevel(cm_main$ecog, ref = "Normal")
#cm_main$gcs_bin <- relevel(cm_main$gcs_bin, ref = "15")
cm_main$gcs_bin <- fct_relevel(cm_main$gcs_bin, "15", "11-14", "<=10")


# Rename and create time variable for treatment effect cph analysis

cm_main <- cm_main %>%
  rename(start_date = inclusiondate_scr,
         death_date = datedeath_dth,
         end_date = date_term_end)

cm_main <- cm_main %>%
  mutate(time_days = end_date - start_date)

cm_main <- cm_main %>%
  mutate(time_days_numeric = as.numeric(time_days))


###### Add meningism data and urea data + create 28d mortality

acta_stiffneck <- read_excel("~/Documents/ID training/ACF/Crytococcal Men project/Data/acta_meningism_Feb24.xls")
ambition_stiffneck <- read_csv("~/Documents/ID training/ACF/Crytococcal Men project/Data/Meningism_AMBITION.csv")

acta_stiffneck <- acta_stiffneck %>%
  rename(neck_stiffness = meningok_pth) %>%
  select(-visit_date_pth) %>%
  mutate(
    neck_stiffness = ifelse(neck_stiffness == 'Checked', 'Yes', 'No') %>%
      factor(levels = c('No', 'Yes'))
  )

ambition_stiffneck <- ambition_stiffneck %>%
  select(sid, meningism) %>%
  rename(study_id = sid,
         neck_stiffness = meningism) %>%
  mutate(
    neck_stiffness = ifelse(neck_stiffness == 1, 'Yes', 'No') %>%
      factor(levels = c('No', 'Yes'))
  )

acta_urea <- read_csv("~/Documents/ID training/ACF/Crytococcal Men project/Data/acta_urea_060324.csv")

cm_main <- cm_main %>%
  left_join(acta_stiffneck, by = 'study_id') %>%
  left_join(ambition_stiffneck, by = 'study_id') %>%
  mutate(
    neck_stiff = case_when(
      !is.na(neck_stiffness.x) ~ neck_stiffness.x,
      !is.na(neck_stiffness.y) ~ neck_stiffness.y,
      is.na(neck_stiffness.x) & is.na(neck_stiffness.y) ~ NA_character_,
      T ~ NA_character_
    ) %>% factor(levels = c('No', 'Yes')),
    death_28d = case_when(died == 1 & study_day_death <= 29 ~ "Died", 
                          died == 1 & study_day_death >29 | died == 0 ~ "Alive",
                          T ~ NA_character_) %>%
      factor(levels = c('Alive', 'Died')),
    creatinine = ifelse(
      creat_umol_brf < 20, NA_real_, creat_umol_brf
    )
  )


### Testing to work out Urea units issue

# ACTA
cm_urea_test <- cm_main %>%
  filter(trial == 'ACTA') %>%
  mutate(
    uc_ratio = creat_umol_brf / urea_mmol_brf
  )

acta_urea_density <- cm_urea_test %>%
  filter(country != 'Cameroon' | (country == 'Cameroon' & uc_ratio >= 5)) %>%
  ggplot(aes(x = urea_mmol_brf, fill = factor(country))) + geom_density(alpha = 0.3)
acta_urea_density

acta_ucratio_density <- cm_urea_test %>%
  ggplot(aes(x = uc_ratio, fill = factor(country))) + geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(0,50)) + geom_vline(xintercept = c(4, 8), linewidth = 0.5)
acta_ucratio_density

acta_urea_point <- ggplot(
  data = cm_urea_test[cm_urea_test$urea_mmol_brf < 20 & cm_urea_test$creat_umol_brf < 200 ,],
  aes(x = urea_mmol_brf, y = creat_umol_brf, color = factor(country))
) +
  geom_point() + geom_smooth(method = 'lm', se = F)
acta_urea_point

#Ambition
cm_ambition_ureatest <- cm_main %>%
  filter(trial == 'Ambition') %>%
  mutate(
    uc_ratio = creat_umol_brf / urea_mmol_brf
  )

ambition_urea_density <- cm_ambition_ureatest %>% 
  ggplot(aes(x = urea_mmol_brf, fill = factor(country))) + geom_density(alpha = 0.3) + coord_cartesian(xlim = c(0,75))
ambition_urea_density

ambition_ratio_density <- cm_ambition_ureatest %>% 
  ggplot(aes(x = uc_ratio, fill = factor(country))) + geom_density(alpha = 0.5) + coord_cartesian(xlim = c(0,40))
ambition_ratio_density

ambition_ratio_density_filtered <- cm_ambition_ureatest %>%
  filter(country != 'Botswana') %>%
  ggplot(aes(x = uc_ratio, fill = factor(country))) + geom_density(alpha = 0.5) + coord_cartesian(xlim = c(0,40))
ambition_ratio_density_filtered

ambition_urea_point <- cm_ambition_ureatest %>%
  filter(urea_mmol_brf < 100 & creat_umol_brf < 300) %>%
  ggplot(aes(x = urea_mmol_brf, y = creat_umol_brf, color = factor(country))) +
  geom_point() + geom_smooth(method = 'lm', se = F)
ambition_urea_point


## Convert ratios ≤ 4, set to missing ratios >4 and < 8 and leave everything > 8 in mmol (in Ambition)
# ACTA left as is, as all in mmol

cm_main <- cm_main %>%
  mutate(
    uc_ratio = creat_umol_brf / urea_mmol_brf
  ) %>%
  mutate(
    urea = case_when(
      is.na(urea_mmol_brf) ~ NA_real_,
      trial == 'ACTA' ~ urea_mmol_brf,
      trial == 'Ambition' & uc_ratio >= 8 ~ urea_mmol_brf,
      trial == 'Ambition' & uc_ratio <= 4 ~ (urea_mmol_brf / 6),
      trial == 'Ambition' & uc_ratio > 4 & uc_ratio < 8 ~ NA_real_,
      T ~ NA_real_
    )
  ) %>%
  mutate(
    urea = round(urea, 2)
  ) %>%
  mutate(
    urea = ifelse(urea > 100, NA_real_, urea)
  )


#####save the processed file as an RDS type for use in future scripts
saveRDS(cm_main, file = "cm_main_2023_10_02.rds")
