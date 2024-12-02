# Cryptococcal meningitis prognostic model - Periskope-CM
# 02: Descriptive Analysis
# Description: Analysis of Missingness, Description of cohort
# Started: 04/10/23
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
library(flextable)
library(webshot2)
library(corrplot)
library(patchwork)
library(DiagrammeR)

cm_main <- read_rds("cm_main_2023_10_02.rds")

develop_cc <- cm_main %>% filter(trial == "ACTA" | (trial == "Ambition" & country != "Malawi"))
valid_cc <- cm_main %>% filter(country == "Malawi" & trial == "Ambition")

#Proxy variable to allow investigation of missingness etc by trial
data_acta <- subset(cm_main, trial == "ACTA")
data_ambition <- subset(cm_main, trial == "Ambition")

samplesize_total <- as.numeric(length(cm_main$study_id))
samplesize_acta <- as.numeric(length(data_acta$study_id))
samplesize_ambition <- as.numeric(length(data_ambition$study_id))
samplesize_dev <- as.numeric(length(develop_cc$study_id))
samplesize_val <- as.numeric(length(valid_cc$study_id))

####### Output missingness for variables of interest
missing_covariates <- cm_main %>%
  select(age, sex, weight, seizure, headache_dur, novision, gcscore_adm, ecog, temp, papilloedema,
         resp_rate, on_arvs, tb_hist, currtb, wcc, haem_gl, cd4, csf_OP, csf_cellcount, csf_gluc,
         csf_prot, csf_qculture, trial, trial_arm, country)

missing_covariates_summ <- missing_covariates %>%
  summarise_all(~ sum(is.na(.))/length(.) * 100)

###View in table
mcs_df <- data.frame(
  Variable = names(missing_covariates_summ),
  Percentage_missingness = as.numeric(missing_covariates_summ)
)
mcs_df$Percentage_missingness <- round(mcs_df$Percentage_missingness, 1)

table_missing <- mcs_df %>%
  gt() %>%
  tab_header(title = "Data missingness") %>%
  cols_label(Variable = "Variable", Percentage_missingness = "Percentage Missingness")
table_missing

###Graph
ggplot(missing_covariates_summ %>%
         gather(variable, missingness, everything()),
       aes(x = variable, y = missingness, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(x = "Variables", y = "Percentage Missingness") +
  theme(axis.text = element_text(angle = 90, hjust = 1))

#current tb is the only variable that does not meet the 60% presence requirement to be taken forward to MI and analysis


####### Label variable names

label(cm_main$death_2week) <- "Mortality at 2 weeks"
label(cm_main$death_10week) <- "Mortality at 10 weeks"
label(cm_main$trial) <- "Trial"
label(cm_main$trial_arm) <- "Treatment arm"
label(cm_main$treatment) <- "Treatment received"
label(cm_main$country) <- "Study site"
label(cm_main$age) <- "Age; years"
label(cm_main$sex) <- "Sex"
label(cm_main$weight) <- "Weight; kg"
label(cm_main$seizure) <- "Seizures"
label(cm_main$headache_dur) <- "Duration of headache; days"
label(cm_main$novision) <- "Visual loss"
label(cm_main$gcs_bin) <- "GCS score"
label(cm_main$ecog) <- "ECOG performance status"
label(cm_main$temp) <- "Temperature; °C"
label(cm_main$papilloedema) <- "Papilloedema"
label(cm_main$resp_rate) <- "Respiratory rate; per min"
label(cm_main$on_arvs) <- "Taking ARVs on admission"
label(cm_main$tb_hist) <- "History of prior TB"
label(cm_main$currtb) <- "Active TB co-infection"
label(cm_main$wcc) <- "White cell count; x10^9/L"
label(cm_main$neut) <- "Neutrophil count; x10^9/L"
label(cm_main$haem_gl) <- "Haemoglobin; g/L"
label(cm_main$cd4) <- "CD4 count; x10^6/L"
label(cm_main$csf_OP) <- "CSF opening pressure; cmH20"
label(cm_main$csf_op_d7) <- "CSF opening pressure day 7; cmH20"
label(cm_main$csf_op_d14) <- "CSF opening pressure day 14; cmH20"
label(cm_main$csf_cellcount) <- "CSF cell count; WBC per mm3"
label(cm_main$csf_gluc) <- "CSF Glucose; mmol/L"
label(cm_main$csf_prot) <- "CSF Protein; mg/mL"
label(cm_main$csf_qculture) <- "CSF Quantitative culture; CFU/mL"
label(cm_main$csf_qculture_d7) <- "CSF quantitative culture day 7; CFU/mL"
label(cm_main$csf_qculture_d14) <- "CSF quantitative culture day 14; CFU/mL"
label(cm_main$haem_d7) <- "Haemoglobin day 7; g/L"
label(cm_main$haem_d14) <- "Haemoglobin day 14; g/L"
label(cm_main$creat) <- "Creatinine; umol/L"
label(cm_main$creat_d7) <- "Creatinine day 7; umol/L"
label(cm_main$creat_d14) <- "Creatinine day 14; umol/L"
label(cm_main$csf_qculture_log) <- "log(CSF quantitative culture)"
label(cm_main$csf_qculture_d7_log) <- "CSF log10(quantitative culture) day 7"
label(cm_main$csf_qculture_d14_log) <- "CSF log10(quantitative culture) day 14"


########### Create Table 1

tbl_summary_1 <- cm_main %>%
  select(death_2week, country, trial, trial_arm, age, sex, weight, seizure, gcs_bin, ecog, 
         wcc, neut, haem_gl, cd4, csf_OP, csf_cellcount, csf_qculture_log
  ) %>%
  tbl_summary(type = list(),
              by = death_2week,
              digits = list(weight ~ c(0,0,0),
                            country ~ c(0,0),
                            trial ~ c(0,0),
                            trial_arm ~ c(0,0),
                            ecog ~ c(0,0),
                            gcs_bin ~ c(0,0),
                            wcc ~ c(2,2,2),
                            haem_gl ~ c(0,0,0)),
              missing = "ifany",
              missing_text = "Missing",
              statistic = list(all_continuous() ~ "{median} ({p25} to {p75})")) %>%
  add_overall() %>%
  as_flex_table()
tbl_summary_1


# Table of mortality by country
mortality_data <- cm_main %>%
  group_by(country) %>%
  mutate(Mortality_Percentage = round(sum(death_2week == "Died") / n() * 100, 1),
         Number_Participants = n()) %>%
  select(country, Mortality_Percentage, Number_Participants) %>%
  distinct(country, .keep_all = T)

mortality_data$Country <- mortality_data$country
mortality_data$country <- NULL


mortality_table <- mortality_data %>% gt() %>%
  cols_label(
    Country = "Country",
    Mortality_Percentage = "Mortality (%)",
    Number_Participants = "Participants"
  )%>%
  cols_move_to_start("Country")
mortality_table


### Table of key predictors by country

tbl_summary_country <- cm_main %>%
  select(death_2week, country, age, sex, weight, seizure, gcs_bin, ecog,
         neut, haem_gl, cd4, csf_OP, csf_cellcount, csf_qculture, 
         trial_arm
  ) %>%
  tbl_summary(type = list(),
              by = country,
              digits = list(weight ~ c(0,0,0),
                            neut ~ c(2,2,2),
                            haem_gl ~ c(0,0,0)),
              missing = "no",
              missing_text = "Missing") %>%
  add_overall() %>%
  as_flex_table()
tbl_summary_country



#### Tabulate key predictors in development vs validation dataset, including missingness

dev_val_dataset <- cm_main %>%
  mutate(
    dataset = case_when(
      trial == "ACTA" | (trial == "Ambition" & country != "Malawi") ~ "Development",
      country == "Malawi" & trial == "Ambition" ~ "Validation"
    )
  )

tbl_summary_devval <- dev_val_dataset %>%
  select(death_2week, age, sex, weight, seizure, trial_arm,
         gcs_bin, ecog, wcc,
         neut, haem_gl, cd4, 
         csf_OP, csf_cellcount, csf_qculture_log,
         dataset,
  ) %>%
  tbl_summary(type = list(),
              by = dataset,
              digits = list(weight ~ c(0,0,0),
                            neut ~ c(2,2,2),
                            haem_gl ~ c(0,0,0),
                            trial_arm ~ c(0,0),
                            ecog ~ c(0,0),
                            gcs_bin ~ c(0,0),
                            death_2week ~ c(0,0),
                            cd4 ~ c(0,0,0),
                            csf_OP ~ c(0,0,0),
                            csf_qculture_log ~ c(2,2,2)),
              missing = "ifany",
              missing_text = "Missing",
              statistic = list(all_continuous() ~ "{median} ({p25} to {p75})")) %>%
  add_overall() %>%
  as_flex_table()
tbl_summary_devval



# Save sample sizes for flow chart
acta_death_2week <- cm_main %>% filter(trial == 'ACTA' & death_2week == 'Died') %>% nrow()
ambition_death_2week <- cm_main %>% filter(trial == 'Ambition' & death_2week == 'Died') %>% nrow()
acta_death_10week <- cm_main %>% filter(trial == 'ACTA' & death_10week == 'Died') %>% nrow()
ambition_death_10week <- cm_main %>% filter(trial == 'Ambition' & death_10week == 'Died') %>% nrow()

total_death_2week <- cm_main %>% filter(death_2week == 'Died') %>% nrow()
total_death_10week <- cm_main %>% filter(death_10week == 'Died') %>% nrow()

dev_death_2week <- cm_main %>% filter((trial == "ACTA" | (trial == "Ambition" & country != "Malawi")) & 
                                        death_2week == 'Died') %>% nrow()
dev_death_10week <- cm_main %>% filter((trial == "ACTA" | (trial == "Ambition" & country != "Malawi")) & 
                                         death_10week == 'Died') %>% nrow()
val_death_2week <- cm_main %>% filter((country == "Malawi" & trial == "Ambition") & death_2week == 'Died') %>% nrow()
val_death_10week <- cm_main %>% filter((country == "Malawi" & trial == "Ambition") & death_10week == 'Died') %>% nrow()


flowchart_data <- list(ss_acta = samplesize_acta,
                       total_acta = (samplesize_acta + 4),
                       acta_miss = 4,
                       amb_miss = 0,
                       d2_acta = acta_death_2week,
                       d10_acta = acta_death_10week,
                       ss_amb = samplesize_ambition,
                       d2_amb = ambition_death_2week,
                       d10_amb = ambition_death_10week,
                       ss_total = samplesize_total,
                       d2_total =total_death_2week,
                       d10_total = total_death_10week,
                       ss_dev = samplesize_dev,
                       d2_dev = dev_death_2week,
                       d10_dev = dev_death_10week,
                       ss_val = samplesize_val,
                       d2_val = val_death_2week,
                       d10_val = val_death_10week)

# Practice flow diagram code

DiagrammeR::grViz("Flowdiagram.dot")


