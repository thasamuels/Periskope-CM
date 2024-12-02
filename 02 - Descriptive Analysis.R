# Cryptococcal meningitis prognostic model
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

#####Cross tabulation of missingness by outcome

#Subset by mortality at 2 weeks, and create two seperate dfs
subset_died <- cm_main %>%
  filter(death_2week == "Died") %>%
  select(age, sex, weight, seizure, headache_dur, novision, gcscore_adm, ecog, temp, papilloedema,
         resp_rate, on_arvs, tb_hist, currtb, wcc, haem_gl, cd4, csf_OP, csf_cellcount, csf_gluc,
         csf_prot, csf_qculture, trial, trial_arm, country, death_2week)

subset_alive <- cm_main %>%
  filter(death_2week == "Alive") %>%
  select(age, sex, weight, seizure, headache_dur, novision, gcscore_adm, ecog, temp, papilloedema,
         resp_rate, on_arvs, tb_hist, currtb, wcc, haem_gl, cd4, csf_OP, csf_cellcount, csf_gluc,
         csf_prot, csf_qculture, trial, trial_arm, country, death_2week)

# in both dfs, create a missingness percentage for each selected variable before combining the two dfs
# The mutate command creates another column with Died for the died_missing data and Alive for the alive_missing
died_missing <- subset_died %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100))

alive_missing <- subset_alive %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100))

combined_missing <- bind_rows(
  mutate(died_missing, Outcome_Variable = "Died"),
  mutate(alive_missing, Outcome_Variable = "Alive")
)

#get rid of the following variable as not needed. 
combined_missing <- combined_missing %>%
  select(-death_2week)

# change the structure of the df from two long rows (died and alive, labelled as 1 and 2), to a df that has
# the key column as variable, with missingness as the second column, and that keeps Outcome_Variable as a separate third column
gathered_missing <- combined_missing %>%
  gather(key = "Variable", value = "Percentage_Missingness", -Outcome_Variable) #

# You can now reshape again, taking the Outcome_Varianble to be the column names and missingness to remain the data. Rows now variable names
reshaped_missing <- gathered_missing %>%
  spread(key = Outcome_Variable, value = Percentage_Missingness)

# rounds all values to 1dp
reshaped_missing_rounded <- reshaped_missing %>%
  mutate(across(c("Alive", "Died"), ~round(., 1)))

#tabulate
missing_table <- reshaped_missing_rounded %>% 
  mutate(
    Variable = case_when(
      Variable == "age" ~ "Age",
      Variable == "cd4" ~ "CD4 count",
      Variable == "country" ~ "Study Site",
      Variable == "csf_cellcount" ~ "CSF Cell Count",
      Variable == "csf_gluc" ~ "CSF Glucose",
      Variable == "csf_OP" ~ "CSF Opening Pressure",
      Variable == "csf_prot" ~ "CSF Protein",
      Variable == "csf_qculture" ~ "CSF Quantitative Culture",
      Variable == "currtb" ~ "Current TB infection",
      Variable == "ecog" ~ "ECOG",
      Variable == "gcscore_adm" ~ "GCS on admission",
      Variable == "haem_gl" ~ "Haemoglobin",
      Variable == "headache_dur" ~ "Duration of headache",
      Variable == "novision" ~ "Visual loss",
      Variable == "on_arvs" ~ "On ARVs",
      Variable == "papilloedema" ~ "Papilloedema",
      Variable == "resp_rate" ~ "Respiratory Rate",
      Variable == "seizure" ~ "Seizures",
      Variable == "sex" ~ "Sex",
      Variable == "tb_hist" ~ "History of TB",
      Variable == "temp" ~ "Temperature",
      Variable == "trial" ~ "Trial",
      Variable == "trial_arm" ~ "Trial arm",
      Variable == "wcc" ~ "White cell count (blood)",
      Variable == "weight" ~ "Weight",
      TRUE ~ as.character(Variable)
    )) %>%
  gt() %>%
  fmt(columns = c("Alive", "Died"),
      fns = function(x) paste0(x, "%")) %>%
  cols_label(
    Variable = "Variable",
    Alive = "Percentage Missing (Alive)",
    Died = "Percentage Missing (Died)"
  ) 
missing_table  

# Save first as html, then webshot to create an image. 
# missing_table_html <- missing_table %>%
#   gtsave(file = "missing_table_labelled.html") 
# 
# webshot2::webshot(
#   url = "missing_table_labelled.html",
#   file = "missing_table_labelled.png",
#   )

###### Looking at correlation of variables 

#Correlate continuous variables.
# predictor_columns <- c('age', 'weight', 'headache_dur', 'gcscore_adm', 'temp', 'resp_rate',
#                          'wcc', 'haem_gl')
# 
# correlation_matrix <- cor(cm_main [, predictor_columns], use = "pairwise.complete.obs")
# print(correlation_matrix)
# 
# library(reshape2) # Reshape the correlation matrix into tidy format
# correlation_tidy <- melt(correlation_matrix)
# 
# ggplot(data = correlation_tidy, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #Very little if any correlation
# 
# # Try the following for all variables:
# library(GGally)
# 
# explanatory_v <- c('ecog', 'papilloedema', 'on_arvs')
# cm_main %>%
#   #remove_labels() %>%
#   ggpairs(columns = explanatory_v)

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
#save_as_image(tbl_summary_1, "Table_1_CM.png")


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
#gtsave(mortality_table, filename = "Graphs and Tables/mortality_bycountry_table.png")


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
#save_as_image(tbl_summary_country, "Graphs and Tables/Country_characteristics.png")

###### Selecting the number of knots for RCS of continuous variables

# This is a plot to show continuous data as a density plot, with the relevant percentiles labelled. Excludes missing values.
# For loop to produce each in turn. 

continuous_vars <- c('age', 'weight', 'headache_dur', 'temp', 'resp_rate', 'wcc', 'haem_gl', 'cd4', 'csf_OP', 'csf_cellcount', 
                     'csf_gluc', 'csf_prot', 'csf_qculture')

for (var in continuous_vars) {
  plot <- ggplot(data = cm_main, aes(x = !!sym(var))) +
    geom_density(fill = "blue") +  
    geom_vline(aes(xintercept = quantile(!!rlang::sym(var), c(0.05), na.rm = TRUE)), color = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = quantile(!!rlang::sym(var), c(0.25), na.rm = TRUE)), color = "green", linetype = "dashed") +
    geom_vline(aes(xintercept = quantile(!!rlang::sym(var), c(0.5), na.rm = TRUE)), color = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = quantile(!!rlang::sym(var), c(0.75), na.rm = TRUE)), color = "green", linetype = "dashed") +
    geom_vline(aes(xintercept = quantile(!!rlang::sym(var), c(0.95), na.rm = TRUE)), color = "red", linetype = "dashed") +
    labs(title = paste("Density Plot of ", var), x = var, y = "Density")
  
  #ggsave(paste("density_plot_", var, ".png", sep = ""), plot, width = 8, height = 8)
}

# Examine the relationship between ECOG and GCS to see if there is redundancy in including both in the model

only_ecogcs <- na.omit(cm_main[c("ecog", "gcscore_adm")])
#gcs_under15 <- only_ecogcs %>% filter(gcscore_adm != 15)

ecog_gcs <- ggplot(only_ecogcs, aes(x = ecog, y = gcscore_adm)) +
  geom_violin(trim = F) +
  labs(title = "GCS vs ECOG",
       x = "ECOG",
       y = "GCS") +
  scale_y_continuous(breaks = seq(3,15, by = 2)) +
  theme_minimal()
ecog_gcs
#ggsave(file.path("Graphs and Tables", "ecog_gcs_violin.png"), plot = ecog_gcs)


# Examine the relationship between WCC, neutrophils and CD4 count. 

#Plot WCC vs Neut
wcc_neut <- ggplot(data = cm_main, aes(x = neut, y = wcc)) +
  geom_point(color = "black", shape = 1) +
  geom_smooth(method = 'lm', color = 'red', size = 0.4, se = F) +
  labs(y = "White Cell Count (x10^9/L)", x = "Neutrophil Count (x10^9/L)") +
  ggtitle("Scatterplot of White Cell Count vs. Neutrophil Count")
wcc_neut

#ggsave(file.path("Graphs and Tables", "WCC_v_neut_scatter.png"), plot = wcc_neut)

# Calculate the correlation coefficient WCC vs Neut - 0.87
correlation_wcc_neut <- cor(cm_main$wcc, cm_main$neut, use = "complete.obs")
print(paste("Correlation Coefficient:", round(correlation_wcc_neut, 2)))

# Plot WCC vs CD4

wcc_cd4 <- ggplot(data = cm_main, aes(x = cd4, y = wcc)) +
  geom_point(color = 'black', shape = 1, size = 0.5) +
  geom_smooth(method = 'lm', color = 'red', size = 0.4, se = F) +
  labs(x = 'CD4 count (x10^6/L', y = 'Total White Cell Count (x10^9/L)') +
  ggtitle('Scatterplot of Total WCC vs CD4 count')
wcc_cd4

#ggsave(file.path("Graphs and Tables", "WCC_v_CD4_scatter.png"), plot = wcc_cd4)

# Correlation WCC vs CD4 - 0.15
correlation_wcc_cd4 <- cor(cm_main$wcc, cm_main$cd4, use = "complete.obs")
print(paste("Correlation Coefficient:", round(correlation_wcc_cd4, 2)))

# 


##### Try to understand why csf_qculture is creating a matrix singularity in section 03 before log transformation
continuous_reseach_vars <- c('age', 'weight', 'neut', 'haem_gl', 'cd4', 'csf_OP', 'csf_gluc', 'csf_prot')
selected_r_vars <- c(continuous_reseach_vars, 'csf_qculture')
develop_cc_r_vars_nomiss <- na.omit(select(develop_cc, all_of(selected_r_vars)))
cor_matrix <- cor(develop_cc_r_vars_nomiss)
cor_matrix
contin_r_vars_corrplot <- corrplot(cor_matrix, method = "color", type = "lower", tl.col = "black", tl.srt = 45)
print(contin_r_vars_corrplot)

#No strong correlation with other continuous predictors. Evaluate factorial variables 

factor_r_vars <- c('sex', 'seizure', 'gcs_bin', 'ecog', 'flucytosine')
boxplot_r_list <- list()

for (var in factor_r_vars) {
  plot <- ggplot(data = develop_cc, aes(x = !!sym(var), y = csf_qculture_log)) +
    geom_boxplot() +  
    labs(x = var, y = "CSF Quant Culture (cfu/ml)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
          axis.title.y = element_text(size = 8))
  
  #ggsave(paste("boxplot_", var, "_v_csfqculture.png", sep = ""), plot, width = 8, height = 8)
  
  boxplot_r_list[[var]] <- plot
}
combined_plot <- wrap_plots(boxplot_r_list, ncol = 3) 
combined_plot <- combined_plot +
  plot_annotation(title = "Boxplots of Factor Indep Variables against CSF Q Culture",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12)))

combined_plot
#ggsave("combined_csfqculture_boxplots.png", combined_plot, width = 12, height = 8)

# Plots for CSFQC vs Death, normal and log
csfqc_v_2wmort <- ggplot(data = develop_cc, aes(x = death_2week, y = csf_qculture)) +
  geom_boxplot() +
  labs(x = "2-week Mortality", y = "CSF Quant Culture (cfu/ml)")
csfqc_v_2wmort

csfqc_v_2wmort_log <- ggplot(data = develop_cc, aes(x = death_2week, y = csf_qculture_log)) +
  geom_boxplot() +
  labs(x = "2-week Mortality", y = "Log10(CSF Quant Culture)")
csfqc_v_2wmort_log

csfqcult_list <- list(csfqc_v_2wmort, csfqc_v_2wmort_log)

csfqcult_combined <- wrap_plots(csfqcult_list, ncol = 2)
csfqcult_combined
#ggsave(file.path("Graphs and Tables", "CSFQculture_v_mortality_boxplots.png"), plot = csfqcult_combined)

develop_cc <- develop_cc %>%
  mutate(
    csf_qculture_nooutliers = ifelse(csf_qculture < 2000000, csf_qculture, NA_real_)
  )

csfqc_v_2wmort_nooutlier <- ggplot(data = develop_cc, aes(x = death_2week, y = csf_qculture_nooutliers)) +
  geom_boxplot() +
  labs(x = "2-week Mortality", y = "CSF Quant Culture (cfu/ml)")
csfqc_v_2wmort_nooutlier


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
#save_as_image(tbl_summary_devval, file.path("Graphs and Tables/Table_DevVal_CM.png"))


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


