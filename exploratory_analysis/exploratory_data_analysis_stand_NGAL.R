# Exploratory data analysis for the following descriptive
# properties of the data:
# 1. Conditional on patients who experienced and did not
#    experience AKI, what is the average change in 
#    NGAL standardized by creatinine
#    between X-clamp off and baseline.

# load package for manipulating tibbles
library(tidyverse)
# load package for reading xlsx files
library(readxl)

# define a function that computes a confidence interval
compute_conf_interval_t <- function(vec) {
  n <- sum(!is.na(vec))
  mu <- mean(vec, na.rm=TRUE)
  se <- sd(vec, na.rm=TRUE)
  t_star <- qt(0.975, df=(n-1))
  
  lower <- mu - t_star * se / sqrt(n)
  upper <- mu + t_star * se / sqrt(n)
  
  return(c(lower, upper))
}

# read the raw creatinine data
creatinine_data <- read_excel("../../../../KDIGO-AKI.xlsx")

# subset to patients that experienced AKI and look only at
# change in creatinine 48 hours post-surgery
creatinine_data_aki_yes <- creatinine_data %>%
  filter(KDIGO_AKI_48hrs == "YES") %>%
  select(PID, `48hrs_delCr`)

creatinine_data_aki_no <- creatinine_data %>%
  filter(KDIGO_AKI_48hrs == "NO") %>%
  select(PID, `48hrs_delCr`)

# read the mediation data
mediation_data <- read.csv("standardized_KIM1_NGAL.csv")

# subset to patients that experienced AKI and look only at
# change in NGAL and NGAL
filtered_data <- mediation_data %>%
    select(PID, delta_NGAL_stand, delta_NGAL_stand_T2, delta_NGAL_stand_T3, delta_NGAL_stand_T4)

# merge the patients from the biomarkers data to make sure
# that the patients experiencing AKI are the same in both datasets
aki_yes_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_yes$PID)

print("cohort yes AKI")

png("plots/NGAL_stand_T1.png")
hist(aki_yes_data$delta_NGAL_stand, main="Change in Standardized NGAL at T1 Among\nAKI Patients",
     xlab="Change in NGAL at T1 Compared to Baseline")

png("plots/NGAL_stand_T2.png")
hist(aki_yes_data$delta_NGAL_stand_T2, main="Change in Standardized NGAL at T2 Among\nAKI Patients",
     xlab="Change in NGAL at T2 Compared to Baseline")

png("plots/NGAL_stand_T3.png")
hist(aki_yes_data$delta_NGAL_stand_T3, main="Change in Standardized NGAL at T3 Among\nAKI Patients",
     xlab="Change in NGAL at T3 Compared to Baseline")

png("plots/NGAL_stand_T4.png")
hist(aki_yes_data$delta_NGAL_stand_T4, main="Change in Standardized NGAL at T4 Among\nAKI Patients",
     xlab="Change in NGAL at T4 Compared to Baseline")

# print their average delta_NGAL and delta_NGAL
print("sample mean NGAL change")
print(mean(aki_yes_data$delta_NGAL_stand, na.rm=TRUE))
print("sample standard deviation NGAL")
print(sd(aki_yes_data$delta_NGAL_stand, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_NGAL_stand))

print("sample mean NGAL change T2-T0")
print(mean(aki_yes_data$delta_NGAL_stand_T2, na.rm=TRUE))
print("sample standard deviation NGAL T2-T0")
print(sd(aki_yes_data$delta_NGAL_stand_T2, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_NGAL_stand_T2))

print("sample mean NGAL change T3-T0")
print(mean(aki_yes_data$delta_NGAL_stand_T3, na.rm=TRUE))
print("sample standard deviation NGAL T3-T0")
print(sd(aki_yes_data$delta_NGAL_stand_T3, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_NGAL_stand_T3))

print("sample mean NGAL change T4-T0")
print(mean(aki_yes_data$delta_NGAL_stand_T4, na.rm=TRUE))
print("sample standard deviation NGAL T4-T0")
print(sd(aki_yes_data$delta_NGAL_stand_T4, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_NGAL_stand_T4))

# now for cohort of patients that did not experience AKI
aki_no_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_no$PID)

print("")
print("cohort of no AKI")

# print(t.test(aki_yes_data$delta_NGAL_stand, aki_no_data$delta_NGAL_stand, alternative="greater"))
# print(t.test(aki_yes_data$delta_NGAL_stand_T2, aki_no_data$delta_NGAL_stand_T2, alternative="greater"))
# print(t.test(aki_yes_data$delta_NGAL_stand_T3, aki_no_data$delta_NGAL_stand_T3, alternative="greater"))
# print(t.test(aki_yes_data$delta_NGAL_stand_T4, aki_no_data$delta_NGAL_stand_T4, alternative="greater"))

png("plots/NGAL_stand_aki_no_T1.png")
hist(aki_no_data$delta_NGAL_stand, main="Change in Standardized NGAL at X-Clamp Off Among\nNo AKI Patients",
     xlab="Change in NGAL at X-Clamp Off Compared to Baseline")

png("plots/NGAL_stand_aki_no_T2.png")
hist(aki_no_data$delta_NGAL_stand_T2, main="Change in Standardized NGAL at T2 Among\nNo AKI Patients",
     xlab="Change in NGAL at T2 Compared to Baseline")

png("plots/NGAL_stand_aki_no_T3.png")
hist(aki_no_data$delta_NGAL_stand_T3, main="Change in Standardized NGAL at T3 Among\nNo AKI Patients",
     xlab="Change in NGAL at T3 Compared to Baseline")

png("plots/NGAL_stand_aki_no_T4.png")
hist(aki_no_data$delta_NGAL_stand_T4, main="Change in Standardized NGAL at T4 Among\nNo AKI Patients",
     xlab="Change in NGAL at T4 Compared to Baseline")

# print their average delta_NGAL and delta_NGAL
print("sample mean NGAL change")
print(mean(aki_no_data$delta_NGAL_stand, na.rm=TRUE))
print("sample standard deviation NGAL")
print(sd(aki_no_data$delta_NGAL_stand, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_NGAL_stand))

print("sample mean NGAL change T2-T0")
print(mean(aki_no_data$delta_NGAL_stand_T2, na.rm=TRUE))
print("sample standard deviation NGAL T2-T0")
print(sd(aki_no_data$delta_NGAL_stand_T2, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_NGAL_stand_T2))

print("sample mean NGAL change T3-T0")
print(mean(aki_no_data$delta_NGAL_stand_T3, na.rm=TRUE))
print("sample standard deviation NGAL T3-T0")
print(sd(aki_no_data$delta_NGAL_stand_T3, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_NGAL_stand_T3))

print("sample mean NGAL change T4-T0")
print(mean(aki_no_data$delta_NGAL_stand_T4, na.rm=TRUE))
print("sample standard deviation NGAL T4-T0")
print(sd(aki_no_data$delta_NGAL_stand_T4, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_NGAL_stand_T4))

dev.off()