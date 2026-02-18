# Exploratory data analysis for the following descriptive
# properties of the data:
# 1. Conditional on patients who experienced and did not
#    experience AKI, what is the average change in 
#    NGAL and KIM-1 standardized by creatinine
#    between X-clamp off and baseline.
# 2. Plot the histogram of creatinine 48-hour post-operation
#    amongst patients who experienced AKI.

# load package for manipulating tibbles
library(tidyverse)

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
creatinine_data <- read_csv("serum_creat_AKI.csv")

# subset to patients that experienced AKI
creatinine_data_aki_yes <- creatinine_data %>%
  filter(AKI_delta == 1) %>%
  select(PID, AKI_delta)

creatinine_data_aki_no <- creatinine_data %>%
  filter(AKI_delta == 0) %>%
  select(PID, AKI_delta)

# read the mediation data
mediation_data <- read.csv("standardized_KIM1_NGAL.csv")

# subset to patients that experienced AKI and look only at KIM-1
filtered_data <- mediation_data %>%
    select(PID, delta_KIM1_stand, delta_KIM1_stand_T2, delta_KIM1_stand_T3, delta_KIM1_stand_T4)

# merge the patients from the biomarkers data to make sure
# that the patients experiencing AKI are the same in both datasets
aki_yes_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_yes$PID)

print("cohort yes AKI")

# save change in KIM-1 as histograms
png("plots/kim1_stand_T1.png")
hist(aki_yes_data$delta_KIM1_stand, main="Change in Standardized KIM-1 at T1 Among\nAKI Patients",
     xlab="Change in KIM-1 at T1 Compared to Baseline")

png("plots/kim1_stand_T2.png")
hist(aki_yes_data$delta_KIM1_stand_T2, main="Change in Standardized KIM-1 at T2 Among\nAKI Patients",
     xlab="Change in KIM-1 at T2 Compared to Baseline")

png("plots/kim1_stand_T3.png")
hist(aki_yes_data$delta_KIM1_stand_T3, main="Change in Standardized KIM-1 at T3 Among\nAKI Patients",
     xlab="Change in KIM-1 at T3 Compared to Baseline")

png("plots/kim1_stand_T4.png")
hist(aki_yes_data$delta_KIM1_stand_T4, main="Change in Standardized KIM-1 at T4 Among\nAKI Patients",
     xlab="Change in KIM-1 at T4 Compared to Baseline")

# print their average delta_KIM-1 and delta_NGAL
print("sample mean KIM-1 change")
print(mean(aki_yes_data$delta_KIM1_stand, na.rm=TRUE))
print("sample standard deviation KIM-1")
print(sd(aki_yes_data$delta_KIM1_stand, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_KIM1_stand))

print("sample mean KIM-1 change T2-T0")
print(mean(aki_yes_data$delta_KIM1_stand_T2, na.rm=TRUE))
print("sample standard deviation KIM-1 T2-T0")
print(sd(aki_yes_data$delta_KIM1_stand_T2, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_KIM1_stand_T2))

print("sample mean KIM-1 change T3-T0")
print(mean(aki_yes_data$delta_KIM1_stand_T3, na.rm=TRUE))
print("sample standard deviation KIM-1 T3-T0")
print(sd(aki_yes_data$delta_KIM1_stand_T3, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_KIM1_stand_T3))

print("sample mean KIM-1 change T4-T0")
print(mean(aki_yes_data$delta_KIM1_stand_T4, na.rm=TRUE))
print("sample standard deviation KIM-1 T4-T0")
print(sd(aki_yes_data$delta_KIM1_stand_T4, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_yes_data$delta_KIM1_stand_T4))

# now for cohort of patients that did not experience AKI
aki_no_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_no$PID)

print("")
print("cohort of no AKI")

png("plots/kim1_stand_aki_no_T1.png")
hist(aki_no_data$delta_KIM1_stand, main="Change in Standardized KIM-1 at T1 Among\nNo AKI Patients",
     xlab="Change in KIM-1 at T1 Compared to Baseline")

png("plots/kim1_stand_aki_no_T2.png")
hist(aki_no_data$delta_KIM1_stand_T2, main="Change in Standardized KIM-1 at T2 Among\nNo AKI Patients",
     xlab="Change in KIM-1 at T2 Compared to Baseline")

png("plots/kim1_stand_aki_no_T3.png")
hist(aki_no_data$delta_KIM1_stand_T3, main="Change in Standardized KIM-1 at T3 Among\nNo AKI Patients",
     xlab="Change in KIM-1 at T3 Compared to Baseline")

png("plots/kim1_stand_aki_no_T4.png")
hist(aki_no_data$delta_KIM1_stand_T4, main="Change in Standardized KIM-1 at T4 Among\nNo AKI Patients",
     xlab="Change in KIM-1 at T4 Compared to Baseline")

# print their average delta_KIM-1 and delta_NGAL
print("sample mean KIM-1 change")
print(mean(aki_no_data$delta_KIM1_stand, na.rm=TRUE))
print("sample standard deviation KIM-1")
print(sd(aki_no_data$delta_KIM1_stand, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_KIM1_stand))

print("sample mean KIM-1 change T2-T0")
print(mean(aki_no_data$delta_KIM1_stand_T2, na.rm=TRUE))
print("sample standard deviation KIM-1 T2-T0")
print(sd(aki_no_data$delta_KIM1_stand_T2, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_KIM1_stand_T2))

print("sample mean KIM-1 change T3-T0")
print(mean(aki_no_data$delta_KIM1_stand_T3, na.rm=TRUE))
print("sample standard deviation KIM-1 T3-T0")
print(sd(aki_no_data$delta_KIM1_stand_T3, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_KIM1_stand_T3))

print("sample mean KIM-1 change T4-T0")
print(mean(aki_no_data$delta_KIM1_stand_T4, na.rm=TRUE))
print("sample standard deviation KIM-1 T4-T0")
print(sd(aki_no_data$delta_KIM1_stand_T4, na.rm=TRUE))
print("confidence intervals")
print(compute_conf_interval_t(aki_no_data$delta_KIM1_stand_T4))

# print(t.test(aki_yes_data$delta_KIM1_stand, aki_no_data$delta_KIM1_stand, alternative="greater"))
# print(t.test(aki_yes_data$delta_KIM1_stand_T2, aki_no_data$delta_KIM1_stand_T2, alternative="greater"))
# print(t.test(aki_yes_data$delta_KIM1_stand_T3, aki_no_data$delta_KIM1_stand_T3, alternative="greater"))
# print(t.test(aki_yes_data$delta_KIM1_stand_T4, aki_no_data$delta_KIM1_stand_T4, alternative="greater"))

dev.off()