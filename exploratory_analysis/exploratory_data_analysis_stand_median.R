# Exploratory data analysis for the following descriptive
# properties of the data:
# 1. Conditional on patients who experienced and did not
#    experience AKI, what are the means and IQRs of change
#    in standardized KIM-1 at different time points?

# load package for manipulating tibbles
library(tidyverse)
# load package for reading xlsx files
library(readxl)
# load package for making plots
library(ggplot2)

# define a function that computes the median and IQR
compute_median_IQR <- function(vec) {
  return(quantile(vec, probs=c(0.25, 0.5, 0.75), na.rm=TRUE))
}

# define a function that makes plots
make_plot <- function(vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8, title, filename) {
  T1_yes <- quantile(vec1, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T2_yes <- quantile(vec2, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T3_yes <- quantile(vec3, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T4_yes <- quantile(vec4, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T1_no <- quantile(vec5, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T2_no <- quantile(vec6, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T3_no <- quantile(vec7, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  T4_no <- quantile(vec8, probs=c(0.25, 0.25, 0.5, 0.75, 0.75), na.rm=TRUE)
  
  png(filename, width=1000)
  boxplot(list(T1_yes=T1_yes, T1_no=T1_no, T2_yes=T2_yes, T2_no=T2_no,
               T3_yes=T3_yes, T3_no=T3_no, T4_yes=T4_yes, T4_no=T4_no), 
          main=title, ylab="Standardized KIM-1")
  
  return()
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
# change in NGAL and KIM-1
filtered_data <- mediation_data %>%
    select(PID, delta_KIM1_stand, delta_KIM1_stand_T2, delta_KIM1_stand_T3, delta_KIM1_stand_T4)

# merge the patients from the biomarkers data to make sure
# that the patients experiencing AKI are the same in both datasets
aki_yes_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_yes$PID)
print("cohort yes AKI")

# median and IQR for delta_KIM-1
print("median and IQR KIM-1 change T1-T0")
print(compute_median_IQR(aki_yes_data$delta_KIM1_stand))

print("median and IQR KIM-1 change T2-T0")
print(compute_median_IQR(aki_yes_data$delta_KIM1_stand_T2))

print("median and IQR KIM-1 change T3-T0")
print(compute_median_IQR(aki_yes_data$delta_KIM1_stand_T3))

print("median and IQR KIM-1 change T4-T0")
print(compute_median_IQR(aki_yes_data$delta_KIM1_stand_T4))

# now for cohort of patients that did not experience AKI
aki_no_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_no$PID)

print("")
print("cohort of no AKI")

# median and IQR for delta_KIM-1
print("median and IQR KIM-1 change T1-T0")
print(compute_median_IQR(aki_no_data$delta_KIM1_stand))

print("median and IQR KIM-1 change T2-T0")
print(compute_median_IQR(aki_no_data$delta_KIM1_stand_T2))

print("median and IQR KIM-1 change T3-T0")
print(compute_median_IQR(aki_no_data$delta_KIM1_stand_T3))

print("median and IQR KIM-1 change T4-T0")
print(compute_median_IQR(aki_no_data$delta_KIM1_stand_T4))

make_plot(aki_yes_data$delta_KIM1_stand, aki_yes_data$delta_KIM1_stand_T2,
          aki_yes_data$delta_KIM1_stand_T3, aki_yes_data$delta_KIM1_stand_T4,
          aki_no_data$delta_KIM1_stand, aki_no_data$delta_KIM1_stand_T2,
          aki_no_data$delta_KIM1_stand_T3, aki_no_data$delta_KIM1_stand_T4,
          "Q1, Median, and Q3 of Standardized KIM-1\nCompared to Baseline for AKI Yes and No Populations",
          "plots/kim1_stand_boxplot.png")

dev.off()