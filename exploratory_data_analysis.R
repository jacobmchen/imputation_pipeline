# Exploratory data anlysis for the following descriptive
# properties of the data:
# 1. Conditional on patients who experienced and did not
#    experience AKI, what is the average change in 
#    NGAL and KIM-1 between X-clamp off and baseline.
# 2. Plot the histogram of creatinine 48-hour post-operation
#    amongst patients who experienced AKI.

# load package for manipulating tibbles
library(tidyverse)
# load package for reading xlsx files
library(readxl)

# read the mediation data
mediation_data <- read.csv("analysis_data_11192025.csv")

# subset to patients that experienced AKI and look only at
# change in NGAL and KIM-1
filtered_data <- mediation_data %>%
    filter(aki == 1) %>%
    select(PID, delta_KIM.1, delta_NGAL)

# save change in KIM-1 and NGAL as histograms
png("plots/kim1.png")
hist(filtered_data$delta_KIM.1, main="Change in KIM-1 at X-Clamp Off Among\nAKI Patients",
     xlab="Change in KIM-1 at X-Clamp Off Compared to Baseline")
png("plots/ngal.png")
hist(filtered_data$delta_NGAL, main="Change in NGAL at X-Clamp Off Among\nAKI Patients",
     xlab="Change in NGAL at X-Clamp Off Compared to Baseline")

# print their average delta_KIM-1 and delta_NGAL
print("sample mean KIM-1 change")
print(mean(filtered_data$delta_KIM.1, na.rm=TRUE))
print("sample standard deviation KIM-1")
print(sd(filtered_data$delta_KIM.1, na.rm=TRUE))
print("sample mean NGAL change")
print(mean(filtered_data$delta_NGAL, na.rm=TRUE))
print("sample standard deviation NGAL change")
print(sd(filtered_data$delta_NGAL, na.rm=TRUE))

# read the raw creatinine data
creatinine_data <- read_excel("../../../KDIGO-AKI.xlsx")

# subset to patients that experienced AKI and look only at
# change in creatinine 48 hours post-surgery
creatinine_data <- creatinine_data %>%
    filter(KDIGO_AKI == "YES") %>%
    select(PID, `48hrs_delCr`)

# merge the patients from the biomarkers data to make sure
# that the patients experiencing AKI are the same in both datasets
filtered_data <- full_join(filtered_data, creatinine_data, by="PID")

print("sample mean change in creatinine")
print(mean(filtered_data$`48hrs_delCr`))
print("sample standard deviation change in creatinine")
print(sd(filtered_data$`48hrs_delCr`))

# plot a histogram of the change in creatinine data
png("plots/creatinine.png")
hist(filtered_data$`48hrs_delCr`, main="Change in Creatinine 48 Hours After Surgery\nAmong AKI Patients",
     xlab="Change in Creatinine 48 Hours After Surgery Compared to Baseline")
dev.off()