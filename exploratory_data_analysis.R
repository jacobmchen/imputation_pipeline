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
mediation_data <- read.csv("mediation_data_final_original.csv")

# subset to patients that experienced AKI and look only at
# change in NGAL and KIM-1
filtered_data <- mediation_data %>%
    filter(aki == 1) %>%
    select(delta_KIM_1, delta_NGAL)

# print their average delta_KIM-1 and delta_NGAL
print("sample mean KIM-1 change")
print(mean(filtered_data$delta_KIM_1))
print("sample standard deviation KIM-1")
print(sd(filtered_data$delta_KIM_1))
print("sample mean NGAL change")
print(mean(filtered_data$delta_NGAL))
print("sample standard deviation NGAL change")
print(sd(filtered_data$delta_NGAL))

# read the raw creatinine data
creatinine_data <- read_excel("../../../KDIGO-AKI.xlsx")

# subset to patients that experienced AKI and look only at
# change in creatinine 48 hours post-surgery
filtered_data <- creatinine_data %>%
    filter(KDIGO_AKI == "YES") %>%
    select(48hrs_delCr)

# plot a histogram of the change in creatinine data
hist(filtered_data$48hrs_delCr)
