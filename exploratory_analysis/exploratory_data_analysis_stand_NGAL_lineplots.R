# Exploratory data analysis for the following descriptive
# properties of the data:
# 1. Conditional on patients who experienced and did not
#    experience AKI, what is the progression of standardized
#    NGAL from T1 to T3?

# load package for manipulating tibbles
library(tidyverse)
# load package for reading xlsx files
library(readxl)
# load package for making plots
library(ggplot2)

# read the raw creatinine data
creatinine_data <- read_csv("serum_creat_AKI.csv")

# subset to patients that experienced AKI and look only at
# change in creatinine 48 hours post-surgery
creatinine_data_aki_yes <- creatinine_data %>%
  filter(AKI_delta == 1) %>%
  select(PID, AKI_delta)

creatinine_data_aki_no <- creatinine_data %>%
  filter(AKI_delta == 0) %>%
  select(PID, AKI_delta)

# read the mediation data
mediation_data <- read.csv("standardized_KIM1_NGAL.csv")

# subset to patients that experienced AKI and look only at
# change in NGAL and KIM-1
# filtered_data <- mediation_data %>%
#     select(PID, delta_NGAL_stand, delta_NGAL_stand_T2, delta_NGAL_stand_T3, delta_NGAL_stand_T4) %>%
#     rename(T1=delta_NGAL_stand, T2=delta_NGAL_stand_T2, T3=delta_NGAL_stand_T3, T4=delta_NGAL_stand_T4)
filtered_data <- mediation_data %>%
  select(PID, delta_NGAL_stand_T4) %>%
  rename(T4=delta_NGAL_stand_T4)

# merge the patients from the biomarkers data to make sure
# that the patients experiencing AKI are the same in both datasets
aki_yes_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_yes$PID)
print("cohort yes AKI")

# transform the data into long form with only three columns: PID, time, and 
# standardized NGAL value
aki_yes_data_long <- pivot_longer(
  aki_yes_data,
  cols = starts_with("T"),
  names_to = "time",
  values_to = "NGAL_stand"
)

# factor the time labels so that R knows what order they are supposed to be in
aki_yes_data_long$time <- factor(aki_yes_data_long$time, levels=c("T4"))

# make the line plot by grouping with PID and coloring points based on PID
png("plots/NGAL_stand_lineplot_akiyes_T4.png")
ggplot(aki_yes_data_long, aes(x=time, y=NGAL_stand, group=PID, color=PID)) +
  geom_line(size=1) +
  geom_point(size=3) +
  theme_minimal() +
  labs(x="Time", y="Standardized NGAL", title="Standardized NGAL Progression Over Time for Yes AKI")

# now for cohort of patients that did not experience AKI
aki_no_data <- filtered_data %>%
  filter(PID %in% creatinine_data_aki_no$PID)

print("")
print("cohort of no AKI")

aki_no_data_long <- pivot_longer(
  aki_no_data,
  cols = starts_with("T"),
  names_to = "time",
  values_to = "NGAL_stand"
)

aki_no_data_long$time <- factor(aki_no_data_long$time, levels=c("T4"))

png("plots/NGAL_stand_lineplot_akino_T4.png")
ggplot(aki_no_data_long, aes(x=time, y=NGAL_stand, group=PID, color=PID)) +
  geom_line(size=1) +
  geom_point(size=3) +
  theme_minimal() +
  labs(x="Time", y="Standardized NGAL", title="Standardized NGAL Progression Over Time for No AKI")

dev.off()