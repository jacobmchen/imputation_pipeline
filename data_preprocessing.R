# file containing some data pre-processing that renames column names, etc.

# read the unscaled datasets
primary_data <- read.csv("primary_dataset_unscaled.csv")
mediation_data <- read.csv("analysis_data_07162025.csv")

head(primary_data)
head(mediation_data)

mediation_data$hct <- mediation_data$hgb * 3

primary_data$aki <- primary_data$aki_by_48h