# load the ImputationPipeline class
source("imputation_pipeline.R")

# read the mediation dataset
mediation_data <- read.csv("merged_data_06222025.csv")

# we need to read in the unscaled data for the primary and mediation datasets
# because reading in the scaled datasets will give us two datasets on different
# scales; after reading two unscaled datasets we will need to scale both of the
# datasets

# read the unscaled datasets
primary_data <- read.csv("primary_dataset_unscaled.csv")
mediation_data <- read.csv("analysis_data_07162025.csv")

head(primary_data)
head(mediation_data)