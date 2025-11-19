# this file will preprocess the biomarker data to create a new column
# representing the average of each biomarker's change compared to baseline
# at the 4 time measured time points of each biomarker
# in math, (time_1 - time_0 + time_2 - time_0 + time_3 - time_0 + time_4 - time_0)/4
# = (time_1 + time_2 + time_3 + time_4 - 4*time_0)/4
# = (time_1+time_2+time_3+time_4)/4 - time_0

# read the mediation data
mediation_data <- read.csv("analysis_data_07162025.csv")

# save a vector of the names of the three biomarkers
biomarker_names <- c("KIM.1", "MCP.1", "NGAL", "YKL.40")

# iterate through each of the biomarker names
for (biomarker in biomarker_names) {
    # declare an empty vector to save the average metric that we are computing
    average_change <- c()

    # get the whole columns for the baseline and 4 time points
    baseline <- mediation_data[[paste0(biomarker, "_0")]]
    time_1 <- mediation_data[[paste0(biomarker, "_1")]]
    time_2 <- mediation_data[[paste0(biomarker, "_2")]]
    time_3 <- mediation_data[[paste0(biomarker, "_3")]]
    time_4 <- mediation_data[[paste0(biomarker, "_4")]]
    average_change <- (time_1 + time_2 + time_3 + time_4)/4 - baseline

    mediation_data[[paste0(biomarker, "_average_change")]] <- average_change
}

# save the new dataset
write.csv(mediation_data, "analysis_data_11192025.csv", row.names=FALSE)
