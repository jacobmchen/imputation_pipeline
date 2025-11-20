# file containing some data pre-processing that renames column names, etc.

# read the unscaled datasets
primary_data <- read.csv("primary_dataset_unscaled.csv")
mediation_data <- read.csv("analysis_data_11192025.csv")

primary_data$nadirDO2 <- primary_data$NadirDO2
primary_data$xclamp_duration <- primary_data$XClampDuration
primary_data$rbc_transfusion <- primary_data$transfuserbcintraoperativetotal

mediation_data$hct <- mediation_data$hgb * 3
mediation_data$PID <- NULL

primary_data$aki <- primary_data$aki_by_48h

# keep only relevant columns of the data
variable_list <- c( c("creatlst", "age", "gender", "bmi", "hypertn",
        "dm", "copd", "chf", "prior_mi", "hct", "hdef",
        "statin", "acearb", "betablocker", "rbc_transfusion"),
        c("nadirDO2", "xclamp_duration"),
        # c("delta_KIM.1", "delta_MCP.1", "delta_NGAL", "delta_YKL.40"),
        c("aki") )

primary_data <- primary_data[, variable_list]

# create a named list to keep track of the different variations
# that we need to preprocess
variations <- list()
variations[["original"]] <- c("delta_KIM.1", "delta_MCP.1", "delta_NGAL", "delta_YKL.40")
variations[["average_change"]] <- c("KIM.1_average_change", "delta_MCP.1_average_change", "delta_NGAL_average_change", "delta_YKL.40_average_change")
variations[["KIM.1"]] <- c("KIM.1")
variations[["MCP.1"]] <- c("MCP.1")
variations[["NGAL"]] <- c("NGAL")
variations[["YKL.40"]] <- c("YKL.40")

# create a vector of the names in the named list that we need to iterate through
variation_names <- c("original", "average_change", "KIM.1", "MCP.1", "NGAL", "YKL.40")

for (name in variation_names) {
    # create the variable list
    variable_list <- c( c("creatlst", "age", "gender", "bmi", "hypertn",
                          "dm", "copd", "chf", "prior_mi", "hct", "hdef",
                          "statin", "acearb", "betablocker", "rbc_transfusion"),
                        c("nadirDO2", "xclamp_duration"),
                        variations[[name]],
                        c("aki") )

    # create a subset of the observed data
    mediation_data <- mediation_data[, variable_list]

    # keep only fully observed rows
    mediation_data <- mediation_data[complete.cases(mediation_data), ]
    # save the dataset
    write.csv(mediation_data, paste0("mediation_data_final_", name, ".csv"), row.names=FALSE)
}

# keep only fully observed rows for primary_data then write to csv
primary_data <- primary_data[complete.cases(primary_data), ]
write.csv(primary_data, "primary_data_final.csv", row.names=FALSE)
