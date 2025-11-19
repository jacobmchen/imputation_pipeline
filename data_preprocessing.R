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

variable_list <- c( c("creatlst", "age", "gender", "bmi", "hypertn",
                      "dm", "copd", "chf", "prior_mi", "hct", "hdef",
                      "statin", "acearb", "betablocker", "rbc_transfusion"),
                    c("nadirDO2", "xclamp_duration"),
                    c("delta_KIM.1", "delta_MCP.1", "delta_NGAL", "delta_YKL.40"),
                    c("KIM.1_average_change", "MCP.1_average_change", "NGAL_average_change", "YKL.40_average_change"),
                    c("aki") )

mediation_data <- mediation_data[, variable_list]

# keep only fully observed rows
primary_data <- primary_data[complete.cases(primary_data), ]
mediation_data <- mediation_data[complete.cases(mediation_data), ]

# write.csv(primary_data, "primary_data_final.csv", row.names=FALSE)
write.csv(mediation_data, "mediation_data_final.csv", row.names=FALSE)