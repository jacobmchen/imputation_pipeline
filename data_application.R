# load the ImputationPipeline class
source("imputation_pipeline.R")

# set the seed
set.seed(0)

# we need to read in the unscaled data for the primary and mediation datasets
# because reading in the scaled datasets will give us two datasets on different
# scales; after reading two unscaled datasets we will need to scale both of the
# datasets

# read the unscaled datasets
primary_data <- read.csv("primary_data_final.csv")
print(paste("primary data size:", nrow(primary_data)))
mediation_data <- read.csv("mediation_data_final.csv")
print(paste("mediation data size:", nrow(mediation_data)))

# define the variable dictionary
variable_dictionary <- list(
    "X" = c("creatlst", "age", "gender", "bmi", "hypertn",
            "dm", "copd", "chf", "prior_mi", "hct", "hdef",
            "statin", "acearb", "betablocker", "rbc_transfusion"),
    "A" = c("nadirDO2", "xclamp_duration"),
    "M" = c("delta_KIM.1", "delta_MCP.1", "delta_NGAL", "delta_YKL.40"),
    "Y" = c("aki")
)

# save the variables that should be excluded from standardization
to_exclude <- c("gender", "hypertn", "dm", "copd", "chf", "prior_mi", "statin",
                "acearb", "betablocker")

# define whether we should read model data from RDS data
read_from_rds <- TRUE

# create the pipeline object
pipeline <- NewImputationPipeline(primary_data, mediation_data, variable_dictionary)

# standardize the data
pipeline <- standardizeData(pipeline, to_exclude)

# learn the mediation densities and update the object
pipeline <- learnMediationDensities(pipeline, read_from_rds, "07172025")

# impute M values for the primary dataset, set the seed to be 0
pipeline <- imputeMediators(pipeline, 0)

# learn the treatment densities and update the object
pipeline <- learnTreatmentDensities(pipeline, read_from_rds, "07172025")

# learn the marginal treatment densities and update the object
pipeline <- learnMarginalTreatmentDensities(pipeline, read_from_rds, "07172025")

# define the interventional values that we're interested in
a_prime_vals <- c(0.5, 5)
a_vals <- c(4, 1.5)

# estimate the mediation term
pipeline <- computePseudoOutcome(pipeline, a_prime_vals, a_vals)
# print("***debugging")
# print(mean(pipeline@Y_p))
# print(pipeline@Y_p)
# print(mean(pipeline@combined_imputed_data[["aki"]]))

# compute the MSM weights
pipeline <- computeMSMWeights(pipeline)

# estimate the counterfactual terms
print("bart estimates")
pipeline <- estimateEffects(pipeline, a_prime_vals, a_vals, adjust=TRUE)
print(pipeline@total_effect)
print(pipeline@indirect_effect)
print(pipeline@direct_effect)

print("bart estimates, no adjustment")
pipeline <- estimateEffects(pipeline, a_prime_vals, a_vals, adjust=FALSE)
print(pipeline@total_effect)
print(pipeline@indirect_effect)
print(pipeline@direct_effect)

print("linear estimators")
pipeline <- estimateEffectsLinear(pipeline, a_prime_vals, a_vals, adjust=TRUE)
print(pipeline@total_effect)
print(pipeline@indirect_effect)
print(pipeline@direct_effect)

print("linear estimators, no adjustment")
pipeline <- estimateEffectsLinear(pipeline, a_prime_vals, a_vals, adjust=FALSE)
print(pipeline@total_effect)
print(pipeline@indirect_effect)
print(pipeline@direct_effect)
