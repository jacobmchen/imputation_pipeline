# use object oriented S4 programming in R to set up an imputation pipeline that
# we can use more easily with the data application

# library for learning continuous densities
library(haldensify)

# define an expit function
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

# the syntax of declaring a new class requires two parameters, the first
# is the name of the class, and the second is a list of attributes with
# the data types
# for attributes that will be populated later, we leave the data type as
# ANY so that it can be a NULL value
setClass("ImputationPipeline", 
         slots=list(primary_data="data.frame", # the primary, larger dataset without biomarker data
                    mediation_data="data.frame", # the smaller dataset with biomarker data
                    variable_dictionary="list", # a dictionary mapping canonical variable names to the
                                                # variable names in the dataframes
                    mediation_models="ANY", # a list of mediator densities learned from the haldensify package
                    treatment_models="ANY" # a list of treatment densities learned from the haldensify package
                    ))

# create a constructor
NewImputationPipeline <- function(primary_data, mediation_data, variable_dictionary) {
  new("ImputationPipeline", primary_data=primary_data, 
                mediation_data=mediation_data,
                variable_dictionary=variable_dictionary,
                mediation_models=NULL,
                treatment_models=NULL)
}

# define the show method for the ImputationPipeline class to define how it should
# look when printed, we don't need to set a generic for the show method since its
# generic is built-in for S4 classes
setMethod("show", "ImputationPipeline", function(object) {
  print("primary_data:")
  print(head(object@primary_data))

  print("secondary_data:")
  print(head(object@mediation_data))

  print("variable_dictionary:")
  print(object@variable_dictionary)

  print("length of mediation_models:")
  print(length(object@mediation_models))
})

# define a helper method for learning the density of a variable given other variables;
# this method is intended to be called only in methods of the class ImputationPipeline
# read_from_rds: boolean variable stating whether we learn the model from data
# or read it from a previously saved RDS file
# exposure: string containing the name of the exposure to model
# adjustment_set: vector of strings containing the proper adjustment set
# analysis_data: the dataset
learn_density <- function(read_from_rds=FALSE, exposure, adjustment_set, analysis_data, filename_addon="") {
  # compute the number of observations in our dataset
  n_obs <- nrow(analysis_data)
  
  # learn the density for f(exposure | adjustment_set)
  # make a copy of the data
  analysis_data_copy <- data.table::copy(analysis_data)
  
  # get the exposure of interest
  A <- analysis_data_copy[[exposure]]
  
  # get the covariates
  L <- analysis_data_copy[, adjustment_set]
  
  # define the name of save file
  filename_model <- paste0("model_", exposure, "_", filename_addon, ".RDS")
    
  if (read_from_rds) {
    model <- readRDS(filename_model)
  
    return(model)
  } else {
    # use parameters specified in the vignette, except modify length of lambda sequence and max_degree
    model <- haldensify(A=A,  W=L,
                          n_bins=c(20,30),
                          #n_bins=round(c(1.5, 2) * sqrt(n_obs)),
                          grid_type="equal_range",
                          lambda_seq=exp(seq(-0.1, -10, length=20)),
                          max_degree=2,
                          reduce_basis=1/sqrt(n_obs))
    
    # save the learned model as an RDS file
    saveRDS(model, filename_model)

    return(model)
  }
}

# set the generic so that the method learnMediationDensities exists in the method table for
# the S4 dispatcher
setGeneric("learnMediationDensities", function(object, read_from_rds, filename)
                            standardGeneric("learnMediationDensities"))

# define a method for learning the mediation models
setMethod("learnMediationDensities", "ImputationPipeline", function(object, read_from_rds, filename) {
  # get the full vector of mediation variables
  mediation_variables <- object@variable_dictionary[["M"]]

  # get the full vector of covariate variables
  covariate_variables <- object@variable_dictionary[["X"]]

  # get the number of mediation variables
  med_length <- length(mediation_variables)

  # declare an empty list to save the learned models in
  models <- list()

  # learn a density for each of the mediation variables
  for (index in 1:med_length) {
    # get current M we want to learn the density for
    cur_M <- mediation_variables[index]

    # chain rule the rest of the densities
    adjustment_M <- mediation_variables[(index+1):med_length]

    # if this is the final density to learn, there is nothing to chain rule
    if (index == med_length) adjustment_M <- c()

    # define the covariate set
    adjustment_set <- c(covariate_variables, adjustment_M)

    # learn the density from the mediation data
    cur_model <- learn_density(read_from_rds, cur_M, adjustment_set, object@mediation_data, filename)

    # save the learned model to the list
    models[[index]] <- cur_model
  }

  # save the mediation models to the slot
  object@mediation_models <- models

  # return the modified pipeline object
  return(object)
})

# Simulation study starts here
# define the size of the primary dataset
primary_n <- 1000
# define the size of the secondary dataset
secondary_n <- 150

# define the standard deviation
sd = 1

# generate variables for the primary dataset
X_primary <- rnorm(primary_n, 0, 1)
A_primary <- X_primary + rnorm(primary_n, 0, 0.1*sd)
M_primary <- A_primary + X_primary + rnorm(primary_n, 0, 0.1*sd)
Y_primary <- rbinom(primary_n, 1, expit(A_primary + X_primary + M_primary))

# save the primary variables into a dataframe
primary_data <- data.frame(
  X=X_primary,
  A=A_primary,
  Y=Y_primary
)

# generate variables for the secondary dataset; it is the same DGP as
# above
X_secondary <- rnorm(secondary_n, 0, 1)
A_secondary <- X_secondary + rnorm(secondary_n, 0, 0.1*sd)
M_secondary <- A_secondary + X_secondary + rnorm(secondary_n, 0, 0.1*sd)
Y_secondary <- rbinom(secondary_n, 1, expit(A_secondary + X_secondary + M_secondary))

# save the variables with the mediation data into the secondary
# dataset
mediation_data <- data.frame(
  X=X_secondary,
  A=A_secondary,
  M=M_secondary,
  Y=Y_secondary
)

# create a named list, which will behave like a Python dictionary
# this named list maps canonical variable names to a vector of what the variable
# names are in the dataframe
variable_dictionary <- list(
  "X" = c("X"),
  "A" = c("A"),
  "M" = c("M"),
  "Y" = c("Y")
)

# call the predefined constructor for creating a new imputation pipeline object
pipeline <- NewImputationPipeline(primary_data, mediation_data, variable_dictionary)

# learn the mediation densities and update the object
pipeline <- learnMediationDensities(pipeline, TRUE, "test")

# impute M values for the primary dataset


show(pipeline)