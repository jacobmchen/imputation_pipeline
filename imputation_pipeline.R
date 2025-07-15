# use object oriented S4 programming in R to set up an imputation pipeline that
# we can use more easily with the data application

# library for learning continuous densities
library(haldensify)

# library for imputing missing values
library(mice)

# library for random forests
library(ranger)

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
                    combined_imputed_data="data.frame", # a placeholder for the combined imputed data
                    mediation_models="ANY", # a list of mediator densities/estimations learned from the haldensify package
                    treatment_models="ANY", # a list of treatment densities/estimations learned from the haldensify package
                    treatment_marginal_models="ANY", # a list of marginal treatment densities/estimations
                    MSM_weights="ANY", # a vector storing the calculation of MSM weights
                    Y_p="ANY", # a vector storing the pseudo-outcome
                    mediation_term="numeric", # a number to store the estimate of the mediation term
                    counterfactual_a="numeric", # a number to store the estimate of intervening on a
                    counterfactual_a_prime="numeric" # a number to store the estimate of intervening on a'
                    ))

# create a constructor
NewImputationPipeline <- function(primary_data, mediation_data, variable_dictionary) {
  new("ImputationPipeline", primary_data=primary_data, 
                mediation_data=mediation_data,
                variable_dictionary=variable_dictionary,
                combined_imputed_data=data.frame(),
                mediation_models=NULL,
                treatment_models=NULL,
                treatment_marginal_models=NULL,
                MSM_weights=NULL,
                Y_p=NULL,
                mediation_term=0,
                counterfactual_a=0,
                counterfactual_a_prime=0)
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

  print("length of treatment_models:")
  print(length(object@treatment_models))

  print("length of treatment_marginal_models:")
  print(length(object@treatment_marginal_models))

  print("combined_imputed_data:")
  print(head(object@combined_imputed_data))

  print("weights:")
  print(head(object@MSM_weights))

  print("pseudo-outcome:")
  print(head(object@Y_p))
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
  filename_density <- paste0("density_", exposure, "_", filename_addon, ".RDS")
    
  if (read_from_rds) {
    model <- readRDS(filename_model)

    density <- readRDS(filename_density)
  
    return(list(model, density))
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

    density <- predict(model, new_A=A, new_W=L)
    saveRDS(density, filename_density)

    return(list(model, density))
  }
}

# helper function for creating a formula
create_formula <- function(outcome, features, two_way=FALSE) {
  if (length(features) == 0) {
    return(paste0(outcome, "~1"))
  } else {
    if (two_way) {
      return(paste0(outcome, "~(", paste(features, collapse="+"), ")^2"))
    }
    return(paste0(outcome, "~", paste(features, collapse="+")))
  }
}

# set the generic so that the method learnMediationDensities exists in the method table for
# the S4 dispatcher
setGeneric("learnMediationDensities", function(object, read_from_rds, filename)
                            standardGeneric("learnMediationDensities"))

# define a method for learning the mediation densities
setMethod("learnMediationDensities", "ImputationPipeline", function(object, read_from_rds, filename) {
  # get the full vector of mediation variables
  mediation_variables <- object@variable_dictionary[["M"]]

  # get the full vector of treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]

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
    adjustment_set <- c(covariate_variables, treatment_variables, adjustment_M)

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

# set the generic for learnTreatmentDensities
setGeneric("learnTreatmentDensities", function(object, read_from_rds, filename)
                                          standardGeneric("learnTreatmentDensities"))

# define a method for learning the treatment densities
setMethod("learnTreatmentDensities", "ImputationPipeline", function(object, read_from_rds, filename) {
  # get the full vector of treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]

  # get the full vector of covariate variables
  covariate_variables <- object@variable_dictionary[["X"]]

  # get the number of mediation variables
  treatment_length <- length(treatment_variables)

  # declare an empty list to save the learned models in
  models <- list()

  # learn a density for each of the mediation variables
  for (index in 1:treatment_length) {
    # get current M we want to learn the density for
    cur_A <- treatment_variables[index]

    # chain rule the rest of the densities
    adjustment_A <- treatment_variables[(index+1):treatment_length]

    # if this is the final density to learn, there is nothing to chain rule
    if (index == treatment_length) adjustment_A <- c()

    # define the covariate set
    adjustment_set <- c(covariate_variables, adjustment_A)

    # learn the density from the mediation data
    cur_model <- learn_density(read_from_rds, cur_A, adjustment_set, object@combined_imputed_data, filename)

    # save the learned model to the list
    models[[index]] <- cur_model
  }

  # save the mediation models to the slot
  object@treatment_models <- models

  # return the modified pipeline object
  return(object)
})

# set the generic for learnTreatmentDensities
setGeneric("learnMarginalTreatmentDensities", function(object, read_from_rds, filename)
                                          standardGeneric("learnMarginalTreatmentDensities"))

# define a method for learning the treatment densities
setMethod("learnMarginalTreatmentDensities", "ImputationPipeline", function(object, read_from_rds, filename) {
  # add the string marginal so that the treatment density and marginal treatment density
  # file names do not conflict
  filename <- paste0("marginal_", filename)

  # get the full vector of treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]

  # get the number of mediation variables
  treatment_length <- length(treatment_variables)

  # add a column of 1's so that the final treatment does not
  # have an empty adjustment set
  object@combined_imputed_data$const <- 1

  # declare an empty list to save the learned models in
  models <- list()

  # learn a density for each of the mediation variables
  for (index in 1:treatment_length) {
    # get current M we want to learn the density for
    cur_A <- treatment_variables[index]

    # chain rule the rest of the densities
    adjustment_A <- treatment_variables[(index+1):treatment_length]

    # if this is the final density to learn, there is nothing to chain rule
    if (index == treatment_length) adjustment_A <- c("const")

    # define the covariate set
    adjustment_set <- c(adjustment_A)

    # learn the density from the mediation data
    cur_model <- learn_density(read_from_rds, cur_A, adjustment_set, object@combined_imputed_data, filename)

    # save the learned model to the list
    models[[index]] <- cur_model
  }

  # save the mediation models to the slot
  object@treatment_marginal_models <- models

  # return the modified pipeline object
  return(object)
})

# set the generic for imputeMediators
setGeneric("imputeMediators", function(object, seed) standardGeneric("imputeMediators"))

# define a method for imputing the mediators in the primary data from the learned densities in 
# the mediators data
setMethod("imputeMediators", "ImputationPipeline", function(object, seed) {
  # take the primary dataset and create columns of NA for each of the mediators
  primary_data <- object@primary_data

  # get the list of mediator variables
  mediation_variables <- object@variable_dictionary[["M"]]

  # create a column of missing values for each mediator
  for (index in 1:length(mediation_variables)) {
    primary_data[[mediation_variables[index]]] <- NA
  }

  # combine the primary and mediation datasets, keeping the mediators in the
  # primary dataset as NAs
  combined_data <- base::rbind(primary_data, object@mediation_data)

  # impute missing values using a random forest model one time
  imputed_data <- mice(combined_data, method="rf", m=1, seed=seed)

  # get the completed dataset from the first (and only) imputation
  completed_data <- complete(imputed_data, 1)

  # save the new dataframe to the object
  object@combined_imputed_data <- data.frame(completed_data)

  return(object)
})

# set the generic for estimating the mixed mediation term
setGeneric("computePseudoOutcome", function(object, a_prime_vals, a_vals) standardGeneric("computePseudoOutcome"))

# define a method for estimating the mediation term
setMethod("computePseudoOutcome", "ImputationPipeline", function(object, a_prime_vals, a_vals) {
  # compute the pseudo-outcome
  Y_p <- object@combined_imputed_data[[object@variable_dictionary[["Y"]]]]

  # get the covariate variables
  covariate_variables <- object@variable_dictionary[["X"]]

  # get the treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]
  treatment_length <- length(treatment_variables)

  # get the mediation variables and its length
  mediation_variables <- object@variable_dictionary[["M"]]
  med_length <- length(mediation_variables)

  # get the models for the mediation analysis
  mediation_models <- object@mediation_models

  # get the number of rows in the dataset we are working with
  row_n <- nrow(object@combined_imputed_data)

  for (i in 1:med_length) {
    # get the current model, which is stored in the first position of the tuple
    cur_model <- mediation_models[[i]][[1]]

    # get the adjustment set for this M
    adjustment_M <- mediation_variables[(i+1):med_length]
    # if this is the final density to learn, there is nothing to chain rule
    if (i == med_length) adjustment_M <- c()
    # define the covariate set
    adjustment_set <- c(covariate_variables, treatment_variables, adjustment_M)
  
    # make a copy of the dataset
    data_copy_prime <- data.frame(object@combined_imputed_data)
    # iterate through all treatment variables and change each treatment
    # value to its prime value
    for (j in 1:length(treatment_variables)) {
      data_copy_prime[[treatment_variables[j]]] <- rep(a_prime_vals[j], row_n)
    }

    # save the outcome and covariate variables to make predictions for
    A <- data_copy_prime[[mediation_variables[i]]]
    L <- data_copy_prime[, adjustment_set ]
    
    # the prime values are in the denominator
    denominator <- predict(cur_model, new_A=A, new_W=L)
    
    # make a copy of the dataset
    data_copy <- data.frame(object@combined_imputed_data)
    # iterate through all treatment variables and change each treatment
    # value to its a value
    for (j in 1:length(treatment_variables)) {
      data_copy[[treatment_variables[j]]] <- rep(a_vals[j], row_n)
    }
    
    # save the outcome and covariate variables to make predictions for
    A <- data_copy[[mediation_variables[i]]]
    L <- data_copy[, adjustment_set ]
    
    # the non-prime values are in the numerator
    numerator <- predict(cur_model, new_A=A, new_W=L)
    
    # update the value of the pseudo-outcome
    Y_p <- Y_p * numerator / denominator
  }

  # save the pseuo-outcome into the slot
  object@Y_p <- Y_p

  # save the pseudo-outcome as a column in the final dataframe
  object@combined_imputed_data[["Y_p"]] <- Y_p

  # return the object
  return(object)
})

# set the generic for computing the MSM weights
setGeneric("computeMSMWeights", function(object) standardGeneric("computeMSMWeights"))

setMethod("computeMSMWeights", "ImputationPipeline", function(object) {
  # get the size of the dataset
  row_n <- nrow(object@combined_imputed_data)

  # get the treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]
  treatment_length <- length(treatment_variables)

  # compute MSM weights
  weights <- rep(1, row_n)
  # iterate through the treatment models
  for (i in 1:treatment_length) {
    # get the estimated densities for the numerator, which are in index 2 of the tuple
    numerator <- object@treatment_marginal_models[[i]][[2]]
    # get the estimated densities for the denominator, which are in index 2 of the tuple
    denominator <- object@treatment_models[[i]][[2]]
    # update the weights
    weights <- weights * numerator / denominator
  }

  # standardize weights
  weights_stand <- weights / sum(weights) * row_n

  # store the weights into the slot
  object@MSM_weights <- weights_stand

  # return the object
  return(object)
})

# set the generic for estimating the mediation term
setGeneric("estimateMediationTerm", function(object, a_prime_vals) standardGeneric("estimateMediationTerm"))

# define method for estimating the mediation term
setMethod("estimateMediationTerm", "ImputationPipeline", function(object, a_prime_vals) {
  # retrieve the combined and imputed dataset
  data <- object@combined_imputed_data

  # save the weights as a column in the dataset
  data$MSM_weights <- object@MSM_weights

  # create a formula for predicting the pseudo-outcome
  formula <- create_formula("Y_p", object@variable_dictionary[["A"]], two_way=TRUE)

  # train a linear model for the pseudo-outcome given the data and weights
  model <- glm(formula, family=gaussian, data=data, weights=MSM_weights)
  print(summary(model))

  # make a copy of the combined_imputed_data
  data_copy <- data.frame(object@combined_imputed_data)

  # get the size of the dataset
  row_n <- nrow(data_copy)

  # get the treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]
  treatment_length <- length(treatment_variables)

  # replace each of the A's with the a_prime_value
  for (i in 1:length(treatment_variables)) {
    data_copy[[treatment_variables[i]]] <- rep(a_prime_vals[i], row_n)
  }

  # use the model to make predictions using the copied dataset
  predictions <- predict(model, newdata=data_copy, type="response")

  # the estimate is the mean of the predictions
  object@mediation_term <- mean(predictions)

  return(object)
})

# set the generic for estimating counterfactual terms
setGeneric("estimateCounterfactual", function(object, a_vals, prime) standardGeneric("estimateCounterfactual"))

# define method for estimating counterfactual terms
setMethod("estimateCounterfactual", "ImputationPipeline", function(object, a_vals, prime) {
  # retrieve the combined and imputed dataset
  data <- object@combined_imputed_data

  # save the weights as a column in the dataset
  data$MSM_weights <- object@MSM_weights

  # create a formula for predicting the pseudo-outcome
  formula <- create_formula(object@variable_dictionary[["Y"]][1], object@variable_dictionary[["A"]], two_way=TRUE)

  # train a linear logistic model for the pseudo-outcome given the data and weights
  # the regular outcome is binary
  model <- glm(formula, family=binomial, data=data, weights=MSM_weights)

  # make a copy of the combined_imputed_data
  data_copy <- data.frame(object@combined_imputed_data)

  # get the size of the dataset
  row_n <- nrow(data_copy)

  # get the treatment variables
  treatment_variables <- object@variable_dictionary[["A"]]
  treatment_length <- length(treatment_variables)

  # replace each of the A's with the a_vals
  for (i in 1:length(treatment_variables)) {
    data_copy[[treatment_variables[i]]] <- rep(a_vals[i], row_n)
  }

  # use the model to make predictions using the copied dataset
  predictions <- predict(model, newdata=data_copy, type="response")

  # the estimate is the mean of the predictions
  # if prime is true then save the estimate to the prime slot,
  # otherwise save the estimate to the non-prime slot
  if (prime) object@counterfactual_a_prime <- mean(predictions)
  else object@counterfactual_a <- mean(predictions)

  return(object)
})

# Simulation study starts here
# set the seed
set.seed(0)

# define the size of the primary dataset
primary_n <- 1000
# define the size of the secondary dataset
secondary_n <- 150

# define the standard deviation
sd = 1

# generate variables for the primary dataset
X_primary <- rnorm(primary_n, 0, sd)
A_primary <- X_primary + rnorm(primary_n, 0, sd)
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
X_secondary <- rnorm(secondary_n, 0, sd)
A_secondary <- X_secondary + rnorm(secondary_n, 0, sd)
M_secondary <- A_secondary + X_secondary + rnorm(secondary_n, 0, sd)
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
pipeline <- imputeMediators(pipeline, 0)

# learn the treatment densities and update the object
pipeline <- learnTreatmentDensities(pipeline, TRUE, "test")

# learn the marginal treatment densities and update the object
pipeline <- learnMarginalTreatmentDensities(pipeline, TRUE, "test")

# estimate the mediation term
pipeline <- computePseudoOutcome(pipeline, c(1), c(-1))

# compute the MSM weights
pipeline <- computeMSMWeights(pipeline)

# estimate the mediation term
pipeline <- estimateMediationTerm(pipeline, c(1))
print(pipeline@mediation_term)

# estimate the counterfactual terms
pipeline <- estimateCounterfactual(pipeline, c(1), TRUE)
print(pipeline@counterfactual_a_prime)

pipeline <- estimateCounterfactual(pipeline, c(-1), FALSE)
print(pipeline@counterfactual_a)

#####
# compute ground-truth values
X_primary <- rnorm(primary_n, 0, sd)
A_primary <- X_primary + rnorm(primary_n, 0, sd)
M_primary <- 1 + X_primary + rnorm(primary_n, 0, 0.1*sd)
Y_primary <- rbinom(primary_n, 1, expit(1 + X_primary + M_primary))
print(paste("Y(1)", mean(Y_primary)))

X_primary <- rnorm(primary_n, 0, sd)
A_primary <- X_primary + rnorm(primary_n, 0, sd)
M_primary <- -1 + X_primary + rnorm(primary_n, 0, 0.1*sd)
Y_primary <- rbinom(primary_n, 1, expit(-1 + X_primary + M_primary))
print(paste("Y(-1)", mean(Y_primary)))

X_primary <- rnorm(primary_n, 0, sd)
A_primary <- X_primary + rnorm(primary_n, 0, sd)
M_primary <- -1 + X_primary + rnorm(primary_n, 0, 0.1*sd)
Y_primary <- rbinom(primary_n, 1, expit(1 + X_primary + M_primary))
print(paste("Y(1, M(-1))", mean(Y_primary)))

# show(pipeline)