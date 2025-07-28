# use object oriented S4 programming in R to set up an imputation pipeline that
# we can use more easily with the data application

# library for learning continuous densities
library(haldensify)

# library for imputing missing values
library(mice)

# library for random forests
library(ranger)

# library for BART
# library(BART)
library(dbarts)

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
                    mediation_term="vector", # a number to store the estimate of the mediation term
                    counterfactual_a="vector", # a number to store the estimate of intervening on a
                    counterfactual_a_prime="vector" # a number to store the estimate of intervening on a'
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
                mediation_term=c(0),
                counterfactual_a=c(0),
                counterfactual_a_prime=c(0))
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

# define a helper function for counterfactuals using BART
# <FUNCTION: compute_po_Y>
# compute E[Y(A=a)] by treating Y as a numeric variable
# <Input>
# data: the dataset in a data.frame format
# Y: string containing the name of the output
# A: string containing the name of the treatment
# weights: vector of MSM weights to adjust for confounding
# data_a: the dataset with values of the treatment set at the 
# counterfactual value
# q: the quantile for which we want to draw the confidence interval
# ndraws: number of MCMC draws for the bart operator
# <Output>
# E[Y(A=a)]
compute_po_Y <- function(data, Y, A, weights, data_a, q = 0.025, ndraws) {
  # get the training features data
  data_copy <- data.frame(data)
  data_copy[[Y]] <- NULL
  x.train <- as.matrix(data_copy[, A])
  
  # get the outcome data
  y.train <- as.matrix(data[[Y]])
  
  # get the testing features data and replacing all values of A with a
  data_copy <- data.frame(data_a)
  data_copy[[Y]] <- NULL
  x.test <- as.matrix(data_copy[, A])
  
  bart_Y <- bart(x.train=x.train, y.train=y.train, x.test=x.test, ndpost=ndraws, 
                  weights=weights, verbose=FALSE)
  
  # calculate the expected value of the potential outcome
  expected_val <- bart_Y$yhat.test.mean
  quants <- as.numeric(quantile(expected_val, c(q, 1-q)))
  bart.expected_val <- mean(expected_val)
  
  return(c(bart.expected_val, quants[1], quants[2]))
}

# define a helper function for counterfactuals using BART
# <FUNCTION: compute_po_Y_binary>
# compute E[Y(A=a)] by treating Y as a binary variable
# <Input>
# data: the dataset in a data.frame format
# Y: string containing the name of the output
# A: string containing the name of the treatment
# weights: vector of MSM weights to adjust for confounding
# data_a: the dataset with values of the treatment set at the 
# counterfactual value
# q: the quantile for which we want to draw the confidence interval
# ndraws: number of MCMC draws for the bart operator
# <Output>
# E[Y(A=a)]
compute_po_Y_binary <- function(data, Y, A, weights, data_a, q = 0.025, ndraws) {
  # get the training features data
  data_copy <- data.frame(data)
  data_copy[[Y]] <- NULL
  x.train <- as.matrix(data_copy[, A])
  
  # get the outcome data
  y.train <- as.matrix(data[[Y]])
  
  # get the testing features data and replacing all values of A with a
  data_copy <- data.frame(data_a)
  data_copy[[Y]] <- NULL
  x.test <- as.matrix(data_copy[, A])
  
  bart_Y <- bart(x.train=x.train, y.train=y.train, x.test=x.test, ndpost=ndraws, 
                  weights=weights, verbose=FALSE)

  # yhat.test is a matrix with dimension ndraws x nrows (ndraws rows and nrows columns),
  # the row number represents the number of the posterior draw
  # the column number represents the row number in the original dataset
  # since the predictor variables in the test data are the same for each row in the original
  # dataset, the estimates
  # across each row are all the same (given a specific posterior draw, predictions will be
  # the same for each row in the original dataset)
  # and we can just take the first column and use its mean as the estimate
  # compute the estimate for the probability of success
  # estimate <- mean(pnorm(bart_Y$yhat.test[,1]))
  
  estimates <- c()
  for (i in 1:nrow(data)) {
    estimates <- c(estimates, mean(pnorm(bart_Y$yhat.test[,i])))
  }
  
  estimate <- mean(estimates)
  quants <- as.numeric(quantile(estimates, c(q, 1-q)))
  
  return(c(estimate, quants))
}

# define a function for obtaining bootstrap confidence intervals for the
# BART estimates
bootstrap_intervals <- function(data, Y, A, weights, data_a, q = 0.025, ndraws,
                                num_bootstraps=200, binary=FALSE) {
  bootstrap_estimates <- c()

  # save the weights into the dataframe
  data$weights <- weights

  # make num_bootstraps iterations
  for (i in 1:num_bootstraps) {
    # get the indices of the bootstrap sample
    bootstrap_sample <- sample(nrow(data), size=nrow(data), replace=TRUE)
    
    # get the bootstrap samples
    bootstrap_data <- data[bootstrap_sample, ]
    bootstrap_data_a <- data_a[bootstrap_sample, ]
    
    # obtain the bootstrap estimate depending on whether the outcome is binary
    if (binary == FALSE) {
      bootstrap_estimates <- c(bootstrap_estimates,
                               compute_po_Y(bootstrap_data, Y, A, bootstrap_data$weights,
                                            bootstrap_data_a, q=q, ndraws)[1])
    } else {
      bootstrap_estimates <- c(bootstrap_estimates,
                               compute_po_Y_binary(bootstrap_data, Y, A, bootstrap_data$weights,
                                            bootstrap_data_a, q=q, ndraws)[1])
    }
  }
  
  # get the bootstrap confidence interval
  quants <- as.numeric(quantile(bootstrap_estimates, c(q, 1-q)))
  
  return(quants)
}

# set the generic for a method that standardizes the primary and mediation data using the
# combined scale of both datasets
setGeneric("standardizeData", function(object, to_exclude) standardGeneric("standardizeData"))

# define a method for standardizing the data
# to_exclude is a vector of variable names that should be excluded from standardization
# such as binary and categorical variables
# note that we do not standardize the Y variable
setMethod("standardizeData", "ImputationPipeline", function(object, to_exclude) {
  # get the A variables and its length
  A_variables <- object@variable_dictionary[["A"]]
  A_length <- length(A_variables)

  # iterate through all of the A variables and standardize them across both datasets
  for (i in 1:A_length) {
    # get the string of the current A
    cur_A <- A_variables[i]

    # skip the excluded variable
    if (cur_A %in% to_exclude) next

    # combine the data for the current A in both datasets
    combined_cur_A <- c(object@primary_data[[cur_A]], object@mediation_data[[cur_A]])

    # compute the scale, which is the standard deviation of the combination of both
    # datasets
    scale <- sd(combined_cur_A)

    # scale the current A variable with the combined standard deviation
    object@primary_data[[cur_A]] <- object@primary_data[[cur_A]] / scale
    object@mediation_data[[cur_A]] <- object@mediation_data[[cur_A]] / scale
  }

  # get the X variables and its length
  X_variables <- object@variable_dictionary[["X"]]
  X_length <- length(X_variables)

  # iterate through all of the X variables and standardize them across both datasets
  for (i in 1:X_length) {
    # get the string of the current X
    cur_X <- X_variables[i]

    # skip the excluded variable
    if (cur_X %in% to_exclude) next

    # combine the data for the current X in both datasets
    combined_cur_X <- c(object@primary_data[[cur_X]], object@mediation_data[[cur_X]])

    # compute the scale
    scale <- sd(combined_cur_X)

    # scale the current X variable
    object@primary_data[[cur_X]] <- object@primary_data[[cur_X]] / scale
    object@mediation_data[[cur_X]] <- object@mediation_data[[cur_X]] / scale
  }

  # get the M variables and its length
  M_variables <- object@variable_dictionary[["M"]]
  M_length <- length(M_variables)

  # iterate through all of the M variables and standardize them across both datasets
  for (i in 1:M_length) {
    # get the string of the current M variable
    cur_M <- M_variables[i]

    # skip the excluded variable
    if (cur_M %in% to_exclude) next

    # only the mediation data has mediation variables, so directly scale it
    # by the observed standard deviation
    object@mediation_data[[cur_M]] <- object@mediation_data[[cur_M]] / sd(object@mediation_data[[cur_M]])
  }

  return(object)
})

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

# function header of computing a potential outcome from BART
# compute_po_Y <- function(data, Y, A, X, data_a, q = 0.025, ndraws)

# define method for estimating the mediation term
setMethod("estimateMediationTerm", "ImputationPipeline", function(object, a_prime_vals) {
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

  bart_estimates <- compute_po_Y(object@combined_imputed_data,
                                 "Y_p",
                                 object@variable_dictionary[["A"]],
                                 object@MSM_weights,
                                 data_copy,
                                 q=0.025,
                                 ndraws=1000)
  
  point_estimate <-bart_estimates[1]
  
  intervals <- bootstrap_intervals(object@combined_imputed_data,
                                   "Y_p",
                                   object@variable_dictionary[["A"]],
                                   object@MSM_weights,
                                   data_copy,
                                   q=0.025,
                                   ndraws=1000,
                                   num_bootstraps=200,
                                   binary=FALSE)
  
  # use adjustment instead
  # bart_estimates <- compute_po_Y(object@combined_imputed_data,
  #                                "Y_p",
  #                                c(object@variable_dictionary[["A"]], object@variable_dictionary[["X"]]),
  #                                rep(1, row_n),
  #                                data_copy,
  #                                q=0.025,
  #                                ndraws=1000)

  object@mediation_term <- c(point_estimate, intervals)

  return(object)
})

# set the generic for estimating counterfactual terms
setGeneric("estimateCounterfactual", function(object, a_vals, prime) standardGeneric("estimateCounterfactual"))

# define method for estimating counterfactual terms
setMethod("estimateCounterfactual", "ImputationPipeline", function(object, a_vals, prime) {
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

  bart_estimate <- compute_po_Y_binary(object@combined_imputed_data,
                                 object@variable_dictionary[["Y"]][1],
                                 object@variable_dictionary[["A"]],
                                 object@MSM_weights,
                                 data_copy,
                                 q=0.025,
                                 ndraws=1000)
  
  intervals <- bootstrap_intervals(object@combined_imputed_data,
                                   object@variable_dictionary[["Y"]][1],
                                   object@variable_dictionary[["A"]],
                                   object@MSM_weights,
                                   data_copy,
                                   q=0.025,
                                   ndraws=1000,
                                   num_bootstraps=200,
                                   binary=TRUE)
  
  # use adjustment instead of weights
  # bart_estimate <- compute_po_Y_binary(object@combined_imputed_data,
  #                                      object@variable_dictionary[["Y"]][1],
  #                                      c(object@variable_dictionary[["A"]], object@variable_dictionary[["X"]]),
  #                                      rep(1, row_n),
  #                                      data_copy,
  #                                      q=0.025,
  #                                      ndraws=1000)

  # the estimate is the mean of the predictions
  # if prime is true then save the estimate to the prime slot,
  # otherwise save the estimate to the non-prime slot
  if (prime) object@counterfactual_a_prime <- c(bart_estimate[1], intervals)
  else object@counterfactual_a <- c(bart_estimate[1], intervals)

  return(object)
})

# this if statement ensures that the code below only runs when this file
# is run directly in the terminal
if (sys.nframe() == 0) {
  # Simulation study starts here
  # set the seed
  set.seed(0)
  
  # define the size of the primary dataset
  primary_n <- 1000
  # define the size of the secondary dataset
  secondary_n <- 1000
  
  # define the standard deviation
  sd = 1
  
  # generate variables for the primary dataset
  X_primary <- rnorm(primary_n, 0, sd)
  A_primary <- X_primary + rnorm(primary_n, 0, sd)
  M_primary <- A_primary + X_primary + rnorm(primary_n, 0, sd)
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

  # standardize the data
  # pipeline <- standardizeData(pipeline, c())
  
  # learn the mediation densities and update the object
  pipeline <- learnMediationDensities(pipeline, FALSE, "test")
  
  # impute M values for the primary dataset
  pipeline <- imputeMediators(pipeline, 0)
  
  # learn the treatment densities and update the object
  pipeline <- learnTreatmentDensities(pipeline, FALSE, "test")
  
  # learn the marginal treatment densities and update the object
  pipeline <- learnMarginalTreatmentDensities(pipeline, FALSE, "test")
  
  # estimate the mediation term
  pipeline <- computePseudoOutcome(pipeline, c(1), c(-1))
  
  # compute the MSM weights
  pipeline <- computeMSMWeights(pipeline)
  
  # estimate the counterfactual terms
  print("counterfactual estimates")
  pipeline <- estimateCounterfactual(pipeline, c(1), TRUE)
  print(pipeline@counterfactual_a_prime)
  
  pipeline <- estimateCounterfactual(pipeline, c(-1), FALSE)
  print(pipeline@counterfactual_a)

  # estimate the mediation term
  print("mediation term estimate")
  pipeline <- estimateMediationTerm(pipeline, c(1))
  print(pipeline@mediation_term)
  
  #####
  # compute ground-truth values
  X_primary <- rnorm(primary_n, 0, sd)
  A_primary <- X_primary + rnorm(primary_n, 0, sd)
  M_primary <- 1 + X_primary + rnorm(primary_n, 0, sd)
  Y_primary <- rbinom(primary_n, 1, expit(1 + X_primary + M_primary))
  print(paste("Y(1)", mean(Y_primary)))
  
  X_primary <- rnorm(primary_n, 0, sd)
  A_primary <- X_primary + rnorm(primary_n, 0, sd)
  M_primary <- -1 + X_primary + rnorm(primary_n, 0, sd)
  Y_primary <- rbinom(primary_n, 1, expit(-1 + X_primary + M_primary))
  print(paste("Y(-1)", mean(Y_primary)))
  
  X_primary <- rnorm(primary_n, 0, sd)
  A_primary <- X_primary + rnorm(primary_n, 0, sd)
  M_primary <- -1 + X_primary + rnorm(primary_n, 0, sd)
  Y_primary <- rbinom(primary_n, 1, expit(1 + X_primary + M_primary))
  print(paste("Y(1, M(-1))", mean(Y_primary)))
}