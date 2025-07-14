# use object oriented S4 programming in R to set up an imputation pipeline that
# we can use more easily with the data application

# library for learning continuous densities
library(haldensify)

# in this file, we will define the class, methods, attributes, etc. and also
# run a simulation study using the class to demonstrate how to use it

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
                    treatment_models="ANY" # a list of treatment models learned from the haldensify package
                    ))

# create a constructor
NewImputationPipeline <- function(primary_data, mediation_data, variable_dictionary) {
  new("ImputationPipeline", primary_data=primary_data, 
                mediation_data=mediation_data,
                variable_dictionary=variable_dictionary,
                treatment_models=NULL)
}

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
M_primary <- A_primary + X_primary + rnorm(secondary_n, 0, 0.1*sd)
Y_primary <- rbinom(primary_n, 1, expit(A_primary + X_primary + M_primary))

# save the primary variables into a dataframe
primary_data <- data.frame(
  X=X_primary,
  A=A_primary,
  Y=Y_primary
)

# generate variables for the secondary dataset
X_secondary <- rnorm(secondary_n, 0, 1)
A_secondary <- X + rnorm(secondary_n, 0, 0.1*sd)
M_secondary <- A + X + rnorm(secondary_n, 0, 0.1*sd)
Y_secondary <- rbinom(secondary_n, 1, expit(A + X + M))

# save the variables with the mediation data into the secondary
# dataset
mediation_data <- data.frame(
  X=X_secondary,
  A=A_secondary,
  Y=Y_secondary
)

pipeline <- NewImputationPipeline(primary_data, mediation_data)

print(pipeline)