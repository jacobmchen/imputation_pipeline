library(dbarts)

# define an expit function
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

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

# experiments to test if an estimation strategy for the mediation term works
if (sys.nframe() == 0) {
    # set up a simple DGP
    set.seed(0)

    n <- 5000

    X <- rnorm(n, 0, 1)
    A <- rbinom(n, 1, expit(X))
    M <- rbinom(n, 1, expit(X + A))
    Y <- X + A + M + rnorm(n, 0, 1)

    df <- data.frame(
        X=X,
        A=A,
        M=M,
        Y=Y
    )

    # compute ground truth values
    M_a <- rbinom(n, 1, expit(X + 0))
    Y_a <- X + 0 + M_a + rnorm(n, 0, 1)
    print(paste("E[Y(a)]", mean(Y_a)))

    M_a_prime <- rbinom(n, 1, expit(X + 1))
    Y_a_prime <- X + 1 + M_a_prime + rnorm(n, 0, 1)
    print(paste("E[Y(a')]", mean(Y_a_prime)))

    Y_mixed <- X + 1 + M_a + rnorm(n, 0, 1)
    print(paste("E[Y(a', M(a))]", mean(Y_mixed)))

    # estimate counterfactuals
    df_a = data.frame(df)
    df_a$A <- rep(0, n)
    # print("E[Y(a)] estimate")
    # print(compute_po_Y(df, "Y", c("A", "X"), rep(1, n), df_a, q = 0.025, ndraws=500))

    df_a_prime = data.frame(df)
    df_a_prime$A <- rep(1, n)
    # print("E[Y(a')] estimate")
    # print(compute_po_Y(df, "Y", c("A", "X"), rep(1, n), df_a_prime, q = 0.025, ndraws=500))

    # estimate mediation term
    # step 1, estimate the probability of M given A and X using a linear logistic regression
    M_model <- glm(M~A+X, family=binomial, df)
    print(summary(M_model))

    M_a_preds <- predict(M_model, newdata=df_a, type="response")
    M_a_prime_preds <- predict(M_model, newdata=df_a_prime, type="response")
    numerator <- df$M * M_a_preds + (1-df$M) * (1-M_a_preds)
    denominator <- df$M * M_a_prime_preds + (1-df$M) * (1-M_a_prime_preds)

    fraction <- numerator / denominator
    
    df$Y_p <- df$Y * fraction

    df_a_prime <- data.frame(df)
    df_a_prime$A <- rep(1, n)

    print("E[Y(a', M(a))] estimate")
    print(compute_po_Y(df, "Y_p", c("A", "X"), rep(1, n), df_a_prime, q=0.025, ndraws=500))
}