library('rlist')

## generating data for a linear mixed model

data_generator <- function(n_subjects, n_observations, p, q, SNR = 2, prop_relevant = 0.2, rho = 0.5, Sigma = NULL){


    # p: number of covariates
    # q: number of random effects (other than the intercept)
    # SNR: signal-to-noise ratio
    # prop_relevant: proportion of covariates with non-zero coefficients
    # rho: parameter for the compound simmetry covariance matrix
    # Sigma: covariance matrix for the random effects, defaults to identity matrix

    n <- n_subjects*n_observations
    subjects = as.factor(gl(n_subjects,n_observations))

    covMat <- bdiag(list(comp_cor(n_observations, rho))[rep(1,n_subjects)])
    X <- t(mvrnorm(p, rep(0, n), covMat))
    colnames(X) <- paste("X", 1:p, sep="")

    beta_0 = 1

    n_rel = floor(p*prop_relevant)

    possible_coefs = c(1, 0.8, -0.8,-1)

    beta = c(sample (possible_coefs, size=n_rel, replace=T), rep(0, p-n_rel))

    etaFix <- beta_0 + X%*%beta #fixed component of the model
    signal <- sd(etaFix)
    sd = signal / SNR

    # the random variables are a subset of the fixed
    # we consider the first q variables to have a random effect as well
    Z <- c()
    if (q>0) Z <- X[,1:q]

    # Sigma = positive definite covariance matrix for random effects
    if(is.null(Sigma)) Sigma <- diag(q+1)
    # bi = realizations of random effects
    bi = MASS::mvrnorm(n_subjects, rep(0,q+1), Sigma)

    e = rnorm(n, 0, sd)

    y = etaFix + rowSums(bi[subjects,]*cbind(rep(1,n),Z)) + e

    return(list(X = X, Z=Z, subjects = subjects, y = y, beta=beta, sd=sd, covMat=covMat))

}






# generating data with time random effect

data_generator_random_time <- function(n_subjects, n_observations, p, SNR = 2, prop_relevant = 0.2, rho = 0.5, Sigma = NULL){


    # p: number of covariates
    # q: number of random effects --> equal to n_observations - 1
    # SNR: signal-to-noise ratio
    # prop_relevant: proportion of covariates with non-zero coefficients
    # rho: parameter for the compound simmetry covariance matrix
    # Sigma: covariance matrix for the random effects, defaults to identity matrix

    n <- n_subjects*n_observations
    q <- n_observations-1
    subjects = as.factor(gl(n_subjects,n_observations))

    covMat <- bdiag(list(ar1_cor(n_observations, rho))[rep(1,n_subjects)])
    X <- t(mvrnorm(p, rep(0, n), covMat))
    colnames(X) <- paste("X", 1:p, sep="")

    beta_0 = 1

    n_rel = floor(p*prop_relevant)

    possible_coefs = c(1, 0.8, -0.8,-1)

    beta = c(sample (possible_coefs, size=n_rel+q, replace=T), rep(0, p-n_rel))

    # the random variables are the intercept and dummies identifying the time (numbered observation)
    Z <- list.rbind(list(rbind(rep(0,n_observations-1),diag(n_observations-1)))[rep(1,n_subjects)])
    colnames(Z) <- paste("t", 1:q, sep="")
    X <- cbind(Z,X)

    etaFix <- beta_0 + X%*%beta #fixed component of the model
    signal <- sd(etaFix)
    sd = signal / SNR

    # Sigma = positive definite covariance matrix for random effects
    if(is.null(Sigma)) Sigma <- diag(q+1)
    # bi = realizations of random effects
    bi = MASS::mvrnorm(n_subjects, rep(0,q+1), Sigma)

    e = rnorm(n, 0, sd)

    y = etaFix + rowSums(bi[subjects,]*cbind(rep(1,n),Z)) + e

    return(list(X = X, Z=Z, subjects = subjects, y = y, beta=beta, sd=sd, covMat=covMat))

}