library(tidyverse)
library(gamm4)
library(lmerTest)
library(selfmade)
library(pracma)
library(glmmLasso)
library(base)
library(glmnet)
library(MASS)
library(nlme)


## generating the compound simmetry covariance matrix

comp_cor <- function(n,rho){
    mat <- matrix(rho,n,n)
    diag(mat) <- rep(1,n)
    return(mat)
}

## generating the autocorrelation covariance matrix

ar1_cor <- function(n, rho) {
exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
    (1:n - 1))
rho^exponent
}

## data generating function

data_generator <- function(n_subjects, n_observations, p, q, SNR, prop_relevant, rho){


    # p: number of covariates
    # q: number of random effects (other than the intercept)
    # SNR: signal-to-noise ratio
    # prop_relevant: proportion of covariates with non-zero coefficients

    n <- n_subjects*n_observations
    subjects = as.factor(gl(n_subjects,n_observations))

    # X <- matrix(0, nrow = n, ncol = p)
    # for(i in 0:(n_subjects-1)) X[(n_observations*i+1):(n_observations*(i+1)),] <- t(mvrnorm(p, rep(0, n_observations), comp_cor(n_observations, rho)))

    cov <- bdiag(list(comp_cor(n_observations, rho))[rep(1,n_subjects)])
    X <- t(mvrnorm(p, rep(0, n), cov))

    beta_0 = 1

    n_rel = floor(p*prop_relevant)

    beta = c(sample (c(1,-1), size=n_rel, replace=T), rep(0, p-n_rel))

    etaFix <- beta_0 + X%*%beta #fixed component of the model
    signal <- sd(etaFix)
    sd = signal / SNR

    # the random variables are a subset of the fixed
    # we consider the first q variables to have a random effect as well
    Z <- c()
    if (q>0) Z <- X[,1:q]

    # Sigma = positive definite covariance matrix for random effects
    Sigma = diag(q+1)
    # bi = realizations of random effects
    bi = MASS::mvrnorm(n_subjects, rep(0,q+1), Sigma)

    e = rnorm(n, 0, sd)

    y = etaFix + rowSums(bi[subjects,]*cbind(rep(1,n),Z)) + e

    return(list(X = X, Z=Z, subjects = subjects, y = y, beta=beta, sd=sd, cov=cov))

}

## function to return metrics (True Positive Rate, False Discovery Rate)

metrics <- function(selected, true){

    tpr <- sum(selected*true)/sum(true)
    fdr <- sum(selected==1 & true==0)/sum(selected)

    return(list(tpr = tpr, fdr = fdr))
}


# from rugamer, generates estimated cov mat
vcov_RI <- function(mod, sigma = getME(mod, "sigma"), tau = getME(mod, "theta"),
                    sigmaT = NULL, tauT = NULL)
{
  
  if(!is.null(sigmaT)) sigma <- sigmaT
  if(!is.null(tauT)) tau <- tauT
  n <- NROW(mod@resp$y)
  nrsubj <- nlevels(mod@flist[[1]])
  bdiag(list(diag(n/nrsubj)*sigma^2 + tau^2)[rep(1,nrsubj)])
  
}


# copied from cAIC4 package
cnms2formula <-
  function(cnms) {
    # A function that builds a random effects formula from the ?component names?,
    # a list that can be extracted from an lmerMod object by .@cnms or
    # getME(., "cnms").
    #
    # Args:
    #   cnms     = List from an lmerMod object by .@cnms or getME(., "cnms").
    #
    # Returns:
    #   reFormula = random effects part of a lmerMod formula
    #
    len      <- unlist(lapply(cnms, length))
    cnms     <- cnms[which(len != 0)]
    charForm <- character(length(cnms))
    
    for(i in 1:length(cnms)) {
      if (cnms[[i]][1] == "(Intercept)") {
        cnms[[i]][1] <- "1"
      } else {
        tpv <- cnms[[i]]
        cnms[[i]] <- append("",tpv)
        cnms[[i]][1] <- "-1"
      }
      charForm[i] <- paste("(", paste(cnms[[i]], collapse = " + "),
                           " | ",names(cnms)[i], ")", sep = "")
    }
    
    reFormula <- paste(charForm, collapse = " + ")
    
    return(reFormula)
  }
