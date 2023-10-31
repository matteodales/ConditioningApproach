library(tidyverse)
library(gamm4)
library(lmerTest)
#library(selfmade)
library(pracma)
library(glmmLasso)
library(base)
library(glmnet)
library(MASS)
library(nlme)

setwd("C:/Users/matteda/OneDrive - Universitetet i Oslo/Skrivebord/phd/ConditioningApproach/selfmade/")
files.sources = list.files()
invisible(sapply(files.sources, source))
setwd("C:/Users/matteda/OneDrive - Universitetet i Oslo/Skrivebord/phd/ConditioningApproach/src/")


# returns the average length of confidence intervals in the results of selective inference

ci_length <- function(res){

  len <- length(res$selinf)
  not_na <- 0
  avg <- 0

  for(i in 1:len){
    ciu <- as.numeric(res$selinf[[i]]['ciu'])
    if(is.na(ciu)) next
    cil <- as.numeric(res$selinf[[i]]['cil'])
    if(is.na(cil)) next
    not_na <- not_na + 1
    avg <- avg + (ciu - cil)
  }

  if(not_na>0) return(avg/not_na)
  else return(NA)
  
}

# function that returns the boolean vector of selected variables after selective inference
# using Benjamini Hockberg procedure

selection_with_selinf <- function(res, sel_without_selinf, fdr_level = 0.1){

    # sel_vec is the boolean vector of selected variables after lasso

    sel <- c(rep(0,length(sel_without_selinf)))
    names(sel) <- names(sel_without_selinf)
    only_sel <- names(sel_without_selinf[sel_without_selinf == 1])

    p_vals <- c(rep(0, length(res$selinf)))

    for(i in 1:length(res$selinf)){
        p_vals[i] <- res$selinf[[i]]['pval']
    }

    #p_vals <- p.adjust(p_vals, method='BH')
    
    for(i in 1:length(res$selinf)){
        if(p_vals[i]<fdr_level) sel[only_sel[i]] <- 1
        }

    return(sel)
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

# adapting the functions in selfmade package
# now we assume the variables that have a random effect, to be always included in the model (and excluded from penalization)
# these variables are also used in the reduced model that computes the estimate of variance


vcov_RI <- function(mod, sigma = getME(mod, "sigma"), tau = getME(mod, "theta"),
                    sigmaT = NULL, tauT = NULL)
{
  
  if(!is.null(sigmaT)) sigma <- sigmaT
  if(!is.null(tauT)) tau <- tauT
  n <- NROW(mod@resp$y)
  nrsubj <- nlevels(mod@flist[[1]])
  getME(mod,'Z')%*%(sigma(mod)^2*getME(mod,'Lambda')%*%getME(mod,'Lambdat'))%*%t(getME(mod,'Z'))+ diag(sigma(mod)^2,n)
  
}
