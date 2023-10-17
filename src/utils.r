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

    adj_p_vals <- p.adjust(p_vals, method='BH')
    
    for(i in 1:length(res$selinf)){
        if(adj_p_vals[i]<fdr_level) sel[only_sel[i]] <- 1
        }

    return(sel)
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
