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
    charForm1 <- character(length(cnms))
    charForm2 <- character(length(cnms))
    
    for(i in 1:length(cnms)) {
      if (cnms[[i]][1] == "(Intercept)") {
        cnms[[i]][1] <- "1"
      } else {
        tpv <- cnms[[i]]
        cnms[[i]] <- append("",tpv)
        cnms[[i]][1] <- "-1"
      }
      charForm1[i] <- paste("(", paste(cnms[[i]], collapse = " + "),
                           " | ",names(cnms)[i], ")", sep = "")
      charForm2[i] <- paste(paste(cnms[[i]], collapse = " + "))
    }
    
    reFormula <- paste(c(charForm2, charForm1), collapse = " + ")
    
    return(reFormula)
  }


minMod <- function(mod)
{
  lhs <- formula(mod)[[2]]
  reFormula <- cnms2formula(mod@cnms)

  #modIC <-
  update(mod, formula. = reformulate(c("1", reFormula), lhs))
  # return(list(sigma = getME(modIC, "sigma"),
  #             tau = getME(modIC, "theta"))
  # )
  
}

