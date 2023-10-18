mocasin <- function(
  mod,
  checkFun,
  this_y = NULL,
  nrSamples = 1000,
  bayesian = TRUE,
  varInTestvec = c("est", "minMod", "varY", "supplied"),
  varForSampling = c("est", "minMod", "varY", "supplied"),
  VCOV_vT = NULL,
  VCOV_sampling = NULL,
  conditional  = TRUE,
  name = NULL,
  nrlocs = 7,
  complete_effect = NULL,
  which = NULL,
  vT = NULL,
  G = NULL,
  efficient = TRUE,
  trace = TRUE
)
{

  ######### get further setting specitivities #########
  # check type
  isMM <- or(inherits(x = mod, what = "merMod"),inherits(x = mod, what = "lme"))
  isLM <- any(class(mod)=="lm")
  isGAMM <- any(class(mod)=="gamm")

  # get response
  if(is.null(this_y)){
    if(isMM) this_y <- mod@resp else if(isLM)
      this_y <- mod$y else
        if(isGAMM) this_y <- mod$lme$data$y else
        this_y <- mod$mer@resp$y
  }

  # match arguments for variance
  varInTestvec = match.arg(varInTestvec)
  varForSampling = match.arg(varForSampling)

  # type of VCOV
  diagCOV <- (conditional & isMM) | isLM | !isMM

  # set bayesian to FALSE for marginal case
  if(!conditional & isMM) bayesian <- FALSE
  #####################################################

  ######### define further models and params ##########
  if(any("minMod" %in% c(varInTestvec, varForSampling))){

    if(!isLM){ # estimate IC mod
      if(isMM) modIC = minMod(mod) else
        modIC = minMod(mod$mer)
    }

    if(isLM){

      sigmaIC2 <- sigma(mod)^2
      tauIC2 <- NULL

    }else{


      sigmaIC2 = sigma(modIC)^2
      tauIC2 = getME(modIC, "theta")

    }

  }


  # define #obs
  n <- length(this_y)

  # estimated variance
  if(isLM | isMM) sigma2 <- sigma(mod)^2 else
    sigma2 <- mod$gam$sig2

  # plugin by Tibshirani et al. 2018
  sigma2_y <- var(this_y)*(n-1)/n

  #####################################################

  ############## define (co-)variances ################

  if(varInTestvec=="supplied" &
     is.null(VCOV_vT))
      stop("Must supply VCOV_vT if varInTestvec=='supplied'")


  if(varForSampling=="supplied" &
    is.null(VCOV_sampling))
      stop("Must supply VCOV_sampling if varForSampling=='supplied'")

  # if no supplied variance for testvector, define

  sigma2_sampling <- NULL

  if(is.null(VCOV_sampling)){

    if(conditional){

      VCOV_sampling <- switch (varForSampling,
                               est = sigma2,
                               varY = sigma2_y,
                               minMod = sigmaIC2
      )

      VCOV_sampling = VCOV_sampling * diag(rep(1,n))

    }else{

      if(varForSampling == "varY")
        stop("Option 'varY' not possible for marginal mixed model inference.")

      VCOV_sampling <- switch(varForSampling,
                              est = mod$sigma^2 * diag(rep(1,n)),
                              minMod = vcov_RI(modIC)
      )

      sigma2_sampling <- switch(varForSampling,
                                varY = sigma2_y,
                                est = mod$sigma^2,
                                minMod = modIC$sigma^2)

    }
  }

  if(is.null(sigma2_sampling)) sigma2_sampling = VCOV_sampling[1,1]

  # if no supplied variance for testvector, define
  if(is.null(VCOV_vT)){

    if(conditional){

      VCOV_vT <- switch (varInTestvec,
                         est = NULL,
                         varY = sigma2_y,
                         minMod = sigmaIC2
      )

    }else{ # marginal



    }

  }

  #####################################################

  #### define testvector corresponding to settings ####

  if(!is.null(which)) wn <- which else wn <- name

  if(isMM){

    if(conditional){

      vT <- testvec_for_mm(mod,
                           marginal = !conditional,
                           G = G,
                           sig2 = VCOV_vT,
                           which = wn)

    }else{ # marginal


      # use VCOV only if supplied, else more efficient matrix inversion
      # in the testvec function
      vT <- testvec_for_mm(mod,
                           marginal = !conditional,
                           VCOV = VCOV_vT,
                           sig2 = sigma2_sampling,
                           which = wn,
                           efficient = efficient)



    }

  }else{ # additive case

    if(is.null(wn)) wn <- attr(mod$gam$terms, "term.labels")
    if(is.null(complete_effect)){ 
      complete_effect <- rep(FALSE, length(wn))
      names(complete_effect) <- wn
    }
    
    if(any(complete_effect) & length(name) == 1) {
      names(complete_effect) <- name
    }
    
    if(length(name) == 2) {
      if(!any(complete_effect)) stop("Single point testing only for univariate effects.")
      wn <- paste(name, collapse = ",")
      names(complete_effect) <- paste(wn, collapse = ",")
    }
    
    vT <- sapply(wn, function(name)
      res <- testvec_for_gamm4(mod,
                               name = name,
                               sigma2 = VCOV_vT,
                               nrlocs = nrlocs,
                               complete_effect = complete_effect[[name]])
    )
    
    if(is.list(vT[[1]])) vT <- unlist(vT, recursive = FALSE)
    
  }
  #####################################################

  ############## compute selective inference ##########
  pbi <- 0
  selinf <- lapply(1:length(vT), function(j){

    vt = vT[[j]]
    if(any(complete_effect)) 
      ce = complete_effect[[which(sapply(wn, function(x) grepl(x, names(vT)[j])))]] else
      ce = FALSE
    if(trace) cat("Computing inference for variable (location) ", pbi+1, "\n\n")
    res <- pval_vT_cov(vT = vt,
                       VCOV = VCOV_sampling,
                       this_y = this_y,
                       nrSamples = nrSamples,
                       checkFun = checkFun,
                       bayesian = bayesian,
                       trace = trace,
                       complete_effect = ce
    )
    if(trace) cat("\n\n")
    pbi <<- pbi + 1
    return(res)
  }
  )
  # if(!isMM) names(selinf) <- wn
  
  ######################################################

  retl <- list(vT = vT, selinf = selinf)
  class(retl) <- "selfmade"
  return(retl)

}

#' @title Generic methods for selfmade objects
#'
#' @description Generic methods which can be used for objects
#' fitted with the \code{mocasin} function
#'
#' @param x selfmade object
#' @param ... further arguments, currently unused.
#'
#' @method print selfmade
#' @rdname methodsSelfmade
#'
print.selfmade <- function(x, ...)
{

  covariate <- names(x$selinf)
  df <- do.call("rbind", x$selinf)
  df$covariate <- covariate
  print(df)
  invisible(df)

}

testvec_for_mm <- function(
  mod,
  which = NULL,
  marginal = TRUE,
  VCOV = NULL,
  G = NULL,
  sig2 = NULL,
  efficient = TRUE
)
{
  
  if(class(mod)=="lm"){
    
    X <- model.matrix(mod)
    vTs <- solve(crossprod(X))%*%t(X)
    if(is.null(which)) which <- 1:ncol(X)
    return(lapply(which, function(j) t(vTs[j,])))
    
  }
  if(is.null(sig2)) sig2 <- sigma(mod)^2
  
  Z       <- getME(mod, "Z")
  X       <- getME(mod, "X")
  n       <- nrow(X)
  Lambda <- getME(mod, "Lambda")
  Lambdat <- getME(mod, "Lambdat")
  
  
  if(!marginal){
    
    if(length(mod@Gp)>3) stop("Please implement for more than 2 REs.")
    C <- cbind(X,Z)
    if(is.null(G) || dim(G)[2] != dim(getME(mod,"ST")[[1]])[2]){
      A <- bdiag(list(matrix(0, ncol=ncol(X), nrow=ncol(X)), Lambda%*%Lambdat)) 
    }else{
      if(all(c(G)==0))
        solveG <- G else solveG <- solve(G)
        A <- bdiag(list(matrix(0, ncol=ncol(X), nrow=ncol(X)),
                        bdiag(list(solveG)[rep(1,ncol(Z)/ncol(G))])))
    }
    if(class(A)=="try-error") # G does not fit to the selected model
    {
      A <- bdiag(list(matrix(0, ncol=ncol(X), nrow=ncol(X)), Lambda%*%Lambdat))
      warning("Could not create A, using the estimated covariance")
    }
    if(is.null(VCOV)) VCOV <- Matrix::crossprod(C)/sig2 + A
    
  }
  
  
  if(marginal){
    
    L       <- getME(mod, "L")
    if(is.null(VCOV)){
      
      V0inv   <- diag(rep(1, n)) - crossprod(solve(L, system = "L") %*%
                                               solve(L, Lambdat, system = "P") %*% t(Z))
      # = solve(diag(n)+Z%*%Lambda%*%Lambdat%*%t(Z))
      V0inv <- V0inv / sig2
      
    }else{
      
      V0inv <- solve(VCOV)
      
    }
    
  }
  
  if(marginal){
    if(efficient)
      vTs <- solve(crossprod(X, V0inv%*%X))%*%t(X)%*%V0inv else
        vTs <- solve(crossprod(X))%*%t(X)
  }else{ # conditional
    V0inv <- solve(VCOV)
    vTs <- V0inv%*%Matrix::t(C)/sig2
  }
  
  if(is.null(which)) if(marginal) which <- 1:ncol(X) else which <- 1:ncol(C)
  vTs <- lapply(which, function(j){
    
    vt <- vTs[j,]
    if(!marginal) attr(vt, "var") <- V0inv[j,j]
    return(vt)
    
  })
  return(lapply(vTs,t))
  
}