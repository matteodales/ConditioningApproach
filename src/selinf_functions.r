## function to select the model given a fixed lambda

selFun_fixed_lambda <- function(X, subjects, y, fixed_form, rand_form, lambda){
    
    dat <- data.frame(X, subjects, y)
        
    suppressWarnings(glm1 <- try(glmmLasso(fix = fixed_form, rnd = list(subjects =~ 1), data=dat,
                        , lambda=lambda,switch.NR=FALSE,final.re=FALSE), silent=TRUE))

    return(list(names = names(glm1$coefficients)[glm1$coefficients!=0], vec = glm1$coefficients!=0))
    
    }


selFun_fixed_lambda_randtime <- function(X, Z, subjects, y, lambda){
    
        lasso_res <- lmmlasso(cbind(rep(1,n),X),y,cbind(rep(1,n),Z),subjects,lambda=lambda,nonpen = c(1:4))

        return(list(names = names(lasso_res$coefficients)[lasso_res$coefficients!=0], vec = lasso_res$coefficients!=0))
    
    }

## function to select the model by choosing the optimal lambda through BIC

selFun_adapting_lambda <- function(X, subjects, y, fixed_form, rand_form, lambda.max.min.ratio = 0.01, n_lambdas = 50, plotting=FALSE, add_lambda = FALSE){
  
  dat <- data.frame(X, subjects, y)
  lambda_max <- max(abs(t(X) %*% y))
  lambdas <- linspace(lambda_max,lambda_max*lambda.max.min.ratio, n=n_lambdas)
  BIC_vec<-rep(Inf,length(lambdas))

  for(j in 1:length(lambdas))
  {
    
    suppressWarnings(glm1 <- try(glmmLasso(fix = fixed_form, rnd = rand_form, data=dat,
                      , lambda=lambdas[j],switch.NR=FALSE,final.re=FALSE), silent=TRUE))
    
    if(!inherits(glm1, "try-error"))
    {  
      BIC_vec[j]<-glm1$bic
    }else{break}
    
  }

  if(plotting) plot(lambdas,BIC_vec)

  opt<-which.min(BIC_vec)

  if(add_lambda) chosen_lambdas_vec <<- append(chosen_lambdas_vec,lambdas[opt])

  print(lambdas[opt])

  glm1_final <- glmmLasso(fix = fixed_form, rnd = list(subjects =~ 1), data=dat,
                          , lambda=lambdas[opt],switch.NR=FALSE,final.re=FALSE)

  return(list(names = names(glm1_final$coefficients)[glm1_final$coefficients!=0], vec = glm1_final$coefficients!=0, lambda = lambdas[opt]))
  
}

## before using selective inference two more functions have to be defined:
## * selFun to be either one of the two created here, based on your needs
## * checkFun to check that the model selected is the same as the original
