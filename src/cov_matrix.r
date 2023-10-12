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