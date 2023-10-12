## function to return metrics (True Positive Rate, False Discovery Rate)

metrics <- function(selected, true){

    tpr <- sum(selected*true)/sum(true)
    fdr <- sum(selected==1 & true==0)/sum(selected)

    return(list(tpr = tpr, fdr = fdr))
}