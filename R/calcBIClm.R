calcBIClm <-
function(fit){
    ## ###########################################
    ## A Function caluclating the BIC of an lm fit
    n   <- length(fit$residuals)
    edf <- n - fit$df.residual
    RSS <- sum(weighted.residuals(fit)^2)
    bic <- n * log(RSS/n) + edf * log(n)
    if (bic==-Inf) {
        stop("BIC is -Inf")
    } else {
        return(c(edf, bic))
    }
}
