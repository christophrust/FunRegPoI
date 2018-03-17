calcGCVlm <-
function(lmObj){
    ## ###########################################
    ## A Function caluclating the GCV of an lm fit
    n   <- length(lmObj$residuals)
    edf <- n - lmObj$df.residual
    RSS <- sum(weighted.residuals(lmObj)^2)
    gcv <- 1/n * RSS/((1-edf/n)^2 ) ## see, for instance eq 7.52 in Elements of Statistical Learning
    ## return
    c(edf=edf, gcv=gcv)
}
