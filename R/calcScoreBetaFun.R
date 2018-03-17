calcScoreBetaFun <-
function(ef, coefs , p){
    ## returns a function object from (evaluated) eigenfunctions and coefficients
    if ( length(coefs) == 1 ){
        y  <- ef * coefs
    } else{
        y  <- ef %*% coefs
    }
    return(splinefun( y =  y, x = seq(0,1 , len=p)))
}
