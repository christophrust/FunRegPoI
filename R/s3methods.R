#' @export
plot.FunRegPoI <-function(x,...){
    coefs <- x$coefficients
    estBeta  <-range(beta <- coefs$betaCurve) + diff(range(beta)) * c(-0.2,+0.2)
    grd <- c(0,1)
    taus <- coefs$tauGrd
    ylims <- range(beta) + diff(range(beta)) * c(-0.2,+0.2)
    
    ## par(oma = c(0,0,0,0), mar=c(6,2,3,2))
    plot(y=estBeta , x=grd , col="white" ,xaxt="n", ...)
    ## Abline grey PoI lines
    if (length(taus)>0){
        for (tau in taus) abline( v=tau , col="darkgrey" , lwd=10)
        ## Dates and Tau_i
        for (i in 1:length(taus)) text( x=taus[i] , y=par("usr")[3] , 
                                       labels = bquote( widehat(tau)[.(i) ]==.(round(taus[i], digits=2)))  ,
                                       xpd = TRUE , pos=1, ...)
        ## Beta_i 
        for (i in 1:length(taus)) text( x=taus[i] , y=par("usr")[4] , 
                                       labels = bquote( widehat(beta)[.(i) ]==.(format(coefs$betaPoI[i], digits=2)))  ,
                                       xpd = TRUE , pos=3, ...)
    }
    ## Abline black beta curve
    lines( y=beta, x=seq(0,1, len = length(beta)),  col = "black", lwd =5)
}

#' @export
print.FunRegPoI <-function(x , ...){
    coefs <- x$coefficients
    estBeta  <-coefs$betaCurve
    grd <- seq(0,1, len = length(estBeta))
    ylim <- range(estBeta) + diff(range(estBeta)) * c(-0.2,+0.2)
    taus <- coefs$tauGrd
    bic <- x$model$bic
    if ( x$call != "CKS" ){
        cat("Number of selected Points of Impact: ")
        cat(length(taus))
        cat("\nMaximum PoIs allowed: ")
        cat(x$model$maxPoI)
        if (length(taus)>0){ 
            cat("\nPoint of Impact location(s) (index) in [0,1]: ")
            cat(paste0(round(taus , digits = 2) , " (" , coefs$tauInd , ")" , collapse = " "))
            cat("\nPoint of Impact influences: ")
            cat(format(coefs$betaPoI , digits = 3))
        }
        if (lenCoef <- length(coefs$betaAddVar) > 0 ){
            cat("\n Influence of additional variables:")
            for (i in 1:lenCoef) {
                cat(paste0("\n", names(coefs$betaAddVar)[i] , ": " , round(coefs$betaAddVar, digits=3)))
            }
        }
        cat(paste0("\nBIC: " , round(bic,digits = 3)))
        cat("\nOptimal k (delta): ")
        cat( paste0(x$model$k , " (" , round(x$model$delta, digits = 2) , ")" , collapse = " "))
        cat("\n")
    }
    else {
        cat("CKS model (without PoIs)\n")
        cat("Selected rho: ")
        cat(format(x$model$rho , digits = 5))
        cat("\nPossible range for rho: ")
        cat(paste0("[",x$model$rho_space[1] , "," , x$model$rho_space[2], "]"))
        cat("\nValue of GCV: ")
        cat(format(x$model$gcv , digits = 3))
        cat("\n")
    }
}

#' @export
summary.FunRegPoI <- function(object, confidence.level = 0.95,...){
    coefs <- object$coefficients
    estBeta <- coefs$betaCurve
    grd <- seq(0,1, len = length(estBeta))
    ylim <- range(estBeta) + diff(range(estBeta)) * c(-0.2,+0.2)
    taus <- coefs$tauGrd
    tausInd <- object$coef$tauInd
    BIC <- object$model$bic
    eff.df <- object$model$eff_df
    call <- object$call
    betaPoI <- object$betaPoI
    betaAddVar <- object$betaAddVar
    selS <- object$selS

    obj <- list(estBeta = estBeta,
                taus = taus,
                tausInd = tausInd,
                BIC = BIC,
                eff.df = eff.df,
                selS = selS,
                call = call,
                model = object$model)
                
    if (call != "KPS"){

        ## point wise standard errors for beta-curve and standard errors for PoIs and add.vars:
        XtX1_Xt <- object$model$XtX_1Xt
        s2.resid <- sum(object$residuals^2)/( length(object$residuals) - eff.df )
        coef.serr.ptw  <- sqrt(diag(tcrossprod(XtX1_Xt)) * s2.resid )
        
        if ( (length(coefs$betaAddVar) + length(coefs$betaPoI))>0){

            ## PoI coef and AddVar coef summaries
            coef <- matrix(NA , ncol = 4, nrow = length(coefs$betaPoI) + length(coefs$betaAddVar))
            colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            rownames(coef) <- c(if (length(coefs$betaPoI)>0) paste0("PoI.",1:length(coefs$betaPoI)) else NULL ,
                                names(coefs$betaAddVar))
            coef[,1] <- c(coefs$betaPoI , coefs$betaAddVar)
            coef[,2] <- coef.serr.ptw[ ( length(estBeta)+1):length(coef.serr.ptw)]
            coef[,3] <- coef[,1] / coef[,2]
            coef[,4] <- 1 - pnorm( abs( coef[,3]))

            obj <- c(obj , list(coefficients = coef))
        } else {
            obj <- c(obj , list(coefficients = NULL))
        }
                
    

        ## beta curve confidence interval
        q.norm <-  pnorm(1-(1- confidence.level)/2 )
        confint.beta <- cbind(lower = estBeta - q.norm * coef.serr.ptw[1:length(estBeta)],
                          upper = estBeta + q.norm * coef.serr.ptw[1:length(estBeta)])

        obj <- c(obj , list( confint.beta = confint.beta))
        
    } else {
        obj <- c(obj , list(betaPoI = betaPoI,
                            betaAddVar = betaAddVar))
    }

    class(obj) <- "summary.FunRegPoI"
    obj
}



#' @export
print.summary.FunRegPoI <-function(x , signif.stars = TRUE,digits = 3,...){
    coefs <- x$coefficients
    estBeta  <-x$estBeta
    grd <- seq(0,1, len = length(estBeta))
    ylim <- range(estBeta) + diff(range(estBeta)) * c(-0.2,+0.2)
    taus <- x$taus
    bic <- x$BIC
    if ( x$call != "CKS" ){
        cat("\nNumber of selected Points of Impact: ")
        cat(length(taus))
        cat("\nMaximum PoIs allowed: ")
        cat(x$model$maxPoI)
        if (length(taus)>0 || length(x$betaAddVar) > 0){ 
            cat("\n\nPoint of Impact location(s) (index) in [0,1]: ")
            cat(paste0(round(taus , digits = 2) , " (" , x$tausInd , ")" , collapse = " "))
            if (x$call != "KPS"){
                cat(paste0("\n\nPoint of Impact influences",
                           if ( length(x$betaAddVar) > 0 ) " and additional Variables:\n\n" else ":\n\n"))
                printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                             na.print = "NA", ...)
            } else {
                cat("\nPoint of Impact influences: ")
                cat(format(x$betaPoI , digits = 3))
                
                if (lenCoef <- length(x$betaAddVar) > 0 ){
                    cat("\n Influence of additional variables:")
                    for (i in 1:lenCoef) {
                        cat(paste0("\n", names(x$betaAddVar)[i] , ": " , round(x$betaAddVar, digits=3)))
                    }
                }
            }
        }
        
        cat(paste0("\nBIC: " , round(bic,digits = 3)))
        cat("\nOptimal k (delta): ")
        cat( paste0(x$model$k , " (" , round(x$model$delta, digits = 2) , ")" , collapse = " "))
        cat("\n")
    }
    else {
        cat("CKS model (without PoIs)\n")
        cat("Selected rho: ")
        cat(format(x$model$rho , digits = 5))
        cat("\nPossible range for rho: ")
        cat(paste0("[",x$model$rho_space[1] , "," , x$model$rho_space[2], "]"))
        cat("\nValue of GCV: ")
        cat(format(x$model$gcv , digits = 3))
        cat("\n")
    }
}


#' @export
plot.summary.FunRegPoI <-function(x,...){
    
    if (x$call != "KPS"){
        beta.rng  <-range(x$confint.beta) + diff(range(x$confint.beta)) * c(-0.2,+0.2)
        estBeta <- x$estBeta
    } else {
        beta.rng  <-range(estBeta <- x$estBeta) + diff(range(x$estBeta)) * c(-0.2,+0.2)
    }
    taus <- x$taus
    if (length(taus)>0) {
        PoIcoefs <- x$coefficients[1:length(taus),1]
    }
    ylims <- range(estBeta) + diff(range(estBeta)) * c(-0.2,+0.2)
    
    ## par(oma = c(0,0,0,0), mar=c(6,2,3,2))
    plot(y=beta.rng , x=c(0,1) , col="white" ,xaxt="n", xlab= "grd", ylab ="beta",...)

    ## plot the pointwise confidence interval:
    if (x$call != "KPS"){
        grd <- seq(0,1, len = length(estBeta))
        y.ci <- c( x$confint.beta[,"lower"], x$confint.beta[length(estBeta):1,"upper"])
        x.ci <- c(grd, grd[length(estBeta):1])
        polygon(x = x.ci, y = y.ci, col = "lightgrey", border = NA)
    }
    
    ## Abline grey PoI lines
    if (length(taus)>0){
        for (tau in taus) abline( v=tau , col="darkgrey" , lwd=10)
        ## Dates and Tau_i
        for (i in 1:length(taus)) text( x=taus[i] , y=par("usr")[3] , 
                                       labels = bquote( widehat(tau)[.(i) ]==.(round(taus[i], digits=2)))  ,
                                       xpd = TRUE , pos=1, ...)
        ## Beta_i 
        for (i in 1:length(taus)) text( x=taus[i] , y=par("usr")[4] , 
                                       labels = bquote( widehat(beta)[.(i) ]==.(format(PoIcoefs[i], digits=2)))  ,
                                       xpd = TRUE , pos=3, ...)
    }
    
        
    ## Abline black beta curve
    lines( y=estBeta, x=seq(0,1, len = length(estBeta)),  col = "black", lwd =5)
}



## predict currently without add.var support for out of sample prediction
#' @export
predict.FunRegPoI <- function(object, newx = NULL){

    ## insample
    if (is.null(newx)) {
        return(object$fitted)
    }

    ## new data
    if (object$call != "CKS"){
        
        betaCurve <- object$coefficients$betaCurve
        betaPoI <- object$coefficients$betaPoI
        tauInd <- object$coefficients$tauInd
        if (object$model$center) {
            data <- object$data
            Xmean <- rowMeans(data$X_mat)
            Ymean <- mean(data$Y)
            
            intercept <- Ymean - sum(Xmean * betaCurve / length(betaCurve) ) - sum(Xmean[tauInd] * betaPoI)
        } else {
            intercept <- 0
        }

        y.pred <- as.vector( t(newx) %*% betaCurve / length(betaCurve) + t(newx[tauInd,]) %*% betaPoI + intercept)
        
    } else {
        betaCurve <- object$coefficients$betaCurve
        if (object$model$center) {
            data <- object$data
            Xmean <- rowMeans(data$X_mat)
            Ymean <- mean(data$Y)
            
            intercept <- Ymean - sum(Xmean * betaCurve / length(betaCurve) )
        } else {
            intercept <- 0
        }
        y.pred <- as.vector( t(newx) %*% betaCurve / length(betaCurve) + intercept)
    }
    ## return
    y.pred
}
