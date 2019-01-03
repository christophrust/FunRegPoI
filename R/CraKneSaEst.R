CraKneSaEst <-
function(Y, X_mat, add.vars, A_m, X_B, rho_rng = c(1e-6,2e2) , rho = NULL, center = TRUE , ...){
    ## wrapper around estBetaCraKneSa

    
    if (center){
        ## Calculate centered Y and centered X
        Y_c   <- scale(Y, scale = FALSE)
        X_c   <- t(scale(t(X_mat), scale = FALSE)) # matplot(x = grd, y = X_c, type = "l"); dim(X_c); dim(X_mat)
        
        betaEstimates  <- estBetaCraKneSa(Y = Y_c, X_mat = X_c, add.vars = add.vars, A_m = A_m, X_B = X_B, rho_rng = rho_rng, rho = rho, PoI = NULL , ...)
        
    } else{
        betaEstimates   <- estBetaCraKneSa(Y, X_mat, add.vars = add.vars, A_m = A_m, X_B = X_B, PoI = NULL, rho_rng = rho_rng,rho = rho, ...)
    }
    

    if (!is.null(add.vars)){
        betaAddVar <- betaEstimates$estAlpha[ (p+1):(p+ncol(add.vars)) ]
        names(betaAddVar) <-
            if (!is.null(colnames(add.vars))) {
                colnames(add.vars)
            } else {
                paste0("add.var", 1:ncol(add.vars))
            }
    } else betaAddVar <- NULL

    ## fitted and residual
    Y_hat <- if (center){
                 as.vector( (t(X_c) %*% betaEstimates$estBeta)/length(betaEstimates$estBeta) ) +
                     if (!is.null(add.vars)) add.vars %*% betaAddVar else 0
             } else{
                 as.vector( (t(X_mat) %*% betaEstimates$estBeta)/length(betaEstimates$estBeta) ) +
                     if (!is.null(add.vars)) add.vars %*% betaAddVar else 0
             }
    
    
    edfsAndBIC <- calPoIBIC( X = t(X_mat), XtX_1Xt = betaEstimates$XtX1Xt, Y, Y_hat,
                            rho = betaEstimates$rho, S = length(betaEstimates$betaAddVar), p=length(betaEstimates$estBeta)) 

    
    estObj <- list(
        coefficients = list(
            betaCurve = as.vector(betaEstimates$estBeta) ,
            betaPoI = NULL ,
            betaAddVar = betaAddVar,
            tauGrd = NULL ,
            tauInd = NULL,
            selS = 0 ),
        residuals =  Y - Y_hat - if(center) mean(Y) else 0,
        fitted =  Y_hat + if(center) mean(Y) else 0,
        data = list(Y = Y, X_mat = X_mat, add.vars = add.vars),
        call = "CKS" ,
        model = list(
            k = NULL ,
            delta = NULL ,
            rho = betaEstimates[["rho"]],
            rho_rng = range(rho_rng) ,
            gcv = betaEstimates[["gcv"]] ,
            cor = NULL,
            eff_df = edfsAndBIC$eff_df,
            bic = edfsAndBIC$BIC,
            XtX_1Xt = betaEstimates$XtX1Xt,
            bic_sq= edfsAndBIC$BIC_sqHat,
            hat_matrix= NULL,
            center = center),
        kSeqRes = NULL
    )
    return(estObj)
}
