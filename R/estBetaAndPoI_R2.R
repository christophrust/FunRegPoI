estBetaAndPoI_R2 <-
function(Y, X_mat, add.vars, N, p, potPoI, searchMethod, rho_rng, A_m, X_B, grd, 
         maxPoI = 10, nbest = nbest, intercept = intercept, plotting = plotting,
         scaleSearchPoI = TRUE){
    
    ## Calculate X_i(potPoI_j)
    PoIXValues <- if (scaleSearchPoI) {
                      scale(t(X_mat))[ , potPoI , drop = FALSE]
                  } else {
                      t(X_mat)[ , potPoI , drop = FALSE]
                  }
    
    ## Calculate beta(t) with all potential Points of Impcat
    allEstimates <- estBetaCraKneSa(Y = Y, X_mat = X_mat, add.vars = add.vars, A_m = A_m, X_B=X_B, rho_rng = rho_rng , PoI = potPoI)

    ## Calculate standardized dependent variable without variability due to X_mat %*% estBeta
    Y_woBeta_st  <- scale( Y - (t(X_mat) %*%  allEstimates$estBeta) / p)
    
    ## Choose from Y- int X(t) beta(t) the final PoIChoice
    PoIChoice <- if (searchMethod == "dirSearch") {
                     selectPoI_dirSearch(Y=Y_woBeta_st, X_mat=NULL, PoIXValues = PoIXValues, 
                                         potPoI = potPoI, k = 0, K = 0,
                                         maxPoI = maxPoI, nbest = 1, intercept = intercept, plotting = plotting)
                 } else {
                     selectPoI_fullSearch(Y=Y_woBeta_st, X_mat=NULL, PoIXValues = PoIXValues, 
                                          potPoI = potPoI, k = 0, K = 0,
                                          maxPoI = maxPoI, nbest = 1, intercept = intercept, plotting = plotting)
                 }
    S_choice  <- length(PoIChoice) 

    ## Estimates of S_choice PoI Beta_j, j = 1, ..., S and beta(t)
    if( S_choice != 0 ){ # If there are PoI
        estPoI  <- allEstimates$estAlpha[p+which(potPoI %in% PoIChoice)]

        ## in case of additional variables:
        if (!is.null(add.vars)){
            add.var.coef <- allEstimates$estAlpha[ (p+S_choice+1):(p+S_choice+ncol(add.vars)) ]
            Y_hat        <- (t(X_mat)[ ,1:p] %*% allEstimates$estBeta)/p + t(X_mat[ PoIChoice , , drop=FALSE]) %*% estPoI + add.vars %*% add.var.coef
        } else {
            add.var.coef <- NULL
            Y_hat        <- (t(X_mat)[ ,1:p] %*% allEstimates$estBeta)/p + t(X_mat[ PoIChoice , , drop=FALSE]) %*% estPoI 
        }
        
    } else { # If there aren't PoI

        estPoI <- NULL
        
        ## in case of additional variables:
        if (!is.null(add.vars)){
            add.var.coef <- allEstimates$estAlpha[ (p+S_choice+1):(p+S_choice+ncol(add.vars)) ]
            Y_hat        <- (t(X_mat)[ ,1:p] %*% allEstimates$estBeta)/p + add.vars %*% add.var.coef
        } else {
            add.var.coef <- NULL
            Y_hat        <- (t(X_mat)[ ,1:p] %*% allEstimates$estBeta)/p
        }    
    }

    ## Calc edf's and BICs
    edfsAndBIC <- calPoIBIC( X = t(X_mat), XtX_1Xt = allEstimates$XtX1Xt, Y, Y_hat,
                            rho = allEstimates$rho, S = S_choice + length(add.var.coef), p=p) 
    
    list(estPoI     = estPoI, 
         estTau     = ind2grd(PoIChoice, grd), # The Value in the domain
         estTauGrd  = PoIChoice,               # The Index of the PoI 
         estBeta    = allEstimates$estBeta,
         betaAddVar = add.var.coef,
         gcv        = allEstimates[["gcv"]],
         edfsAndBIC = edfsAndBIC,
         XtX_1Xt    = allEstimates$XtX1Xt[c(1:p, p+which(potPoI %in% PoIChoice)),],
         rho        = allEstimates$rho, 
         selS       = length(PoIChoice), 
         Y_hat      = Y_hat, 
         BIC        = edfsAndBIC$BIC, 
         BIC_sqHat  = edfsAndBIC$BIC_sqHat, 
         eff_df     = edfsAndBIC$eff_df)
}
