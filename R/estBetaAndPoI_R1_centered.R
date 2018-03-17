estBetaAndPoI_R1_centered <-
function(Y, X_mat, add.vars, N, p, potPoI, searchMethod, rho_rng, A_m, X_B, grd, 
         maxPoI = 10, nbest = nbest, intercept = intercept, plotting = plotting){
    ## #####################################################################################
    ## Function performing the R1-estimate 
    
    ## Calculate centered Y and centered X
    Y_c   <- scale(Y, scale = FALSE)
    X_c   <- t(scale(t(X_mat), scale = FALSE)) # matplot(x = grd, y = X_c, type = "l"); dim(X_c); dim(X_mat)
    
    ## Calculate X_i(potPoI_j)
    PoIXValues <- scale(t(X_mat))[ , potPoI , drop = F]
    
    ## Choose from potPoI the PoIChoice searchMethod <- "dirSearch" searchMethod <- "fullSearch"
    PoIChoice <- if (searchMethod == "dirSearch"){
                     selectPoI_dirSearch(Y=scale(Y), X_mat=NULL, PoIXValues = PoIXValues, 
                                         potPoI = potPoI, k = 0, K = 0,
                                         maxPoI = maxPoI, nbest = 1, intercept = intercept, plotting = plotting)
                 } else {
                     selectPoI_fullSearch(Y=scale(Y), X_mat=NULL, PoIXValues = PoIXValues, 
                                          potPoI = potPoI, k = 0, K = 0,
                                          maxPoI = maxPoI, nbest = 1, intercept = intercept, plotting = plotting)
                 }
    S_choice  <- length(PoIChoice) 
    
    ## Estimate Beta(t) and Beta_i for PoI_i
    allEstimates <- estBetaCraKneSa(Y = Y_c, X_mat = X_c, add.vars = add.vars, A_m=A_m, X_B=X_B, rho_rng = rho_rng, PoI = PoIChoice)
    
    ## Estimates of S_choice PoI Beta_j, j = 1, ..., S and beta(t)
    if( S_choice != 0 ){ # If there are PoI
        estPoI  <- allEstimates$estAlpha[(p+1):(p+S_choice)]

        ## in case of additional variables:
        if (!is.null(add.vars)){
            add.var.coef <- allEstimates$estAlpha[ (p+S_choice+1):(p+S_choice+ncol(add.vars)) ]
            Y_hat_c <- (t(X_c)[ ,1:p] %*% allEstimates$estBeta)/p + t(X_c[ PoIChoice , , drop=F]) %*% estPoI + add.vars %*% add.var.coef
        } else {
            add.var.coef <- NULL
            Y_hat_c <- (t(X_c)[ ,1:p] %*% allEstimates$estBeta)/p + t(X_c[ PoIChoice , , drop=F]) %*% estPoI 
        }
        
    } else { # If there aren't PoI

        estPoI <- NULL
        
        ## in case of additional variables:
        if (!is.null(add.vars)){
            add.var.coef <- allEstimates$estAlpha[ (p+S_choice+1):(p+S_choice+ncol(add.vars)) ]
            Y_hat_c <- (t(X_c)[ ,1:p] %*% allEstimates$estBeta)/p + add.vars %*% add.var.coef
        } else {
            add.var.coef <- NULL
            Y_hat_c <- (t(X_c)[ ,1:p] %*% allEstimates$estBeta)/p
        }
        
        
    }

    ## Calc edf's and BICs
    edfsAndBIC      <- calPoIBIC( X = t(X_c), XtX_1Xt = allEstimates$XtX1Xt, Y = Y_c , Y_hat = Y_hat_c,
                                 rho = allEstimates$rho, S = S_choice + length(add.var.coef) , p=p) 

    
    list( estPoI    = estPoI, 
         estTau     = ind2grd(PoIChoice, grd),         # The Value in the domain
         estTauGrd  = PoIChoice,                       # The Index of the PoI 
         estBeta    = allEstimates$estBeta,
         betaAddVar = add.var.coef,
         edfsAndBIC = edfsAndBIC,
         XtX_1Xt    = allEstimates$XtX1Xt,
         rho        = allEstimates$rho, 
         selS       = length(PoIChoice), 
         Y_hat      = Y_hat_c, 
         gcv        = allEstimates[["gcv"]],
         BIC        = edfsAndBIC$BIC, 
         BIC_sqHat  = edfsAndBIC$BIC_sqHat, 
         eff_df     = edfsAndBIC$eff_df)
    
}
