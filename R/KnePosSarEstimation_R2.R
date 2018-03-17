KnePosSarEstimation_R2 <-
function(Y, X_mat, add.vars, N, potPoI, searchMethod, grd, dom, plotting = FALSE,
         maxK = 40, maxPoI=10, efSelection = NULL, intercept = FALSE){
    ## #####################################################################
    ## Function performing the R2-estimate using the Eigenfunction-Expansion
    ## given delta (decorrelation width)

    ## todo (additional variables currently only implemented in CraKneSaPoI):
    add.var.coef <- NULL
    
    ## Calculate centered Y and centered X
    p     <- nrow(X_mat)
    Y_c   <- scale(Y, scale = FALSE)
    X_c   <- scale(t(X_mat), scale = FALSE) # matplot(x = grd, y = X_c, type = "l"); dim(X_c); dim(X_mat)

    ## Calculate X_i(potPoI_j)
    PoIXValues <- scale(t(X_mat))[ , potPoI , drop = F]


    ## Eigenfunctions, Basis Coefficients and Eigenvalues via PCA
    if (is.null(efSelection)){
        efSelection  <- selectEFwithCumPropVar(X_mat = X_mat)
    }
    scores       <- X_c %*% efSelection$evecs * 1/p

    
    ## Estimate Beta(t) and Beta_i for PoI_i
    ## select number of eigenfunctions wrt to gcv
    estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ ,1:maxK], PoI_vals = X_c[ ,potPoI, drop=FALSE])

    ## compute gcv's
    gcvs <- vapply(0:maxK , function(k){
        if (k==0){
            if (is.null(potPoI)){
                y <- estDataFrame[,1]
                n <- length(y)
                RSS <- sum(y^2)
                return(1/n * RSS)  # here: edf = 0
            } else {
                lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, c((maxK+2):ncol(estDataFrame))]))
            }
        } else {
            if (is.null(potPoI)) {
                lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, 2:(k+1) ]))
            } else {
                lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, c(2:(k+1),(maxK+2):ncol(estDataFrame))]))
            }
        }
        calcGCVlm(lmEst)[2]
    }, 0)
    
    K_choice <- which.min(gcvs)-1
    
    ## the corresponding estimate:
    idx <- c( if (K_choice==0) NULL else 2:(K_choice+1), if (is.null(potPoI)) NULL else (maxK+2):ncol(estDataFrame))
    if (length(idx)>0){
        lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, idx ] ))
    } else {
        lmEst <- list(residuals = estDataFrame[,1], coefficients = rep(0, times = 2))
    }

    
    ## 'clean' Y for effect of beta(t)
    Y_woBeta_st  <- scale( Y - scores[,1:K_choice, drop = FALSE] %*%  lmEst$coefficient[ 1:K_choice ])

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
    
    ## get the final model estimates names(estDataFrame)
    if(K_choice != 0){
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ , 1:K_choice, drop=FALSE], PoI_vals = X_c[ , PoIChoice, drop=FALSE])
    } else {
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = NULL, PoI_vals = X_c[ , PoIChoice, drop=FALSE])
    }
    
    ## estimate Y againgst eigenfunctions and PoI choices
    LMandBIC    <- calcPoIandScoreLM(estDataFrame, s = S_choice, k = K_choice, K = K_choice, S = maxPoI)
    lm.fit      <- LMandBIC$lmFit
    coefnames   <- names(lm.fit$coeffic)
    score_names <- coefnames[substr(coefnames , 1,3)!="PoI"]
    poi_names   <- coefnames[substr(coefnames , 1,3)=="PoI"]
    estEF       <- efSelection$evecs[ , score_names]
    estPsi      <- lm.fit$coefficients[score_names]
    estPoI      <- lm.fit$coefficients[poi_names]
    
    ## Calculate beta_hat(t), PoI_hat, Y_hat
    beta_hat   <- calcScoreBetaFun(ef = estEF, coefs = estPsi, p=p)
    if (S_choice > 0) {
        ## Y_i = integral X_i(t) beta_hat(t) dt + sum_k beta_hat_k * X_i(tau_hat_k)
        Y_hat  <- (X_c %*% beta_hat(grd))/ p + X_c[ ,PoIChoice , drop=FALSE] %*% estPoI 
        
    } else {
        Y_hat  <- (X_c %*% beta_hat(grd))/ p
    }
    hat_matrix <- NULL
    list(selS = S_choice, 
         K_choice = K_choice,
         estPoI = estPoI,
         betaAddVar = add.var.coef,
         estTau = ind2grd(PoIChoice , grd),
         estTauGrd = PoIChoice,
         estPsi = estPsi,
         scores = scores ,
         estEF = estEF,
         Y_hat = Y_hat, 
         estBeta = beta_hat(grd),
         BIC = LMandBIC$BIC,
         hat_matrix = hat_matrix
         )
}
