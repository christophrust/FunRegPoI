KnePosSarEstimation_dirSearch <-
function(Y, X_mat, add.vars, k, grd, dom, threshold_cpv = 0.95, plotting = FALSE, maxK = 40, maxPoI=10){
    ## ########################
    ## A Function performing the KnePosSar algorithm of doing the PoI esimation by directed BIC search
    ## Input:
    ##  - Y: scalar dependent variables (N)
    ##  - X_mat: data matrix (p x N)
    ##  - k: the number of grid indices to to use for potential PoISearch; k <- 1
    ##  - grd: the grid of length p where the functions are observed
    ##  - dom: the domain, where the function is observed
    ##  - plotting=FALSE: should be plots provided?
	##  - threshold_cpv = FALSE
    ## Output:
    ##  - k: the k which was used
    ##  - selS: the chosen number of PoIs
    ##  - estPoI: the estimated Parameters of chosen PoIs
    ##  - estTau: the estimated taus on the domain
    ##  - estTauGrd: the estimated taus on the grid
    ##  - Y_hat: the predicted values
    ##  - estBeta: beta_hat(t) function
    ##  - BIC: the minimal BIC after which the model was chosen
    
    ## todo (additional variables currently only implemented in CraKneSaPoI):
    add.var.coef <- NULL
    
    ## Get potential points of impact
    p      <- nrow(X_mat)
    Y_st   <- scale(Y)
    X_st   <- scale(X_mat) 
    potPoI <- searchPotPoi(k, X_st, Y_st, dom=dom, plotting=plotting) 
    
    ## Take minimum of up to maxPoI or potPoi number PoI into your choice
    maxS   <- min(length(potPoI$ind.x), maxPoI)
    if(maxS > 0){
        potPoI <- potPoI$ind.x[1:maxS]
    } else{
        potPoI <- NULL
    }
    
    ## Centering of Y and X data
    X_c	<- apply( X_mat, 1 , FUN = function(row) return(row - mean(row)) )
    Y_c	<- Y - mean(Y)
    
    ## Eigenfunctions, Basis Coefficients and Eigenvalues via PCA
    efSelection  <- selectEFwithCumPropVar(X_mat = X_mat, threshold_cpv = threshold_cpv)
    scores       <- X_c %*% efSelection$evecs * 1/p  
    
    ## If CPV is supplied, estimate the number of k via CPV using the selection by the selectEFwithCumPropVar function
    if( is.numeric(threshold_cpv) ) { # 
        ## If CPV is supplied, estimate the number of k via CPV
        K_choice <- efSelection$K
        
        ## Prepare Data for BIC Choice
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ ,1:K_choice], PoI_vals = X_c[ ,potPoI, drop=FALSE])
        
        # A directed S search over all number of S; choose model via minimal BIC
        LMandBICs <- lapply(0:maxS, function(s) calcPoIandScoreLM(estDataFrame = estDataFrame, s = s, k = K_choice, K = K_choice, S = maxS))
        S_choice  <- (0:maxS)[ which.min(sapply(LMandBICs, function(x) x[["BIC"]])) ] # S_choice = 5
        BICChoice <- min(sapply(LMandBICs, function(x) x[["BIC"]]))
        
	## If no cpv is supplied, estimate the number of K by dirSearch	and the number S_choice by fullSearch
    } else if( isTRUE(!threshold_cpv) ) { 
        
        ## Prepare Data for directed BIC Search over all K and all S head(estDataFrame)
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ ,1:maxK], PoI_vals = X_c[ ,potPoI, drop=FALSE])
        
        ## A directed S search over all number of S and all number of K up to maxK
        ## For each K generate a list of all BICs per s in S

        ## BICs <- lapply(0:maxK, function(k) { sapply(0:maxS, function(s) calcPoIandScoreLM(estDataFrame = estDataFrame, s = s, k = k, K = maxK, S = maxS)$BIC)})
        ## way faster:
        BICs <- lapply(0:maxK , function(k){
            vapply( 0:maxS, function(s) {
                if (k==0){
                    if (s==0){
                        y <- estDataFrame[,1]
                        n <- length(y)
                        RSS <- sum(y^2)
                        return(n * log(RSS/n))  # here: edf = 0
                    } else {
                        lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, c((maxK+2):(maxK+s+1))]))
                    }
                } else {
                    if (s==0) {
                        lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, 2:(k+1) ]))
                    } else {
                        lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, c(2:(k+1),(maxK+2):(maxK+s+1))]))
                    }
                }
                calcBIClm(lmEst)[2]
            }, 0)
        })

        ## Generate for each k optimal included 
        BICs_per_k <- cbind(K_ind = 0:maxK, t(sapply(BICs, function(bic_k) return(list(S_ind = which.min(bic_k), K_BIC = min(bic_k))))))
        optSoptK   <- BICs_per_k[  which.min(BICs_per_k[ , "K_BIC"]) , c("K_ind", "S_ind", "K_BIC")]
        BICChoice  <- unlist(optSoptK)["K_BIC"]
        S_choice   <- (0:maxS)[ optSoptK[["S_ind"]] ] 		# If optSoptK["S_ind"]=1 --> (0:maxS)[ (optSoptK[["S_ind"]]) ] = 0 -> adjustment
        K_choice   <- (0:maxK)[ (optSoptK[["K_ind"]]+1) ] 	# K_choice = 3
    } else {
        stop("Problem using CPV")
    }
    
    if(S_choice > 0){
        PoIChoice <- potPoI[1:S_choice]
    } else {
        PoIChoice <- NULL
    }
    
    ## get the final model estimates names(estDataFrame)
    if(K_choice != 0){
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ , 1:K_choice, drop=FALSE], PoI_vals = X_c[ , PoIChoice, drop=FALSE])
    } else {
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = NULL, PoI_vals = X_c[ , PoIChoice, drop=FALSE])
    }
    
    ## estimate Y againgst eigenfunctions and PoI choices
    lm.fit     <- calcPoIandScoreLM(estDataFrame, s = S_choice, k = K_choice, K = K_choice, S = maxS)$lmFit
    coefnames  <- names(lm.fit$coeffic)
    score_names  <- coefnames[substr(coefnames , 1,3)!="PoI"]
    poi_names  <- coefnames[substr(coefnames , 1,3)=="PoI"]
    estEF      <- efSelection$evecs[ , score_names]
    estPsi     <- lm.fit$coefficients[score_names]
    estPoI     <- lm.fit$coefficients[poi_names]
    
    ## Calculate beta_hat(t), PoI_hat, Y_hat
    beta_hat   <- calcScoreBetaFun(ef = estEF, coefs = estPsi, p=p)
    if (S_choice > 0) {
        ## Y_i = integral X_i(t) beta_hat(t) dt + sum_k beta_hat_k * X_i(tau_hat_k)
        Y_hat  <- (X_c %*% beta_hat(grd))/ p + X_c[, PoIChoice ,drop=FALSE] %*% estPoI
        
    } else {
        Y_hat  <- (X_c %*% beta_hat(grd))/ p
    }
    hat_matrix <- NULL
    list(k  = k, 
         selS = S_choice, 
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
         BIC = BICChoice,
         hat_matrix = hat_matrix
         )
}
