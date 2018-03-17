KnePosSarEstimation_fullSearch <-
function(Y, X_mat, add.vars, k, grd, dom, threshold_cpv = 0.95, plotting = FALSE, maxK = 40,maxPoI=10 , intercept = FALSE , nbest = 1){
    ## ########################
    ## A Function doing the KnePosSar algorithm of doing the PoI esimation by directed BIC search
    ## Input:
    ##  - Y: scalar dependent variables (N)
    ##	- X_mat: data matrix (p x N)
    ##  - k: the number of grid indices to to use for potential PoISearch; k <- 15
    ##  - grd: the grid of length p where the functions are observed
    ##  - dom: the domain, where the function is observed
    ##  - plotting=FALSE: should be plots provided?
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
    X_st   <- scale(X_mat) # k<- 30; plotting = FALSE
    potPoI <- searchPotPoi(k, X_st, Y_st, dom=dom, plotting=plotting) 
	
	## Take minimum of up to maxPoI or potPoi number PoI into your choice
    maxS   <- min(length(potPoI$ind.x), maxPoI)
    if(maxS > 0){
        potPoI <- potPoI$ind.x[1:maxS]
    } else{
        potPoI <- NULL
    }
	
    ## Centering of Y and X data
    X_c <- apply( X_mat, 1 , FUN = function(row) return(row - mean(row)) )
    Y_c <- Y - mean(Y)
    
    ## Eigenfunctions, Basis Coefficients and Eigenvalues via PCA
    efSelection <- selectEFwithCumPropVar(X_mat = X_mat, threshold_cpv = threshold_cpv)
    scores  <- X_c %*% efSelection$evecs * 1/p  
    
    Y_scores_df        <- data.frame(cbind(Y_c, scores))
    names(Y_scores_df) <- c("Y", colnames(scores))
    
    maxS               <- min(maxPoI, maxS)
    
    ## If CPV is supplied, estimate the number of k via CPV, then adjust Y and choose the PoI via fulsearch
    if( is.numeric(threshold_cpv) ) { 
        K_choice <- efSelection$K
        
        Y_wo_scores_formula <- as.formula(paste0( "Y ~ -1 +" , paste( colnames(scores)[1:K_choice], collapse=" + ")))
        Y_wo_scores         <- residuals(lm(Y_wo_scores_formula, data = Y_scores_df ))
        
        ## Regress Residuals against potPoI; mean of residuals = 0
        estDataFrame       <- createEstDFforDirectSearch(Y = Y_wo_scores, X = NULL, PoI_vals = X_c[ ,potPoI, drop = FALSE])
        regss              <- regsubsets(x = as.matrix(estDataFrame[ , -1]), y = estDataFrame[ , "Y"],
                                         intercept = intercept,
                                         ## ###################################################
                                         nbest     = nbest,     # return only the best selections with 1, 2, ... , PoI_max predictors. 
                                         really.big = TRUE, 
                                         nvmax = maxPoI)
        
        
        ## Use the BIC in order to extract the overall best model 
        ## from the PoI_max different best models with 1, 2, ... , PoI_max predictors. 
        chosen_X_bool <- summary(regss)$which[which.min(summary(regss)$bic), ]  
        PoIChoice     <- potPoI[ chosen_X_bool==TRUE ]
        S_choice      <- length(PoIChoice)
        BICChoice     <- min(summary(regss)$bic)
        
        ## If no cpv is supplied, estimate the number of K by dirSearch	and the number S_choice by fullSearch,
	## i.e. for each directed k, adjust Y and then choose the poi via full search
    } else if( isTRUE(!threshold_cpv) ) {  
        
        ## Estimate for each number of ks the best amount of PoIs, i.e. direct search over k, full search over PoI
        ## For each k = 0, ..., maxK
        ## one has all the minimum BICs over all possible subsets of S
        BICs <- lapply(0:maxK, function(k) { # k = 1
            ## ###################################################
            ## Regress y against k scores to get residuals without k = 1
            Y_wo_scores_formula <- as.formula(paste0( "Y ~ -1 +" , paste( colnames(scores)[1:k], collapse=" + ")))
            Y_wo_scores <- residuals(lm(Y_wo_scores_formula, data = Y_scores_df	))
            
            ## Regress Residuals against potPoI; mean of residuals = 0
            estDataFrame <- createEstDFforDirectSearch(Y = Y_wo_scores, X = NULL, PoI_vals = X_c[ ,potPoI, drop = FALSE])
            regss        <- regsubsets(x = as.matrix(estDataFrame[ , -1]), y = estDataFrame[ , "Y"],
                                       intercept = intercept,
                                       ## ###################################################
                                       nbest     = nbest,     # return only the best selections with 1, 2, ... , PoI_max predictors. 
                                       really.big = TRUE, 
                                       nvmax     = maxPoI)
            return(list(minBICperK = min(summary(regss)$bic), regList = summary(regss)$which[which.min(summary(regss)$bic), ], sum = summary(lm(Y_wo_scores_formula, data = Y_scores_df))))
        })
        
        ## Take the minium of BIC off all k=0, ..., maxK given the minimum over S to choose K_choice, get afterwards S_choice and therefore PoIChoice
        BICChoice   <- min(sapply(BICs, function(x) x$minBICperK))
        K_choice    <- (0:maxK)[which.min(sapply(BICs, function(x) x$minBICperK))]
        S_choice    <- sum(BICs[[K_choice+1]]$regList) # S_choice = 0
        if(S_choice > 0){
            PoIChoice <- unlist(potPoI[ BICs[[K_choice+1]]$regList ])
        } else {
            PoIChoice <- NULL
        }
        
    } else {
        stop("Problem using CPV")
    }
    
    ## get the final model estimates, i.e. generate data.frame for estimation 
    if( K_choice > 0 ){ # if eigenfunctions are chosen
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ , 1:K_choice, drop = FALSE], PoI_vals = X_c[ , PoIChoice, drop = FALSE])
    } else { # if BIC prefers not using eigenfunctions
        estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = NULL, PoI_vals = X_c[ , PoIChoice, drop = FALSE])
    }	
    
    ## estimate Y againgst eigenfunctions and PoI choices
    lm.fit     <- calcPoIandScoreLM(estDataFrame, s = S_choice, k = K_choice, K = K_choice, S = S_choice)$lmFit
    coefnames  <- names(lm.fit$coeffic)
    score_names  <- coefnames[substr(coefnames , 1,3)!="PoI"]
    poi_names  <- coefnames[substr(coefnames , 1,3)=="PoI"]
    estEF      <- efSelection$evecs[ , score_names]
    estPsi     <- lm.fit$coefficients[score_names]
    estPoI     <- lm.fit$coefficients[poi_names]
	
	
    ## Calculate beta_hat(t), PoI_hat, Y_hat
    beta_hat  <- calcScoreBetaFun(ef = estEF, coefs = estPsi)
    
    if (S_choice > 0) { # if model uses PoIs
        ## Y_i = integral X_i(t) beta_hat(t) dt + sum_k beta_hat_k * X_i(tau_hat_k)
        Y_hat  <- (X_c %*% beta_hat(grd))/ p + X_c[ ,PoIChoice, drop=FALSE] %*% estPoI  
    } else {
        PoIChoice <- NULL
        Y_hat     <- (X_c %*% beta_hat(grd))/ p
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
