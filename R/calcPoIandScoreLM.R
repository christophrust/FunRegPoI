calcPoIandScoreLM <-
function(estDataFrame, s, k, K, S){
    ## ########################
    ## A Function calculating for arbitrary s PoIs the corresponding lm models and their BICs
    ## Input: # k =1; s = 1;K = maxK; S = maxS
    ##  - s: the number of included PoIs
    ##	- k: the number of included eigenfunctions
    ##  - estDataFrame: a named data.frame with 1 + K + S columns of Y, the K scores and S potential PoIs
    ##  - K: the maximum number of Eigenfunctions
    ##  - S: the maximum number of PoIs
    ##  - model_out: should the lm fit be returned?
    ## Output:
    ##  - The MLR BIC
    dep_name  <- names( estDataFrame )[1]
    
    ## print(paste0("s = ", s, " k = ", k , " K = ", K, " S = ", S))
    if(s > 0 && s <= S && k <= K){
        ## Get the dependent variable s name
        ef_names   <- NULL 
        poi_names  <- NULL 
        if ( k > 0 ){ ## If other regressors are in estDataFrame, get the name of them
            ef_names <- (names( estDataFrame )[-1])[1:k]	
        } 
        if ( s > 0 ){ ## Get PoI names
            poi_names <- (names( estDataFrame )[-(1:(K+1))])[1:s]
        } 
        ## Generate a formula depending on the number of chosen PoI and estimate it via LM
        s_formula    <- as.formula(paste0( dep_name, " ~ -1 +" , paste( c(ef_names, poi_names) , collapse=" + ")))
        lmfit        <- lm(formula = s_formula , data = estDataFrame)
        BIC          <- calcBIClm(lmfit)[2]
        #GCV          <- calcGCVlm(lmfit)[2]
        return(list(BIC = BIC, #GCV = GCV,
                    lmFit = lmfit))
		
    } else if ( s == 0 ) {
        ef_names    <- NULL 
        if ( k > 0 ){ ## If other regressors are in estDataFrame, get the name of them
            ef_names  <- (names( estDataFrame )[-1])[1:k]
            s_formula <- as.formula(paste0( dep_name, " ~ -1+" , paste(ef_names , collapse = " + ")))
        } 
        else {
            s_formula <- as.formula(paste0( dep_name, " ~ -1" ))
        }
        lmfit  <- lm(formula = s_formula , data = estDataFrame)

        BIC    <- calcBIClm(lmfit)[2]
        #GCV    <- calcGCVlm(lmfit)[2]
        return(list(BIC = BIC, #GCV = GCV,
                    lmFit = lmfit))
        
    } else{
        stop("Something is wrong with estDataFrame with s and k")
    }
}
