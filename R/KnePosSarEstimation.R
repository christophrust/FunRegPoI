KnePosSarEstimation <-
function(Y, X_mat, grd, add.vars, maxPoI = 8 ,threshold_cpv=0.95, estOrder = NULL, searchMethod = "dirSearch" ,
         k_seq = 1:floor(length(grd)/6) , maxK = 40, exPost = FALSE, ...){
    ## Function calculating the eigenfunction-based estimate of FunRegPoI. If estOrder ==NULL, the standard
    ## KPS-Estimation is done, in this case, exPost has no effect
    ##
    ##
    ## tocheck: effect of dom in KnePosSarEstimation, adjust in KnePosSarEstimation_R{1,2}!
    ## centering??
    ##
    ## if estOrder =="R2" and exPost =="FALSE": "PES"-estimation (with EF instead of Splines)
    ## if estOrder =="R2" and exPost =="TRUE": "PESES"-estimation (with EF instead of Splines)
    if (is.null(estOrder)){
        ## wrapper around KnePosSarEstimation_fullSearch and KnePosSarEstimation_dirSearch, already over all k
        dom <- range(grd)
        ## for each k in k_seq calculate the KPS
        KnePosSarcomResults <- lapply(k_seq, function(k) {
            if (searchMethod=="dirSearch") {
                KnePosSarEstimation_dirSearch(Y = Y, X_mat = X_mat, add.vars = add.vars,k = k, grd = grd, dom = dom,
                                              threshold_cpv = threshold_cpv, plotting=FALSE, maxK = maxK )
            } else {
                KnePosSarEstimation_fullSearch(Y = Y, X_mat = X_mat, add.vars = add.vars, k = k, grd = grd, dom = dom,
                                               threshold_cpv = threshold_cpv, plotting=FALSE, maxK = maxK)
            }
        })
        
        bics <- lapply(KnePosSarcomResults , function(x) x$BIC)
        opt  <- which.min(unlist(bics))
        optEntry <- KnePosSarcomResults[[opt]]
        
    } else {
        
        Y_st   <- scale(Y)
        X_st  <- t(scale(t(X_mat)))  # dim(X_st)
        possible_deltas  <- ind2grd(k_seq, grd)

        ## get empirical eigenfunctions:
        efSelection  <- selectEFwithCumPropVar(X_mat = X_mat, threshold_cpv = threshold_cpv)
        
        
        ## for each k in k_seq calculate the KPS
        KnePosSarcomResults <- lapply(k_seq, function(k) {
            dom         <- range(grd)
            PoISearch   <- searchPotPoi(k, X_st, Y_st, dom=dom, plotting = FALSE) # with standardized values; might be better
            potPoI      <- PoISearch[["ind.x"]]
            cor         <- PoISearch[["cor"]]
            maxS        <- length(potPoI)
            
            if (estOrder=="R1"){
                KPSest <- KnePosSarEstimation_R1(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=potPoI,
                                                 searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                                 maxK = 40, maxPoI=10, efSelection = efSelection)
            } else {
                KPSest <- KnePosSarEstimation_R2(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=potPoI,
                                                 searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                                 maxK = 40, maxPoI=10, efSelection = efSelection)
            }
            KPSest[["cor"]]      <- PoISearch[["cor"]]
            KPSest[["k"]]        <- k
            KPSest[["delta"]]    <- possible_deltas[which(k == k_seq)]
            KPSest
        })
        names(KnePosSarcomResults)       <- possible_deltas

        if (exPost){
            exPostEstimation <- lapply(KnePosSarcomResults, function(entry) {
                
                if (estOrder=="R1"){
                    KPSest <- KnePosSarEstimation_R1(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=entry[["estTauGrd"]],
                                                     searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                                     maxK = 40, maxPoI=10, efSelection = efSelection)
                } else {
                    KPSest <- KnePosSarEstimation_R2(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=entry[["estTauGrd"]],
                                                     searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                                     maxK = 40, maxPoI=10, efSelection = efSelection)
                }
                KPSest
            })
        }

        ## nest ExPost results in KnePosSarcomResults
        if (exPost){
            for( delta in possible_deltas){ # delta <- possible_deltas[1]
                delta_c <- which(delta == possible_deltas)
                KnePosSarcomResults[[delta_c]][["exPostEstimation"]] <- exPostEstimation[[delta_c]]
            }
        }
        ## Get Optima via BIC
        optimum   <- which.min(sapply(KnePosSarcomResults, function(x) x$BIC))
        
        if (exPost){
            exPostOptimum  <- which.min(sapply(KnePosSarcomResults, function(x) x[["exPostEstimation"]]$BIC))
        }
        
        ## save the selected entries
        KnePosSarcomResults[["opt"]]        <- list(KnePosSarcomResults[[optimum]], X = X_mat,Y = Y)
        names(KnePosSarcomResults[["opt"]]) <- c("res", "X", "Y")
        
        if (exPost){
            KnePosSarcomResults[["exPostOptimum"]] <-
                KnePosSarcomResults[[exPostOptimum]][["exPostEstimation"]]
        }
        ## get the selected estimation result:
        optEntry <- if (exPost) {
                        KnePosSarcomResults[["exPostOptimum"]]
                    } else {
                        KnePosSarcomResults[["opt"]]$res
                    }
        
    }
        


    if (!is.null(add.vars)){
        names(KnePosSarcomResults[[opt]]$betaAddVar) <-
            if (!is.null(names(add.vars))) {
                colnames(add.vars)
            } else {
                paste0("add.var", 1:col(add.vars))
            }
    }
    
    ## collect and make the return object
    estObj <- list(
        coefficients = list(
            betaCurve = as.vector( optEntry$estBeta) ,
            betaPoI = optEntry$estPoI ,
            betaAddVar = optEntry$betaAddVar,
            tauGrd = optEntry$estTau ,
            tauInd = optEntry$estTauGrd,
            selS = optEntry$selS) ,
        residuals = Y - as.vector(optEntry$Y_hat) - mean(Y) ,
        fitted = as.vector(optEntry$Y_hat) + mean(Y),
        data = list(Y = Y, X_mat = X_mat),
        call = "KPS" ,
        model = list(
            k = optEntry$k ,
            delta = ind2grd(optEntry$k , grd) ,
            cpv = threshold_cpv,
            maxPoI = maxPoI,
            K =  optEntry$K_choice,
            K_space= if(!missing(threshold_cpv)) NULL else c(1,maxK),
            estEigenFct = optEntry$estEF ,
            estScores = optEntry$scores,
            eff_df = optEntry$K_choice + optEntry$selS + length(add.vars) ,
            bic = optEntry$BIC , 
            bic_sq= NULL, ##
            hat_matrix= optEntry$hat_matrix,
            center = TRUE),
        kSeqRes = KnePosSarcomResults
    )
    return(estObj)
}
