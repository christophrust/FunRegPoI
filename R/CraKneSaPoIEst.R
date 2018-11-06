CraKneSaPoIEst <- function(Y , X_mat , grd , add.vars, A_m, X_B, maxPoI=8, dom = range(grd) ,
                           k_seq = 1:floor(length(grd)/6), rho_rng = c(1e-6,2e2),
                           estOrder = "R2" , searchMethod = "dirSearch" , nbest=1 ,
                           intercept =FALSE , opt.crit = "BIC_sqHat", exPost = TRUE,
                           maxStepsES = 3, center = TRUE , ...) {
    
    ## ###########################################
    ## A Function caluclating the PESES estimator
    ## Input:
    ##  - Y: dependent scalar values in R^n
    ##  - X_mat: independent functional values; p x N
    ##  - grd: a grid sequence of length p where functions are observed
    ##     - dom: the corresponding domain [a,b] with the same length
    ##  - k_seq: the k_delta sequence in terms of indices instead of absolute value
    ## one can interchange between k_delta and delta via:
    ## k_seq <- grd2ind(possible_deltas, grd) and possible_deltas <- ind2grd(k_seq, grd)
    ##  - rho_seq: Over which rho can optimized? a domain in R^2_+ is enough 
    ##  - estOrder: "R1" or "R2" described in the Appendix; R2 belongs to PESES
    ##  - searchMethod: "fullSearch" or "dirSearch"
    ##  - A_m: The A_m Matrix described in the paper
    ##  - X_B: The natural Spline Basis
    ##  - maxPoI: number of maximal allowed PoIs (a typical choice is 8) 
    ##  - nbest = 1 for regss function
    ##  - intercept = FALSE for regss function
    ## Output:
    ##  - beta_hat: The CraKneSa - Smoothing spline estimators for functional linear regressions Page 7, Formula 2.6    
    

    N    <- length(Y)
    p    <- length(grd)
    ## Standardize and center y and X
    Y_st   <- scale(Y)
    
    ## X_st              <- scale(X_mat)  # dim(X_st)
    X_st  <- t(scale(t(X_mat)))  # dim(X_st)
    
    
    possible_deltas  <- ind2grd(k_seq, grd)
    
    ## ######################################################################
    ## Main procedure
    res  <- lapply(k_seq, FUN = function(k) { # k = k_seq[1]
        k_c    <- which(k == k_seq)   # counter for filling containers
        
        ## Compute potential PoI: 
        ## Get all possible PoI candidates
        PoISearch   <- searchPotPoi(k, X_st, Y_st, dom=dom, plotting = FALSE) # with standardized values; might be better
        potPoI      <- PoISearch[["ind.x"]]
        cor         <- PoISearch[["cor"]]
        maxS        <- length(potPoI)
        ## Estimate BIC, Beta and PoI Effects in the estOrder (R1 or R2)
        ## paste0("estBetaAndPoI_", estOrder)
        if (center) {
            BICandEstimates <- if (estOrder == "R1") {
                                   estBetaAndPoI_R1_centered(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=potPoI, searchMethod=searchMethod, 
                                                             rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                             nbest = nbest, intercept = intercept, plotting = FALSE)
                               } else {
                                   estBetaAndPoI_R2_centered(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=potPoI, searchMethod=searchMethod, 
                                                             rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                             nbest = nbest, intercept = intercept, plotting = FALSE)
                               }
            
        } else {
            BICandEstimates <- if (estOrder =="R1"){
                                   estBetaAndPoI_R1(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=potPoI, searchMethod=searchMethod, 
                                                    rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                    nbest = nbest, intercept = intercept, plotting = FALSE)
                               } else {
                                   estBetaAndPoI_R2(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=potPoI, searchMethod=searchMethod, 
                                                    rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                    nbest = nbest, intercept = intercept, plotting = FALSE)
                               }
        }
        ## Collect all data for one "k" in a list and add the true DGP Parameters
        BICandEstimates[["cor"]]      <- cor
        BICandEstimates[["k"]]        <- k_seq[k_c]
        BICandEstimates[["delta"]]    <- possible_deltas[k_c]
        BICandEstimates[["nStepsES"]] <- 1
        return(BICandEstimates)
    }) # End of for all k /  deltas save <- res[[r]] 
    names(res)       <- possible_deltas
    
    ## Ex Post Estimation
    ## 2018-11-05: arbitrary repetition of ES step added
    if (exPost){
        if (center){
            exPostEstimation   <- lapply(res, function(entry) {
                if (estOrder == "R1") {
                    nStepsES <- entry[["nStepsES"]]
                    
                    while (nStepsES < maxStepsES){
                        
                        CurSetPotPoI <- entry[["estTauGrd"]]
                        
                        entry <- estBetaAndPoI_R1_centered(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=entry[["estTauGrd"]], searchMethod=searchMethod, 
                                                           rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                           nbest = nbest, intercept = intercept, plotting = FALSE)
                        if (identical(sort(CurSetPotPoI), sort(entry[["estTauGrd"]]) ) ) {
                            ## then additional ES step has nothing changed
                            return(entry)
                        }
                        nStepsES <- nStepsES + 1
                        entry[["nStepsES"]] <- nStepsES
                        oEntry <- entry
                    }
                    return(entry)
                    
                } else {
                    nStepsES <- entry[["nStepsES"]]
                    
                    while (nStepsES < maxStepsES){
                        
                        CurSetPotPoI <- entry[["estTauGrd"]]
                        entry <- estBetaAndPoI_R2_centered(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=entry[["estTauGrd"]], searchMethod=searchMethod, 
                                                           rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                           nbest = nbest, intercept = intercept, plotting = FALSE)
                        if (identical(sort(CurSetPotPoI), sort(entry[["estTauGrd"]]) ) ) {
                            ## then additional ES step has nothing changed
                            return(entry)
                        }
                        nStepsES <- nStepsES + 1
                        entry[["nStepsES"]] <- nStepsES
                    }
                    return(entry)
                }
            })
        }else{
            exPostEstimation   <- lapply(res, function(entry) {
                if (estOrder == "R1") {
                    nStepsES <- entry[["nStepsES"]]
                    
                    while (nStepsES < maxStepsES){
                        
                        CurSetPotPoI <- entry[["estTauGrd"]]
                        entry <- estBetaAndPoI_R1(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=entry[["estTauGrd"]], searchMethod=searchMethod, 
                                                  rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                  nbest = nbest, intercept = intercept, plotting = FALSE)
                        if (identical(sort(CurSetPotPoI), sort(entry[["estTauGrd"]]) ) ) {
                            ## then additional ES step has nothing changed
                            return(entry)
                        }
                        nStepsES <- nStepsES + 1
                        entry[["nStepsES"]] <- nStepsES
                    }
                    return(entry)
                    
                } else {
                    nStepsES <- entry[["nStepsES"]]
                    
                    while (nStepsES < maxStepsES){
                        
                        CurSetPotPoI <- entry[["estTauGrd"]]
                        entry <- estBetaAndPoI_R2(Y=Y, X_mat=X_mat, add.vars = add.vars, N=N, p=p, potPoI=entry[["estTauGrd"]], searchMethod=searchMethod, 
                                                  rho_rng = rho_rng, A_m=A_m, X_B = X_B, grd=grd, maxPoI = maxPoI, 
                                                  nbest = nbest, intercept = intercept, plotting = FALSE)
                        if (identical(sort(CurSetPotPoI), sort(entry[["estTauGrd"]]) ) ) {
                            ## then additional ES step has nothing changed
                            return(entry)
                        }
                        nStepsES <- nStepsES + 1
                        entry[["nStepsES"]] <- nStepsES
                    }
                    return(entry)
                    
                }
            })
        }
        
        for( delta in possible_deltas){ # delta <- possible_deltas[1]
            delta_c <- which(delta == possible_deltas)
            res[[delta_c]][["exPostEstimation"]] <- exPostEstimation[[delta_c]]
        }
    }
    ## Get Optima of via BIC, BIC_sqHat, GCV
    optimum         <- which.min(sapply(res, function(x) x$BIC))
    ## Choose the best (acc. to BIC_sqHat) and save it as "opt_BIC_sqHat"
    optimum_BIC_sqHat   <- which.min(sapply(res, function(x) x$BIC_sqHat))
    ## Choose the best (acc. to gcv) and save it as "opt_gcv"
    optimum_gcv     <- which.min(sapply(res, function(x) x$gcv))
    
    ## Get Ex Post Optima of via BIC, BIC_sqHat, GCV
    ## Choose the best (acc. to ExPost_BIC) and save it as "exPostOptimum"
    if (exPost){
        exPostOptimum           <- which.min(sapply(res, function(x) x[["exPostEstimation"]]$BIC))
        ## Choose the best (acc. to ExPost_BIC_sqHat) and save it as "exPostOptimum_BIC_sqHat"
        exPostOptimum_BIC_sqHat <- which.min(sapply(res, function(x) x[["exPostEstimation"]]$BIC_sqHat))
        exPostOptimum_gcv       <- which.min(sapply(res, function(x) x[["exPostEstimation"]]$gcv))
    }
    ## Save everything
    res[["optimum"]]                <- list(res[[optimum]], X = X_mat,Y = Y)
    names(res[["optimum"]])         <- c("res", "X", "Y")
    res[["optimum_BIC_sqHat"]]      <- res[[optimum_BIC_sqHat]][1:14]
    res[["optimum_gcv"]]            <- res[[optimum_gcv]][1:14]
    if (exPost){
        res[["exPostOptimum"]]            <- res[[exPostOptimum]][["exPostEstimation"]]
        res[["exPostOptimum_BIC_sqHat"]]  <- res[[exPostOptimum_BIC_sqHat]][["exPostEstimation"]]
        res[["exPostOptimum_gcv"]]        <- res[[exPostOptimum_gcv]][["exPostEstimation"]]
    }

    ## #########################################################################################
    ## build the object
    opt <- paste0( if(exPost) "exPostOptimum" else "optimum" , if(opt.crit=="BIC") "" else paste0("_",opt.crit))    
    optInd    <- get(opt)
    optEntry  <- if(opt=="optimum") res[[opt]]$res else res[[opt]]
    betaCurve <- drop(optEntry$estBeta)


    PoI <- optEntry$estTauGrd

    Mat        <- CraKneSaSplineMatrices(X_mat=X_mat, PoI=PoI, add.vars = add.vars, A_m = A_m, p = p)
    X          <- Mat[["X"]]
    A_m_poi    <- Mat[["A_m_poi"]] # dim(A_m)
    
    ##NPXTX      <- 1/(N*p) * t( X) %*%  X
    NPXTX      <- 1/(N*p) * crossprod( X, X )
    hat_matrix <- 1/(N*p) * X %*% chol2inv(chol(NPXTX + optEntry$rho *  A_m_poi )) %*% t( X )

    if (!is.null(add.vars)){
        names(optEntry$betaAddVar) <-
            if (!is.null(colnames(add.vars))) {
                colnames(add.vars)
            } else {
                paste0("add.var", 1:col(add.vars))
            }
    }

    estObj <- list(
        coefficients = list(
            betaCurve = betaCurve,
            betaPoI = optEntry$estPoI,
            betaAddVar = optEntry$betaAddVar,
            tauGrd = optEntry$estTau ,
            tauInd = optEntry$estTauGrd,
            selS = optEntry$selS ) ,
        residuals = Y - drop(optEntry$Y_hat) - if(center) mean(Y) else 0 ,
        fitted = drop(optEntry$Y_hat) + if(center) mean(Y) else 0 ,
        data = list(Y = Y, X_mat = X_mat),
        call = "PESES" ,
        model = list(
            k = res[[ optInd ]]$k ,
            delta = ind2grd(res[[ optInd ]]$k , grd) ,
            rho = optEntry$rho,
            rho_rng = range(rho_rng) ,
            gcv = optEntry$gcv,
            maxPoI = maxPoI,
            cor = res[[ optInd ]]$cor,
            eff_df = optEntry$eff_df,
            bic = optEntry$edfsAndBIC$BIC,
            XtX_1Xt = optEntry$XtX_1Xt,
            bic_sq = optEntry$edfsAndBIC$BIC_sqHat,
            hat_matrix = hat_matrix,
            center = center),
        kSeqEst = res
        )
    estObj
}
