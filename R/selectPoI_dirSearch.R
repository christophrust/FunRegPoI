selectPoI_dirSearch <-
function(Y, X_mat, PoIXValues, potPoI, k = 0, K = 0, maxPoI = 5, nbest = 1, intercept = FALSE, plotting = FALSE){
    
    maxS          <- min(length(potPoI),maxPoI) # potPoI = NULL
    estDataFrame  <- createEstDFforDirectSearch(Y, X_mat, PoI_vals = PoIXValues)
    LMandBICs     <- sapply(0:maxS, function(s) calcPoIandScoreLM(estDataFrame, s = s, k = 0, K = 0, S = maxS)$BIC)
    
    S_choice      <- (0:maxS)[ which.min(LMandBICs) ] 

    if(S_choice == 0 ){
        PoIChoice <- NULL
    } else{
        PoIChoice <- potPoI[1:S_choice]
    }
    PoIChoice
}
