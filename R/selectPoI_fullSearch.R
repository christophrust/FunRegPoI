selectPoI_fullSearch <-
function(Y, X_mat, PoIXValues, potPoI, k = 0, K = 0, maxPoI = 5, nbest = 1, intercept = FALSE, plotting = FALSE , grd = grd){
    ## ########################
    ## A Function returning the choice of PoI given potential candidates
    ## Input:
    ##  - X_mat: N x p Matrix
    ##  - Y: Nx1 vector of dependent variables
    ##  - pot_tau_ind: Potential point of impact indices
    ##  - maxPoI: Maximum size of selected variables
    ##  - nbest: Number of chosen models
    ##  - intercept: Intercept including?
    ##  - plotting: Should a plot be provided?
    ## Output:
    ##  - PoI_Ind: Finally chosen PoIs Indeces
    maxS       <- length( PoIXValues)
    Xy         <- as.data.frame( cbind(PoIXValues, Y) )# put all in a data-frame 
    
    
    if( ncol(Xy) <= 2) {
        if(ncol(Xy) == 2 ){ # if one poi, choose via BIC between one and none
            colnames(Xy)  <- c(paste('TAU',1:length(potPoI),sep=""),"Y")
            minBic        <- which.min( c(one= BIC(lm("Y ~ -1 + TAU1", data = Xy)) , none = BIC(lm("Y ~ -1", data = Xy)) ) )
        } 
        if(ncol(Xy) == 1 || names(minBic) == "none"){     # If none, empty PoI_Ind
            PoIChoice     <- NULL 
        } else if( names(minBic) == "one" ){ # If one, one PoIChoice via potPoI
            PoIChoice     <- potPoI
        }
        
        ## best subset selection only if number of tau-candidates is not too large
    } else if( (ncol(Xy)-1) <= 100 ){
        ## Set colnames
        colnames(Xy)         <- c(paste('TAU',1:length(potPoI),sep=""),"Y")
        ## Best subset selection:
        regss                 <- regsubsets(as.formula(paste("Y ~ ", paste(colnames(Xy)[-dim(Xy)[2]], collapse= "+"))),
                                            data      = Xy,
                                            ## ###################################################
                                            ## Ich habe oben die konstante durch Y-mean(Y) und X-rowMean(X) entfernt.  
                                            ## Daher muessen wir hier das default-intercept entfernt. 
                                            intercept = intercept, # intercept = FALSE
                                            ## ###################################################
                                            nbest     = nbest,         # return only the best selections with 1, 2, ... , PoI_max predictors. 
                                            nvmax     = maxPoI)         # select at most 'PoI_max' taus
        
        ## Use the BIC in order to extract the overall best model 
        ## from the PoI_max different best models with 1, 2, ... , PoI_max predictors. 
        chosen_X_bool   <- summary(regss)$which[which.min(summary(regss)$bic), ]  
        
        PoIChoice       <- potPoI[ chosen_X_bool==TRUE ]
        Xy              <- as.data.frame( cbind(PoIXValues[,chosen_X_bool, drop=F], Y) )# put all in a data-frame 
        colnames(Xy)    <- c(paste('TAU',1:length(PoIChoice),sep=""),"Y")
        
        ## Check again via BIC the best of regsubsets vs. an intercept
        intcOrNoIntc    <- which.min(c(BIC(lm("Y~1", data = Xy)), BIC(lm("Y~-1+.", data = Xy))))
        
        if(plotting){abline(v = grd[ potPoI[chosen_X_bool == TRUE] ], col="darkgreen", lty = 3, lwd = 3)} 
        
        ## If the intercept model is better, return NULL as PoI_Ind
        if(intcOrNoIntc == 1){
            PoIChoice   <- NULL
        }
    }
    
    PoIChoice
}
