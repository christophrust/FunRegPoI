createEstDFforDirectSearch <-
function(Y, X, PoI_vals){
    ## ########################
    ## A Function creating the DF for the DirectSearch Function
    ## Input:
    ## -  Y
    ## -  X
    ## -  PoI_vals
    Y_names <- "Y"
    
    ## if no PoIs
    if( length(PoI_vals) == 0 ){ 	
        if( is.null(X) ){ 			# and if no regressors
            estDataFrame  <- data.frame(Y)
            return(estDataFrame)
        } else { 					# and if regressors
            X_names             <- paste("PC", 1:dim(X)[2], sep="")
            estDataFrame        <- data.frame( cbind(Y , X) )
            names(estDataFrame) <- c( Y_names , X_names  )
            return(estDataFrame)
        }
	## if exactly one PoI
    } else if (length(PoI_vals) == length(Y) ){ 
        if( is.null(X) ){ 			# and if no regressors
            PoI_names           <- paste("PoI.1")
            estDataFrame        <- data.frame(cbind(Y , PoI_vals))
            names(estDataFrame) <- c( Y_names , PoI_names  )
            return(estDataFrame)
        } else {  					# and if regressors
            X_names             <- paste("PC", 1:dim(X)[2], sep="")
            PoI_names           <- paste("PoI.1")
            estDataFrame        <- data.frame(cbind(Y , X, PoI_vals))
            names(estDataFrame) <- c( Y_names , X_names, PoI_names  )
            return(estDataFrame)
        }
        ## if multiple PoIs 
    } else if (dim(PoI_vals)[1] == length(Y) ){ 
        if( is.null(X) ){ 			# and if no regressors
            PoI_names           <- paste("PoI", 1:dim(PoI_vals)[2], sep=".")
            estDataFrame        <- data.frame(cbind(Y , PoI_vals))
            names(estDataFrame) <- c( Y_names , PoI_names  )
            return(estDataFrame)
        } else {  					# and if regressors
            if(length(X) == length(Y)){ # if 1 regressors
                X_names <- "PC1"
            } else{
                X_names  <- paste("PC", 1:dim(X)[2], sep="")
            }
            PoI_names    <- paste("PoI", 1:dim(PoI_vals)[2], sep=".")
            estDataFrame <- data.frame(cbind(Y , X, PoI_vals))
            names(estDataFrame) <- c( Y_names , X_names, PoI_names  )
            return(estDataFrame)
        }	
    } else {
        stop("Something is wrong with the dimensions in X and Y")
    }
}
