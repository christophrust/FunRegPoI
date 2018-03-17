selectEFwithCumPropVar <-
function(X_mat, threshold_cpv = 0.95){
    ## ########################
    ## A Function definining given a DGP_name the PoI Setup und the Beta(t) Setup
    ## Input:
    ##  - X_mat: independent functional values of N functions observed at p points; p x N 
    ##  - threshold_cpv: To choose a number K of eigenfunctions, how much cum. proportion of variance do they have to cover
    ## Output:
    ##  - K: Number of Eigenvectors to choose with the threshold
    ##  - evals: the eigenvalues, p
    ##  - evecs: the eigenvectors, p x p
    ##  - chosenEvecs: the first K eigenvectors, p x K
    
    ## PCA
    p   <- nrow(X_mat)
    pca	<- prcomp(t(X_mat) )
    ## Adjust Eigenvectors/Eigenfunctions
    evecs <- pca$rotation * sqrt(p) # adjustment (approx) in order to have length 1 in L^2
    
    ## To check for the length of Eigenvectors
    ## apply(evecs, 1, FUN = function(row) return(sum((row-mean(row))^2)/p))
    
    ## Calc eigenvalues
    evals <- pca$sdev^2
    
    cpv   <- sapply(1:length(evals), FUN=function(k) return( sum(evals[1:k])/sum(evals)) )
    if( is.numeric(threshold_cpv) ) {
        K  <- which.max( cpv > threshold_cpv) 
    } else if( isTRUE(!threshold_cpv) ){
        K  <- dim(evecs)[2]
    }
    list( K = K, evals = evals, evecs = evecs, chosenEvecs = evecs[ ,1:K])
}
