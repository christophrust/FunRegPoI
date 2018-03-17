CraKneSaSplineMatrices <-
function(X_mat, add.vars, PoI, A_m, p){
    ## CraKneSa - Generates the Matrices of Page 16, Formula 3.10 
    ## Input:
    ##  - X_mat: data matrix (p x N)
    ##  - PoI_Ind: Point of Impact Index vector with S = number of points
    ##  - A_m: Penalty Matrices of used Splines (p x p)
    ##  - p: number of realization points
    ## Output:
    ##  - X: the PoI expanded and transposed matrix ( (p+S) x N)
    ##  - A_m_poi: the PoI zero expanded matrix ( (p+S) x (p+S)) with a S x S Zero Matrix bottom right

    nav <- if (!is.null(add.vars)) dim(add.vars)[2] else 0

    if (0 < (S <- length(PoI))) { 
        ## Expand X by PoI Influence; scale it by p because it will be downscaled again
        X       <- if (is.null(add.vars)){
                       cbind(t(X_mat) , p * t(X_mat[PoI, , drop = FALSE])) # N x (p + S)-Matrix
                   } else {
                       cbind(t(X_mat) , p * t(X_mat[PoI, , drop = FALSE]), p * add.vars) # N x (p + S)-Matrix
                   }
        
        ## Expand Smoothing Matrix by S Zero Rows and Columns
        A_m_med <- cbind(A_m , matrix(0 , nrow = p , ncol = S +nav )) # A_m a p x p; A_m_med a p x (p+S) 
        A_m_poi <- rbind(A_m_med , matrix(0 , nrow = S + nav, ncol = (S + p + nav))) #  A_m_poi a (p+S) x (p+S) 
    } else {
        X       <- if (is.null(add.vars)) {
                       t(X_mat)
                   } else {
                       cbind(t(X_mat), p * add.vars)
                   }
        A_m_poi <- if (is.null(add.vars)) {
                       A_m
                   } else {
                       A_m_med <- cbind(A_m , matrix(0 , nrow = p , ncol = nav )) # A_m a p x p; A_m_med a p x (p+S) 
                       A_m_poi <- rbind(A_m_med , matrix(0 , nrow = nav, ncol = (p + nav)))
                   }
    } # dim(A_m_med ); dim(A_m_poi); dim(X); dim(X_mat); dim(t(X)); rho = 1 
    list(X = X, A_m_poi = A_m_poi)
}
