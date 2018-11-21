estBetaCraKneSa <-
function(Y, X_mat, add.vars, A_m, X_B, PoI = NULL, rho_rng = c(1e-6,2e2) , ...){

    ## ########################
    ## A Function estimating the beta(t)-Curve for a functional linear regression model
    ## Input:
    ##  - rho: smoothing parameter; rho <- 5e-03
    ##  - S: number of PoIs
    ##  - Y: dependent scalar values
    ##  - X_mat: independent functional values; p x N
    ##  - A_m: Second derivatives of natural splines basis depending on order m of splines
    ##  - X_B: X matrices of Splines basis depending on order m of splines
    ##  - N: number of curves
    ##  - p: eval points
    ##  - PoI = PoIChoice
    ## Output:
    ##  - beta_hat: The CraKneSa - Smoothing spline estimators for functional linear regressions Page 7, Formula 2.6

       
    N <- length(Y)
    p <- nrow(X_mat)
    
    Mat     <- CraKneSaSplineMatrices(X_mat=X_mat, PoI= PoI, add.vars = add.vars, A_m=A_m, p=p)
    X       <- Mat[["X"]]
    A_m_poi <- Mat[["A_m_poi"]] # dim(A_m)
    
    ## NPXTX   <- 1/(N*p) * t( X) %*%  X
    NPXTX   <- 1/(N*p) * crossprod( X , X)
       
    ## Use faster Cholesky inversion since NPXTX + x *  A_m_poi is symmetric for all x
    optGCVviaRho     <- function(x) {
        ##H <- 1/(N*p) * X %*% chol2inv(chol(NPXTX + x *  A_m_poi )) %*% t( X ) # dim(H) N x N
        H <- 1/(N*p) * tcrossprod( X %*% backsolve(chol(NPXTX + x *  A_m_poi ), diag(1,ncol(A_m_poi))) )
        return((1/N *  sum ( (Y - H %*% Y)^2 ))/ (( 1 - sum( diag(H)) / N)^2 ))
    }


    
    optRhoGivenPoI   <- optim(1e-3, optGCVviaRho, method = "Brent", lower = range(rho_rng)[1], upper = range(rho_rng)[2])
    
    ## Extract optimal rho
    rho              <- optRhoGivenPoI$par
    
    ## Complete estimator
    ## Alpha Estimation (alpha = discretized version of beta(t) + PoI coefficients (last rows)
    ## CraKneSa - Smoothing spline estimators for functional linear regressions
    ## Page 7, Formula 2.6
    ## alpha_hat is a p+S vector with p points for beta(t) and S estimations of the PoIs
    ## XtX_1Xt is the projection of Y onto alpha_hat WITHOUT PoI    
    ## 2. b. Generate corresponding matrices
    betaEstimates          <- projEstimatorMatrix( X = X, X_B, Y, A_m_poi = A_m_poi, rho, N , p)
    betaEstimates[["rho"]] <- rho
    betaEstimates[["gcv"]] <- optRhoGivenPoI$value

    betaEstimates
}
