projEstimatorMatrix <-
function( X, X_B, Y, A_m_poi, rho, N , p) {
    
    ##XtX_1Xt   <- 1/N * solve( 1/(N * p) * t(X) %*% X + rho * A_m_poi ) %*% t(X) # nur auf X fuer beta kurve und nicht auf die anderen
    XtX_1Xt    <- 1/N * chol2inv( chol( 1/(N * p) * crossprod(X,X) + rho * A_m_poi ) ) %*% t(X) # nur auf X fuer beta kurve und nicht auf die anderen
    alpha_hat  <- XtX_1Xt %*% Y
    
    ## Beta(t): Projection of alpha_hat onto Spline Basis
    ## P(P^T P)^-1 P^T alpha_hat
    ##beta_hat   <- X_B %*% solve( t(X_B) %*% X_B) %*% t(X_B) %*% alpha_hat[1:p]
    ##beta_hat   <- X_B %*% chol2inv( chol( crossprod(X_B, X_B) )) %*% t(X_B) %*% alpha_hat[1:p]
    beta_hat   <- tcrossprod( X_B %*% backsolve(chol(crossprod(X_B, X_B)), diag(1,ncol(X_B))) ) %*% alpha_hat[1:p]
    
    list(estBeta=beta_hat, XtX1Xt = XtX_1Xt, estAlpha = alpha_hat)
}
