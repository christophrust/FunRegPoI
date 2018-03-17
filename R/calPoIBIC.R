calPoIBIC <-
function( X, XtX_1Xt, Y, Y_hat, rho, S, p) {
    ## #################################################
    ## Calculates the BIC of a PoI-Estimate using the
    ## trace of both Hat matrix and the squared Hat matrix
    ##
    
    ##
    ## 4.a. Fitted values and RSS:
    # Y_i_hat = int ( X_i(t) beta(t) dt + sum beta_j X_i(tau_j)
    U_hat  <- Y - Y_hat
    RSS    <- sum(U_hat^2)
    N      <- length(Y)

    # BIC = N * log(SSR/N) + log(N) * trace(SmoothingMatrix)
    #     = N * log( <Y-Y_hat, Y-Y_hat> / N) + log(N) * (Tr(H(X,rho))) where Y_hat = H(X,rho) Y
    # and therefore H(X, rho) = 
    # - the first p columns:    int X(t) beta(t) dt +
    # - the last S columns for one parameter each: S
    eff_df           <- ( sum(diag(( X[ ,1:p] %*% XtX_1Xt[1:p, ])))/ p ) + S
    BIC              <- N * log(RSS/N) + log(N) * eff_df
    eff_df_sqHat     <- ( sum(diag(( X[ ,1:p] %*% XtX_1Xt[1:p, ] %*% X[ ,1:p] %*% XtX_1Xt[1:p, ])))/ p^2 ) + S
    BIC_sqHat        <- N * log(RSS/N) + log(N) * eff_df_sqHat
    
    list(eff_df = eff_df,
         eff_df_sqHat = eff_df_sqHat,
         BIC = BIC,
         BIC_sqHat = BIC_sqHat)
}
