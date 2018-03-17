#' calSecDerNatSpline
#'
#' Natural Spline Basis and the Regularization Matrix \code{A_m}
#'
#' \code{calSecDerNatSpline} calculates evaluated natural cubic splines
#' and the matrix \code{A_m} as described in Crambes, Kneip and Sarda
#' (2016).
#'
#' @param grd A vector of length \eqn{p} indicating the grd values where the
#'        spline functions are evaluated.
#'
#' @details
#' \code{calSecDerNatSpline} calculates an evaluated natural cubic
#' spline basis corresponding to \code{grd} and a \eqn{p \times p} matrix
#' used for regularization in the smoothing splines estimator for the
#' functional linear model (Crambes, Kneip and Sarda,2016). \code{A_m} is
#' composed of a non-classical projection matrix and the classical
#' regularization matrix of squared second derivatives.
#'
#' @return
#' A list with entries
#' \item{A_m}{
#'   The \eqn{p \times p} matrix \eqn{A} as described in (Crambes, Kneip and Sarda,
#'   2016)
#' }
#' \item{X_B}{
#'   A \eqn{p \times p} matrix of natural cubic spline basis evaluated at \code{grd}.
#' }
#'
#' @references
#' Crambes, C., Kneip, A., Sarda, P. (2009) Smoothing Splines
#' Estimators for Functional Linear Regression. The Annals of
#' Statistics, \emph{37}(1), 35-72.
#'
#' @author Dominik Liebl, Stefan Rameseder and Christoph Rust
#'
#' @seealso FunRegPoI
#'
#' @examples
#' \dontrun{
#' library(FunRegPoI)
#'
#' ## Define Parameters for Simulation
#' DGP_name    <- "Easy"
#' k_seq       <- c(1,seq(2, 60, 4))
#' error_sd    <- 0.125
#' N           <- 500
#' p           <- 300
#' domain      <- c(0,1)
#' grd         <- seq(domain[1],domain[2],length.out = p)
#'
#' set.seed(1)
#' ## simulate some data:
#' PoISim      <- FunRegPoISim(N = N , grd, DGP = DGP_name , error_sd = error_sd)
#'
#' NaturalSplines <- calSecDerNatSpline(grd)
#'
#' # PESES
#' EstPeses  <- FunRegPoI(Y = PoISim$Y , X_mat = PoISim$X , grd, estimator = "PESES",
#'                        k_seq = k_seq , A_m = NaturalSplines$A_m , X_B = NaturalSplines$X_B)
#' 
#' }
#'
#' @keywords models
#' 
#' @export
calSecDerNatSpline <-
function(grd) {
    ## ########################
    ## CraKneSa - Smoothing spline estimators for functional linear regressions
    ## P. 6 ff, 
    ## Function for calulcating second derivates of natural m splines
    ## will be used in modelSelect
    ## Input:
    ##  - m: spline order according to CraKneSa
    ##  - grd: Grid of X(t)
    ##  - p: eval points
	
    ## Output:
    ##  - A_m: matrix of second derivatives
    p <-length(grd)
    if (any(!(range(grd)!=c(0,1)))) grd <- seq(0,1,len = length(grd)) # transform if grd is not a seq over [0,1]
    order 	<- 4
    grid.len.2d <- 10* p
    
    ## Generate the B-spline basis matrix for a natural spline  ns() function of splines package
    ## X_B1 	<- ns( x = grd , knots = grd[2:(p-1)] ,intercept = TRUE)
    
    ## Generate polynomial basis of degree m-1
    polynom <- paste( paste("grd^" , c(0:1) , sep="") , collapse =" , ")
    polynom <- paste("cbind(" , polynom , ")" , sep="")
    polymat <- eval(parse( text=polynom))
    
    ## Hat Matrix with evaluated polynomial entries
    P_mat 	<- polymat %*% solve( t(polymat) %*% polymat ) %*%  t(polymat)
    
    ## fine grid for evaluation of second derivatives
    grd_2nd <- seq(0,1,length=grid.len.2d)[1:(grid.len.2d-1)] + 1/(2*grid.len.2d)

    ## all knots, for splineDesign
    Aknots <- c( rep(range(grd) , order) , grd[2:( length(grd)-1)])

    ## evaluated b-spline-basis
    B    <- splineDesign(knots = Aknots , x = grd , ord = order , derivs = 0)
    B_d2 <- splineDesign(knots = Aknots , x = grd_2nd , ord = order , derivs = rep(2,length(grd_2nd)) )

    ## impose knot-condition (zero second derivative at boundary knots)
    C     <- splineDesign(knots = Aknots , x = range(grd_2nd) , ord = order , derivs = c(2,2))
    qr.c  <- qr(t(C))
    X_B   <- as.matrix( (t( qr.qty( qr.c , t(B))))[,-(1:2) , drop = FALSE])
    BB_d2 <- as.matrix( (t( qr.qty( qr.c , t(B_d2))))[,-(1:2) , drop = FALSE])
    
    
    G_mat <- t(BB_d2) %*% BB_d2 * 1/length(grd_2nd)
    
    ## btb      <- solve( t(X_B) %*% X_B )
    btb      <- chol2inv( chol( crossprod(X_B , X_B )))
    
    A_star_m <- X_B %*% btb %*% G_mat %*% btb %*% t(X_B)

    A_m      <- P_mat + p * A_star_m

    list(A_m = A_m, X_B = X_B)
}
