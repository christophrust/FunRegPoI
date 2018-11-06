#' FunRegPoISim
#'
#' Data-Simulation following the Functional Linear Model with Points of Impact.
#' \code{FunRegPoISim} generates data following the functional linear model with
#' scalar response with additional points of impact (see Kneip et
#' al. 2016). The random curves are draws of a Brownian Motion process.
#'
#' @param N Size of drawn sample (hence, N curves are drawn as well as N
#'        corresponding observations of the scalar response variable).
#' @param grd The grid where the functions are evaluated.
#' @param DGP One of "Easy" , "Complicated" , "NoPoI" , "OnlyPoI". See details below.
#' @param error_sd Standard deviation of the error term in the functional linear model.
#' @param mean_X Numeric specifying the mean of the brownian motions.
#' @param sd_X Numeric specifying the standard deviation of the Brownian Motion
#'        increments.
#'
#' @details
#' \code{FunRegPoISim} simulates a sample of curves and corresponding
#'   observations of a scalar response variable following the functional
#'   linear model with points of impact. Such data can be used in the
#'   function \code{FunRegPoI} (see example below). The model is given by
#'
#'   \eqn{Y = \int_0^1 \beta(t) X(t) \mathrm{d}t + \sum_{s=1}^S \beta_s X(\tau_s) +
#'     \epsilon.}
#'
#'   The curves are draws from a Brownian Motion and the implementation
#'   currently permits to simulate four differend DGPs, named "Easy", "Complicated",
#'   "NoPoI" and "OnlyPoI". Their specification is as follows:
#'   \describe{
#'     \item{Easy}{
#'       The beta function of the functional part is a polynomial of order 3:
#'       \eqn{\beta(t) = -(t-1)^2 + 2} together with two points of impact
#'       which are located at 0.3 and 0.6 with coefficients -3 and 3.
#'     }
#'     \item{Complicated}{
#'       The beta function of the functional part is a polynomial of order 4:
#'       \eqn{\beta(t) = -5(t-0.5)^3 - t + 1} together with three points of
#'       impact are located at 0.3, 0.4 and 0.6 with coefficients -3, 3 and 3.
#'     }
#'     \item{NoPoI}{
#'       The beta function of the functional part is the same as in
#'       \code{Easy} but the DGP does not include any further points of impact.
#'     }
#'     \item{OnlyPoI}{
#'       Only the points of impact of \code{Easy} are included.
#'     }
#'   }
#'
#' @return
#' A list with entries
#'   \item{Y}{
#'     A numeric vector of length \code{N} with the scalar outcome.
#'   }
#'   \item{X}{
#'     A \eqn{p \times N} matrix holding trajectories of the curves.
#'   }
#'   \item{trueBetaPoI}{
#'     A vector with the true values of the PoI-coefficients.
#'   }
#'   \item{trueBetaCurve}{
#'     A vector of length \eqn{p} with the evaluated beta-curve.
#'   }
#'   \item{trueTauGrd}{
#'     A vector with true values of the PoI-locations.
#'   }
#'   \item{trueTauInd}{
#'     A vector with indices of the true PoI-locations in \code{grd}.
#'   }
#'
#' @references
#' Kneip A., Poss, D., Sarda, P. (2016) Functional Linear
#' Regression with Points of Impact. The Annals of Statistics, \emph{44}(1), 1-30.
#'
#' @author Dominik Liebl, Stefan Rameseder and Christoph Rust
#'
#' @seealso FunRegPoI
#' @examples
#' \dontrun{
#' library(FunRegPoI)
#'
#' ## Define Parameters for Simulation
#' DGP_name    <- "Easy"
#' k_seq       <- c(1,seq(2, 62, 4))
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
#' # PESES
#' EstPeses    <- FunRegPoI(Y = PoISim$Y , X_mat = PoISim$X , grd, estimator = "PESES", k_seq = k_seq)
#'
#' # KPS
#' EstKPS      <- FunRegPoI(Y = PoISim$Y , X_mat = PoISim$X , grd, estimator = "KPS", k_seq = k_seq)
#'
#' # CKS
#' EstCKS      <- FunRegPoI(Y = PoISim$Y , X_mat = PoISim$X , grd, estimator = "CKS")
#'
#' par(mfrow=c(3,1))
#' plot(EstPeses, main = "PESES")
#' lines(x = grd , y = PoISim$trueBetaCurve , col ="red")
#' abline( v = PoISim$trueTauGrd , col ="red" , lwd = 2)
#' plot(EstKPS, main = "KPS")
#' lines(x = grd , y = PoISim$trueBetaCurve , col ="red")
#' abline( v = PoISim$trueTauGrd , col ="red" , lwd = 2)
#' plot(EstCKS, main = "KPS")
#' lines(x = grd , y = PoISim$trueBetaCurve , col ="red")
#' abline( v = PoISim$trueTauGrd , col ="red" , lwd = 2)
#' dev.off()
#' }
#' 
#' @keywords datagen
#'
#' @export
FunRegPoISim <- function(N, grd, DGP = c("Easy", "Complicated", "NoPoI", "OnlyPoI", "EasyHighNoise"),
                         error_sd = if (identical(DGP, "EasyHighNoise")) 0.5 else 0.125, mean_X = 0 , sd_X = 1) {
    
    ## grd params and deviance
    p        <- length(grd)
    a        <- min(grd)
    b        <- max(grd)
    mean_Y   <- 0		# Mean of Y
    sd_Y     <- error_sd        # Sd of scalar Y errors
    mean_X   <- 0		# Mean of X
    sd_X     <- 1    	        # Sd of functional X 

    ## set the PoI-parameters (tau, PoI-coefs , beta function)
    if (identical(DGP, "EasyHighNoise")) DGP <- "Easy"
    PoIDGP     <- setPoIDGP(DGP)
    beta_fun   <- function(t){ eval(parse(text=PoIDGP$fct_text))}
    beta_fun_grd <- beta_fun(t=grd)
    tau_ind_true <- grd2ind(PoIDGP$tau, grd)
    
    
    X_mat <- sapply( 1:N , function(i) {
        dt <-(b-a)/(p-1)
        W  <- cumsum(c(rnorm(1,0,0), rnorm(p-1 + floor(0.2*p), 0, sqrt(dt) )))
        return(W[(length(W)-(p-1)):length(W)])
    })
    
    ## Evaluate the X functions at the true taus
    X_tau <- t(X_mat[tau_ind_true, 1:N, drop=F]) # N x length(tau)
    
    ## PoI Model: Y_i = beta_0 + int X_i(t) beta(t) + sum_j beta_j X_i(tau_j) + eps_i
    Y     <-  t(X_mat) %*% beta_fun_grd/p + X_tau %*% PoIDGP$beta_tau + rnorm(N, mean = mean_Y, sd = sd_Y)
    
    simObj <- list(Y = Y ,
                   X = X_mat,
                   trueBetaPoI = PoIDGP$beta_tau,
                   trueBetaCurve = beta_fun_grd,
                   trueTauGrd = PoIDGP$tau,
                   trueTauInd = tau_ind_true
                   )
    class(simObj) <- "FunRegPoISim"
    simObj
}
