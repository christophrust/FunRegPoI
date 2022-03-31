#' FunRegPoI
#' 
#' Functional Linear Regression with Points of Impact
#'
#' \code{FunRegPoI} estimates the functional linear regression model with
#' points of impact. It offers the usage of different estimators, such as
#' the originally suggested FPCA-based approach or an estimation
#' approach based on smoothing splines. In particular, the PESES approach
#' as proposed in Liebl, Rameseder, Rust (2017), is implemented.
#'
#' @param Y A vector of length \eqn{N} containing the scalar dependent variable.
#' @param X_mat A \eqn{p \times N}-matrix holding the functional predictors, where \eqn{p}
#'        is the number of grid points provided by \code{grd} and \eqn{N} is
#'        the number of observations.
#' @param grd A vector of length \eqn{p} giving the location of the (equidistant)
#'        grid where the functional observations are observed.
#' @param add.vars Additional explanatory variables (optional, e.g. seasonal dummies or further
#'        scalar valued covariates, supplied as a matrix.
#' @param estimator Character, one of "PESES" , "KPS" (Kneip, Poss and Sarda, 2016), "CKS"
#'        (Crambes, Kneip and Sarda, 2016), see details.
#' @param maxPoI Maximum number of PoIs.
#' @param A_m An optional \eqn{p \times p} matrix holding the scalar product of
#'        the basis function's second derivatives. Specifying this argument
#'        togehter with the \code{X_B} only makes sense in situations, where
#'        the estimation is performed many times on the same domain
#'        (\code{grd}) in order to speed up the calculation, as for instance
#'        in simulation studies. This matrix can be obtained via the function
#'        \code{calSecDerNatSpline}.
#' @param X_B An optional \eqn{p \times p} matrix of the natural cubic spline
#'        basis evaluated at \code{grd}, provided together with
#'        \code{A_m}. This matrix can be obtained via the function
#'        \code{calSecDerNatSpline}.
#' @param \dots Further arguments passed to \code{optim()}.
#'
#' @details
#' \code{FunRegPoI} estimates the functional linear regression model with points
#'  of impact (Kneip, Poss, Sarda, 2010). In addition to the proposed
#'  FPCA-based estimator, an estimation procedure based on  smoothing
#'  splines is available. This approach is based on the work of Crambes, Kneip and Sarda
#'  (2009) which is extended to the situation with points of impact in
#'  Liebl, Rameseder, Rust (2017). The default estimator in this situation
#'  is called \code{"PESES"} and implements the proposed estimation
#'  procedure which consists in a sequential selection procedure, see
#'  Liebl, Rameseder, Rust (2017) for a more detailed description:
#'  \enumerate{
#'    \item{Preselect potential points of impact.}
#'    \item{Estimate a \eqn{\beta}-curve and PoI coefficients given the potential points of impact.}
#'    \item{Sub-Select relevant points of impact.}
#'    \item{re-Estimate a final \eqn{\beta}-curve and PoI cofficients given the selected points of impact.}
#'    \item{re-Select the final set of points of impact.}
#'  }
#'
#' @return
#'  Object of class \code{"FunRegPoI"}:
#' \item{coefficients}{
#'    A list holding the estimated \eqn{\beta}-curve (\code{betaCurve}), the PoI coefficients
#'    (\code{betaPoI}), coefficients for additional variables if provided (\code{betaAddVar}),
#'    the PoI's location in the domain (\code{tauGrd})
#'    transformed to the unit interval with corresponding
#'    \code{grd}-indices (\code{tauInd}) and the number of PoIs (\code{selS}).
#'  }
#'  \item{residuals}{
#'    A vector of length \eqn{p} with residuals of the fit.
#'  }
#'  \item{fitted}{
#'    A vector of length \eqn{p} with fitted values of the estimate.
#'  }
#'  \item{call}{
#'    Type of  the estimator (PESES , CKS, or KPS).
#'  }
#'  \item{model}{
#'    A list containing the relevant model parameters, such as \code{k},
#'    \code{rho}, \code{gcv}, covariance function between \code{Y} and
#'    decorrelated \code{X}, effective degrees of freedom, \code{BIC} and
#'    the hat matrix.
#'  }
#'  \item{kSeqRes}{
#'    A list of same length as \code{k_seq} holding the estimation result for
#'    each \code{k} in \code{k_seq}. These estimation results are used
#'    during estimation procedure to find the optimal fit with respect to
#'    the selected criterium.
#'  }
#'
#' @references
#' Crambes, C., Kneip, A., Sarda, P. (2009) Smoothing Splines
#' Estimators for Functional Linear Regression. The Annals of
#' Statistics, \emph{37}(1), 35-72.
#'
#' Kneip A., Poss, D., Sarda, P. (2016) Functional Linear
#' Regression with Points of Impact. The Annals of Statistics, \emph{44}(1), 1-30.
#'
#' Liebl, D., Rameseder, S., Rust, C. (2018) Improving Estimation in
#'  Functional Linear Regression with Points of Impact: Insights into
#'  Google AdWords. Available at www.dliebl.com/files/Liebl_Rameseder_Rust_2017.pdf
#' 
#' @author Dominik Liebl, Stefan Rameseder and Christoph Rust
#'
#' @note
#' The estimation procedure retransforms any \code{grd} covering a domain
#' different from \eqn{[0,1]} onto the unit interval. Also the estimation
#' output is defined on this interval and has to be interpreted
#' accordingly. Of course, this has no effect on the estimates of the
#' points of impact coefficients, but their location also is transformed
#' to the unit interval.
#'
#' @seealso \code{FunRegPoISim}, \code{calSecDerNatSpline}
#'
#' @examples
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
#'
#' abline( v = PoISim$trueTauGrd , col ="red" , lwd = 2)
#' dev.off()
#'
#' @keywords models
#'
#' 
#' @importFrom leaps regsubsets
#' @importFrom stats lm as.formula BIC optim prcomp residuals sd splinefun weighted.residuals rnorm
#' @importFrom graphics abline plot lines text par
#' @import splines
#'
#' 
#' @export
FunRegPoI <-
  function(Y, X_mat, grd, add.vars = NULL,
           estimator = c("PESES", "KPS", "CKS"), maxPoI= 8, A_m, X_B, ...) {
### Functional Linear Regression with Points of Impact
### wrapper around different estimators: PESES, KPS , CKS
###
### returns object of S3-class FunRegPoI, with methods plot and summary

    if (missing(estimator)) estimator <- "PESES"
    ## natural cubic splines basis and squared second derivatives if not provided
    if (estimator %in% c("PESES", "CKS") && (missing(A_m) | missing(X_B))) {
      NaturalSplines <- calSecDerNatSpline(grd)
      A_m            <- NaturalSplines[["A_m"]]
      X_B            <- NaturalSplines[["X_B"]]
    }

    Y <- as.vector(Y)

    ## check for correct input
    if (!is.matrix(X_mat)) stop("'X_mat' has to be a matrix!")
    if (!is.null(add.vars) && !is.matrix(add.vars) && !is.vector(add.vars))
      stop("'add.vars' has to be a matrix or vector!")
    if (length(Y) != dim(X_mat)[2])
      stop("dimension of 'X_mat' and length of 'Y' don't match!")
    if (abs(diff(range(diff(grd)))) > 1e-10)
      stop("'grd' does not appear to be a regularly spaced grid!")


    ## check which estimator
    if (estimator == "PESES") {
      ## Apply the PESES-Estimator
      estObj <- CraKneSaPoIEst(Y, X_mat, grd, add.vars, A_m, X_B, maxPoI, ...)
      ##
    } else if (estimator == "KPS") {
      estObj <- KnePosSarEstimation(Y, X_mat, grd, add.vars, maxPoI, ...)
    } else if (estimator == "CKS") {
      estObj <- CraKneSaEst(Y, X_mat, add.vars, A_m, X_B, ...)
    } else {
      stop(sprintf("Estimator %s not supported", estimator))
    }
    class(estObj) <- "FunRegPoI"
    estObj
  }
