#' Bootstrapping estimated change points
#' 
#' Obtain bootstrap replicates of changepoint estimators 
#' @param mcpts an object of class \code{mosum.cpts}, \code{multiscale.cpts}
#' or \code{multiscale.bottomUp.cpts}
#' @param N_reps number of bootstrap replications
#' @param alpha numeric value in (0,1), such that the (1-alpha_CI)-confidence bootstrap intervals are computed
#' @return object of class \code{cpts.bootstrap}, containing the following fields:
#'    \item{CI}{data frame of five columns, containing the estimated changepoints (column \code{cpt}),
#'    the pointwise alpha confidence intervals (columns \code{pCI_l} and \code{pCI_r})
#'    and the uniform alpha confidence intervals (columns \code{uCI_l} and \code{uCI_r})}
#'    \item{bootstrap_replicates}{matrix of dimension N_reps times q (where q denotes the estimated 
#' number of changes in mcpts), containing
#' N_reps bootstrapped replicates of the estimated changepoint locations}
#' @details See the referenced literature for further details
#' @references A. Meier, C. Kirch and H. Cho.
#' \emph{mosum: A Package for Moving Sums in Change Point Analysis.}
#' Unpublished manuscript, 2018+.
#' @references H. Cho and C. Kirch.
#' \emph{Multiple change-point detection via multiscale MOSUM procedure with localized pruning.}
#' Unpublished manuscript, 2018+.
#' @examples 
#' set.seed(1337)
#' x <- piecewiseStationary_timeSeries(lengths=rep(100,3), means=c(0,3,1), sds=rep(1,3))
#' mcpts <- mosum.cpts(x, G=40)
#' cpts_boot <- cpts.bootstrap(mcpts, 5000)
#' print(cpts_boot$CI)
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export
cpts.bootstrap <- function(mcpts, N_reps, alpha=0.05) {
  if (class(mcpts)=="mosum.cpts") {
    bstrp <- mosum.cpts.bootstrap(mcpts, N_reps, alpha)
  } else if (class(mcpts)=="multiscale.cpts") {
    bstrp <- multiscale.cpts.bootstrap(mcpts, N_reps, alpha)
  } else {
    stop("mcpts object must be of class mosum.cpts or multiscale.cpts.")
  }
  return(bstrp)
}

#' Bootstrapping MOSUM changepointss
#' @keywords internal
mosum.cpts.bootstrap <- function(mcpts, N_reps, alpha) {
  stopifnot(class(mcpts) == "mosum.cpts")
  if (mcpts$bootstrap) {
    return(mcpts$cpts_bootstrap)
  } else {
    x <- mcpts$m$x 
    G <- mcpts$m$G
    cpts <- mcpts$cpts
    q <- length(cpts)
    cpts.info = matrix(NA, ncol=3, nrow=q)
    cpts.info[,1] <- cpts
    cpts.info[,2] <- rep(G, q)
    if (mcpts$m$symmetric) {
      cpts.info[,3] <- rep(G, q)
    } else {
      cpts.info[,3] <- rep(mcpts$m$G.right, q)
    }
    cpts_bootstrap(cpts.info, x, N_reps, alpha, mcpts)
  }
}

#' Bootstrapping multiscale MOSUM changepoints
#' @keywords internal
multiscale.cpts.bootstrap <- function(ms, N_reps, alpha) {
  stopifnot(class(ms)=="multiscale.cpts")
  if (ms$bootstrap) {
    return(ms$cpts_bootstrap)
  } else {
    x <- ms$x
    cpts_bootstrap(ms$cpts.info, x, N_reps, alpha, ms)
  }
}

#' Helping/wrapper fuction for C++ calls
#' @importFrom stats quantile
#' @keywords internal
cpts_bootstrap <- function(cpts_info, x, N_reps, alpha, mcpts) {
  q <- nrow(cpts_info)
  if (q==0) {
    cpts <- pCI_l <- pCI_r <- uCI_l <- uCI_r <- integer(0)
    k_star <- matrix(NA, nrow=N_reps, ncol=0)
  }
  if (q>0) {
    n <- length(x)
    cpts <- cpts_info[,1]
    # Get bootstrap replicates from C++
    tmp <- cpts_bootstrap_help(cpts_info, x, N_reps)
    k_star <- tmp$k_star
    # Pointwise confidence intervals
    C_value_j <- apply(abs(tmp$k_star1), 2, quantile, 1-alpha/2)
    pCI_l <- pmax(1, cpts - C_value_j)
    pCI_r <- pmin(n, cpts + C_value_j)
    # Uniform confidence intervals
    uCI_help <- apply(abs(tmp$k_star2), 1, max)
    C_value <- quantile(uCI_help, 1-alpha)
    uCI_l <- uCI_r <- rep(NA, q)
    for (j in seq_len(q)) {
      uCI_l[j] <- max(1,
                      cpts[j] - C_value * tmp$sigma2_hat[j] / tmp$d_hat[j]^2)
      uCI_r[j] <- min(n,
                      cpts[j] + C_value * tmp$sigma2_hat[j] / tmp$d_hat[j]^2)
    }
  }
  CI <- data.frame(cpt=cpts, 
                   pCI_l=pCI_l,
                   pCI_r=pCI_r,
                   uCI_l=uCI_l,
                   uCI_r=uCI_r)
  return(structure(list(mcpts=mcpts,
                        N_reps=N_reps,
                        alpha=alpha,
                        CI=CI,
                        bootstrap_replicates=k_star),
              class="cpts.bootstrap"))
}
