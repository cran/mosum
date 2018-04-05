#' Multiscale MOSUM algorithm (bottom-up)
#' 
#' Multiscale MOSUM procedure with bottom-up bandwidth based merging
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G a symmetric bandwidth grid;
#' either an object of type \code{multiscale.grid} or an integer vector
#' of bandwidths
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, \code{mosum.criticalValue} is used for each bandwidth (with significance
#' level \code{alpha}). Alternatively, it is possible to parse a user-defined function 
#' with \code{threshold.function};
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold="critical.value"}
#' @param threshold.function function object of form \code{function(G)}, to compute a
#' threshold of significance for different bandwidths G; use iff \code{threshold="custom"}
#' @param eta see \link[mosum]{mosum.cpts}
#' @param bootstrap flag indicating whether bootstrap replicates of estimated changepoints
#' should be computed
#' @param N_bootstrap number of bootstrap replicates to be generated (iff bootstrap)
#' @param alpha_CI numeric value in (0,1), such that the (1-alpha_CI)-confidence bootstrap intervals are computed (iff bootstrap)
#' @param ... further arguments to be passed to the \code{mosum} calls
#' for respective bandwidths; see \link[mosum]{mosum}
#' @return S3 \code{multiscale.cpts} object, which contains the following fields:
#'    \item{x}{the numeric input vector provided}
#'    \item{cpts}{estimated changepoints after merging}
#'    \item{cpts.info}{data frame containing information about estimated changepoints}
#'    \item{pooled.cpts}{set of changepoint candidates that have been considered during the algorithm}
#'    \item{G,eta}{input parameter}
#'    \item{threshold,threshold.function}{input parameter}
#'    \item{bootstrap}{input parameter}
#'    \item{cpts_bootstrap}{bootstrap replicates and CIs, object of class \link[mosum]{cpts.bootstrap} (iff bootstrap)}
#' @details See Algorithm 1 in the first referenced paper for a comprehensive
#' description of the procedure and further details.
#' @references A. Meier, C. Kirch and H. Cho.
#' \emph{mosum: A Package for Moving Sums in Change Point Analysis.}
#' Unpublished manuscript, 2018+.
#' @references H. Cho and C. Kirch.
#' \emph{Multiple change-point detection via multiscale MOSUM procedure with localized pruning.}
#' Unpublished manuscript, 2018+.
#' @references M. Messer et al. 
#' \emph{A multiple filter test for the detection of rate changes in renewal processes with varying variance.}
#' Annals of Applied Statistics, 8(4):2027-2067, 2014.
#' @examples # many-bandwidth example
#' x <- piecewiseStationary_timeSeries(lengths=c(50,50,200,300,300), 
#'   means=c(0,1,2,3,2.3), sds=rep(1,5))
#' G <- (5:20)*5
#' mcpts <- multiscale.bottomUp.cpts(x, G, alpha=0.1)
#' mcpts$pooled.cpts
#' mcpts$cpts
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export
multiscale.bottomUp.cpts <- function(x, G, 
                                     threshold = c("critical.value", "custom")[1],
                                     alpha=0.05, threshold.function = NULL, eta=0.8, 
                                     bootstrap = F, N_bootstrap, alpha_CI=0.05, ...) {
  n <- length(x)

  if (class(G) == "integer" || class(G) == "numeric") {
    grid <- multiscale.grid(G, method="cartesian")
  } else {
    stopifnot(class(G) == "multiscale.grid")
    if (any(apply(G$grid,1,diff) != 0)) {
      stop("Expecting a grid of symmetric bandwidths")
    }
    grid <- G
  }
  args <- list(...)
  
  abs.bandwidth <- all(grid$grid>=1)
  if (abs.bandwidth) {
    GRID_THRESH <- max(20, 0.05*n)
  } else {
    GRID_THRESH <- 0.05
  }
  
  if (min(grid$grid) < GRID_THRESH) {
    print("Warning: Smallest bandwidth in grid is relatively small (in comparison to n)")
    # warning("Smallest bandwidth in grid is relatively small (in comparison to n)")
  }
  if (threshold != "critical.value" && threshold != "custom") {
    stop("threshold must be either 'critical.value' or 'custom'")
  }

  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(eta >= 0)

  # Retreive change point candidates from all bandwidths.
  cpts.complete <- numeric(0)
  bandwidths.complete <- integer(0)
  pValues.complete <- numeric(0)
  m <- array(list(),nrow(grid$grid))
  
  for (i in seq_len(nrow(grid$grid))) {
    G <- grid$grid[[i,1]]
    m[[i]] <- mosum(x, G, ...)
    if (threshold == "critical.value") {
      cpts <- (mosum.cpts(x, G, ... , threshold="critical.value", alpha=alpha, 
                          criterion="eta", eta=eta))$cpts
    } else {
      threshold_val <- threshold.function(G)
      cpts <- (mosum.cpts(x, G, ... , threshold="custom", threshold.custom=threshold_val,
                          criterion="eta", eta=eta))$cpts
    }
    cpts.complete <- c(cpts.complete, cpts)
    if (!abs.bandwidth) {
      G <- floor(n*G)
    }
    cv <- mosum.criticalValue(length(x),G,G,alpha)
    bandwidths.complete <- c(bandwidths.complete, rep(G, length(cpts)))
    pValues.complete <- c(pValues.complete, mosum.pValue(m[[i]]$stat[cpts],length(x),G))
  }

  # Merge candidates.
  points <- numeric(0)
  bandwidths <- numeric(0)
  pValues <- numeric(0)
  cptsInOrder <- seq_len(length(cpts.complete))
  for (i in cptsInOrder) {
    p <- cpts.complete[[i]]
    G <- bandwidths.complete[[i]]
    pVal <- pValues.complete[[i]]
    if (suppressWarnings(min(abs(p-points))) >= eta*G) { # Note: min(empty_list) = Inf
      points <- c(points, p)
      bandwidths <- c(bandwidths, G)
      pValues <- c(pValues, pVal)
    }
  }
  cpts.merged <- data.frame(point = points, 
                            bandwidth_l = bandwidths, 
                            bandwidth_r = bandwidths,
                            bandwidth = 2*bandwidths,
                            pValue = pValues)
  cpts <- as.matrix(cpts.merged[order(cpts.merged$point),])
  
  ret <- structure(list(x=x,
                        G=G,
                        threshold=threshold,
                        alpha=alpha,
                        threshold.function=threshold.function,
                        rule="bottomUp_merging",
                        criterion="eta",
                        eta=eta,
                        epsilon=NULL,
                        cpts=as.numeric(cpts[,1]), 
                        cpts.info=cpts,
                        pooled.cpts=unique(cpts.complete), 
                        bootstrap=F), # note 
                   class="multiscale.cpts")
  
  if (bootstrap) {
    ret$cpts_bootstrap <- cpts.bootstrap(ret, N_bootstrap, alpha_CI)
    ret$bootstrap <- T
  }
  ret
}
