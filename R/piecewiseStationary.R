#' Piecewise stationary model time series
#' 
#' Create a piecewise stationary model time series (signal+noise)
#' @param model A string indicating which model to be returned; 
#' possible values are "blocks", "fms", "mix", "stairs10", "teeth10"
#' (for the referenced model signals) or "custom" (for user specification 
#' in terms of \code{lengths}, \code{means} and \code{sds})
#' @param lengths (use iff model="custom"): integer vector for the lengths of the piecewise stationary parts
#' @param means (use iff model="custom"): numeric vector for the means of the piecewise stationary parts
#' @param sds (use iff model="custom"): numeric vector for the deviation scaling of the piecewise stationary parts.
#' The values are multiplied with the outcome of \code{rand.gen}. Note that it corresponds to the standard
#' deviation in case of standard normal innovations (\code{rand.gen=rnorm})
#' @param rand.gen optional: a function to generate the noise/innovations
#' @param seed optional: a seed value (if !is.null(seed), then set.seed(seed) is called beforehand)
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return numeric vector containing a realization of the time series, given as signal+noise
#' @importFrom stats rnorm
#' @export
piecewiseStationary_timeSeries <- function(model=c("custom", "blocks", "fms", "mix", 
                                                   "stairs10", "teeth10")[1],
                                           lengths=NULL, means=NULL, sds=NULL,
                                           rand.gen=rnorm, 
                                           seed=NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  signal <- piecewiseStationary_signal(model=model, lengths=lengths,
                                       means=means, sds=sds)
  n <- length(signal$mu_t)
  ts <- signal$mu_t + rand.gen(n, ...)*signal$sigma_t
  return(ts)
}

#' Piecewise stationary model time series signal
#' 
#' Get mean and sd vector of piecewise stationary time series
#' @param model A string indicating which model to be returned; 
#' possible values are "blocks", "fms", "mix", "stairs10", "teeth10"
#' (for the referenced model signals) or "custom" (for user specification 
#' in terms of \code{lengths}, \code{means} and \code{sds});
#' @param lengths (use iff model="custom"): integer vector for the lengths of the piecewise stationary parts;
#' @param means (use iff model="custom"): numeric vector for the means of the piecewise stationary parts;
#' @param sds (use iff model="custom"): numeric vector for the deviation scaling of the piecewise stationary parts.
#' The values are multiplied with the outcome of \code{rand.gen}. Note that it corresponds to the standard
#' deviation in case of standard normal innovations (\code{rand.gen=rnorm});
#' @return a list containing the following entries:
#'   \item{mu_t}{mean vector of piecewise stationary model time series}
#'   \item{sigma_t}{deviation scaling vector of piecewise stationary model time series}
#' @details See Appendix B in the first reference for details about the model time series.
#' @references P. Fryzlewicz.
#' \emph{Wild Binary Segmentation for Multiple Change-Point Detection.}
#' The Annals of Statistics, Vol. 42, No. 6, 2243-2281, 2014.
#' @references A. Meier, C. Kirch and H. Cho.
#' \emph{mosum: A Package for Moving Sums in Change Point Analysis.}
#' Unpublished manuscript, 2018+.
#' @export
piecewiseStationary_signal <- function(model=c("custom", "blocks", "fms", "mix", 
                                               "stairs10", "teeth10")[1],
                                       lengths=NULL, means=NULL, sds=NULL) {
  if (model=="blocks") {
    res <- modelSignal.blocks()
  } else if (model=="fms") {
    res <- modelSignal.fms()
  } else if (model=="mix") {
    res <- modelSignal.mix()
  } else if (model=="stairs10") {
    res <- modelSignal.stairs10()
  } else if (model=="teeth10") {
    res <- modelSignal.teeth10()
  } else if (model=="custom") {
    stopifnot(!is.null(lengths))
    stopifnot(!is.null(means))
    stopifnot(!is.null(sds))
    stopifnot(length(lengths)==length(means))
    stopifnot(length(lengths)==length(sds))
    mu_t <- numeric(0)
    sigma_t <- numeric(0)
    for(i in 1:length(lengths)) {
      mu_t <- c(mu_t, rep(means[i], lengths[i]))
      sigma_t <- c(sigma_t, rep(sds[i], lengths[i]))
    }
    res <- list(mu_t=mu_t,
                sigma_t=sigma_t)
  } else {
    stop("Unkonwn model string")
  }
  return(res)
}
