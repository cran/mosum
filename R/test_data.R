#' Test data with piecewise constant mean
#' 
#' Generate piecewise stationary time series with change-points in the mean.
#' @param model a string indicating from which model a realisation is to be generated;
#' possible values are 'blocks', 'fms', 'mix', 'stairs10', 'teeth10'
#' (for the referenced model signals) or 'custom' (for user specification 
#' in terms of \code{lengths}, \code{means} and \code{sds})
#' @param lengths use iff \code{model='custom'}; an integer vector for the lengths of the piecewise stationary segments
#' @param means use iff \code{model='custom'}; a numeric vector for the means of the piecewise stationary segments
#' @param sds use iff \code{model='custom'}; a numeric vector for the deviation scaling of the piecewise stationary segments.
#' The values are multiplied to the outcome of \code{rand.gen}, coinciding with the standard
#' deviation in case of standard normal innovations (\code{rand.gen=rnorm})
#' @param rand.gen optional; a function to generate the noise/innovations
#' @param seed optional; if a seed value is provided (\code{!is.null(seed)}), 
#' then \code{set.seed(seed)} is called beforehand)
#' @param ... further arguments to be parsed to \code{rand.gen}
#' @return numeric vector containing a realisation of the time series, given as signal+noise
#' @details See Appendix B in the reference for details about the model time series.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @importFrom stats rnorm
#' @export
testData <- function(model=c('custom', 'blocks', 'fms', 'mix', 'stairs10', 'teeth10')[1],
                                           lengths=NULL, means=NULL, sds=NULL,
                                           rand.gen=rnorm, seed=NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  signal <- testSignal(model=model, lengths=lengths, means=means, sds=sds)
  n <- length(signal$mu_t)
  ts <- signal$mu_t + rand.gen(n, ...)*signal$sigma_t
  return(ts)
}

#' Piecewise constant test signal
#' 
#' Produce vectors of mean and dispersion values for generating piecewise stationary time series.
#' @param model a string indicating which model is to be used; 
#' possible values are 'blocks', 'fms', 'mix', 'stairs10', 'teeth10'
#' (for the referenced model signals) or 'custom' (for user specification 
#' in terms of \code{lengths}, \code{means} and \code{sds});
#' @param lengths use iff \code{model='custom'}; an integer vector for the lengths of the piecewise stationary segments
#' @param means use iff \code{model='custom'}; a numeric vector for the means of the piecewise stationary segments
#' @param sds use iff \code{model='custom'}; a numeric vector for the deviation scaling of the piecewise stationary segments.
#' The values are multiplied to the outcome of \code{rand.gen}, coinciding with the standard
#' deviation in case of standard normal innovations (\code{rand.gen=rnorm})
#' @return a list containing the following entries:
#'   \item{mu_t}{mean vector of piecewise stationary model time series}
#'   \item{sigma_t}{deviation scaling vector of piecewise stationary model time series}
#' @details See Appendix B in the reference for details about the model time series.
#' @references P. Fryzlewicz (2014)
#' Wild Binary Segmentation for Multiple Change-Point Detection.
#' \emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
#' @export
testSignal <- function(model=c('custom', 'blocks', 'fms', 'mix', 'stairs10', 'teeth10')[1],
                                       lengths=NULL, means=NULL, sds=NULL) {
  if (model=='blocks') {
    res <- modelSignal.blocks()
  } else if (model=='fms') {
    res <- modelSignal.fms()
  } else if (model=='mix') {
    res <- modelSignal.mix()
  } else if (model=='stairs10') {
    res <- modelSignal.stairs10()
  } else if (model=='teeth10') {
    res <- modelSignal.teeth10()
  } else if (model=='custom') {
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
    stop('Unkonwn model string')
  }
  return(res)
}
