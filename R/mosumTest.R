#' Help function: asymptotic scaling.
#' @keywords internal
mosum.asymptoticA <- function(x, K) {
  return(sqrt(2*log(x)))
}
#' Help function: asymptotic shift
#' @keywords internal
mosum.asymptoticB <- function(x, K) {
  return(2*log(x) + 0.5*log(log(x)) + log((K^2+K+1)/(K+1)) - 0.5*log(pi))
}
#' MOSUM asymptotic critical value
#' 
#' Computes the asymptotic critical value for the MOSUM test.
#' @param n an integer value for the length of the input data
#' @param G.left,G.right integer values for the left moving sum bandwidth (G.left,G.right)
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}
#' @return a numeric value for the asymptotic critical value for the MOSUM test
#' @export
mosum.criticalValue <- function(n, G.left, G.right, alpha) {
  G.min <- min(G.left, G.right)
  G.max <- max(G.left, G.right)
  K <- G.min / G.max
  return((mosum.asymptoticB(n/G.min,K) - log(log(1/sqrt(1-alpha))))/mosum.asymptoticA(n/G.min,K))
}

#' MOSUM asymptotic p-value
#' 
#' Computes the asymptotic p-value for the MOSUM test.
#' @param z a numeric value for the observation
#' @param n an integer value for the length of the input data
#' @param G.left,G.right integer values for the left moving sum bandwidth (G.left,G.right)
#' @return a numeric value for the asymptotic p-value for the asymmetric MOSUM test
#' @export
mosum.pValue <- function(z, n, G.left, G.right=G.left) {
  G.min <- min(G.left, G.right)
  G.max <- max(G.left, G.right)
  K <- G.min / G.max
  return(1-exp(-2*exp(mosum.asymptoticB(n/G.min,K) - mosum.asymptoticA(n/G.min,K)*z)))
}

#' MOSUM test
#' 
#' Computes the result of the MOSUM test for changes in the mean.
#' @param m a \code{mosum} object
#' @return S3 \code{mosum.test} object, which contains the following fields:
#'    \item{m}{the input \code{mosum} object}
#'    \item{statistic}{the value of the MOSUM statistic}
#'    \item{p.value}{the p-value of the test}
#'    \item{estimate}{the changepoint estimation}
#' @examples x <- piecewiseStationary_timeSeries(lengths=rep(100,3), means=c(0,2,-1), sds=rep(1,3))
#' m <- mosum(x, G=40)
#' mTest <- mosum.test(m)
#' mTest$estimate
#' mTest$p.value
#' @export 
mosum.test <- function(m) {
  stopifnot(class(m)=="mosum")
  n <- length(m$stat)
  estimate <- which.max(m$stat)
  statistic <- m$stat[estimate]
  G.left <- m$G
  if (!m$symmetric)
    G.right <- m$G.right
  else
    G.right <- G.left
  p.value <- mosum.pValue(statistic, n, G.left, G.right)
  structure(
    list(m=m,
         statistic=statistic,
         p.value=p.value,
         estimate=estimate),
    class="mosum.test")
}