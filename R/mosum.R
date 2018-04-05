#' MOSUM statistic
#' 
#' Computes the statistical values for the MOSUM test for changes in the mean.
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G an integer value for the length of the moving sum window; 
#' alternatively a number between 0 and 0.5 describing the moving sum bandwidth
#' relative to n (the length of x);
#' @param G.right iff \code{!is.na(G.right)}, the asymmetric bandwidth (G,G.right)
#' will be used;
#' @param var.est.method how the variance should be estimated;
#' possible values are
#' \itemize{
#'    \item{"custom" if a vector is to be parsed by the user;
#'     use \code{var.custom} in this case to to so}
#'    \item{"mosum" for the both-sided MOSUM variance estimator}
#'    \item{"mosum.min" for the minimum of the sample variance estimates from the left and right summation windows}
#' }
#' @param var.custom a numeric vector (of the same length as \code{x}) containing
#' local estimates of the variance or long run variance; use iff \code{var.est.method=custom}
#' @param boundary.extension a logical value indicating whether the boundary
#' values should be filled-up with CUSUM values
#' @return S3 \code{mosum} object, which contains the following fields:
#'    \item{x}{the numeric input vector provided}
#'    \item{G,symmetric,G.right}{input parameters}
#'    \item{var.est.method,var.custom,boundary.extension}{input parameters}
#'    \item{stat}{statistical result values; note that the first and last G values are \code{NA} iff \code{boundary.extension=F}}
#'    \item{rollsums}{equals \code{stat*sqrt(var.estimation)}}
#'    \item{var.estimation}{the local variance estimation according 
#'    to \code{var.est.method}}
#' @details This class only contains the values for the MOSUM statistic.
#' For statistical evaluation and changepoint extraction, use
#' \link[mosum]{mosum.test} and \link[mosum]{mosum.cpts}.
#' See also \link[mosum]{multiscale.cpts}.
#' @references A. Meier, C. Kirch and H. Cho.
#' \emph{mosum: A Package for Moving Sums in Change Point Analysis.}
#' Unpublished manuscript, 2018+.
#' 
#' @references B. Eichinger and C. Kirch.
#' \emph{A MOSUM procedure for the estimation of multiple random change points.}
#' Bernoulli, Volume 24, Number 1 (2018), 526-564.
#' @examples # Statistical values and variance est. for independent data
#' x <- piecewiseStationary_timeSeries(lengths=rep(100,3), means=c(0,5,-2), sds=rep(1,3))
#' m <- mosum(x, G=40)
#' par(mfcol=c(3,1))
#' plot.ts(x)
#' plot(m)
#' plot.ts(m$var.estimation)
#' par(mfcol=c(1,1))
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export 
mosum <- function(x, G, G.right=NA, var.est.method="mosum", 
                  var.custom=NULL, boundary.extension=T) {
  n <- length(x)
  
  symmetric <- is.na(G.right)
  
  abs.bandwidth <- (G>=1)
  if (!abs.bandwidth) {
    G <- floor(n * G)
    if (!symmetric) G.right <- floor(n * G.right)
  }
  
  # Consistency checks on input
  stopifnot(NCOL(x) == 1) 
  stopifnot(class(x)=="ts" || class(x)=="numeric")
  stopifnot(G > 0 && G < n)
  stopifnot(symmetric || !is.na(G.right))
  stopifnot(symmetric || (G.right > 0 && G.right < n))
  
  G.left <- G
  G.right.inputParam <- G.right
  if (symmetric) {
    G.right <- G
  }
  G.min <- min(G.right, G.left)
  G.max <- max(G.right, G.left)
  K <- G.min / G.max
  
  # Calculate value of statistic.
  sums.left <- rolling_sum(x, G.left) # zoo::rollsum(x, k=G.left, fill = NA, align="left")
  if (G.left == G.right) {
    sums.right <- sums.left
  } else {
    sums.right <- rolling_sum(x, G.right) # zoo::rollsum(x, k=G.right, fill = NA, align="left")
  }
  unscaledStatistic <- c(rep(NA,G.left-1), (G.min/G.right*sums.right[(G.left+1):n] - G.min/G.left*sums.left[1:(n - G.left)]), NA) / (sqrt((K+1)*G.min))
  
  # Calculate variance estimation.
  if (!is.null(var.custom) && var.est.method != "custom") {
    stop("Please use var.est.method=custom when parsing var.custom.")
  }
  if (var.est.method == "custom") {
    if (is.null(var.custom)) {
      stop("Expecting var.custom to be not NULL for var.est.method=custom")
    }
    if (length(var.custom) != n) {
      stop("Expecting var.custom to be of length n=length(x)")
    }
    var <- var.custom
  } else if (var.est.method == "global") {
    # Note: This is Deprecated
    var <- rep((sum(x^2) - (sum(x)^2)/n)/n,n)
  } else { # MOSUM-based variance estimators
    summedSquares.left <- rolling_sum(x^2, G.left) # zoo::rollsum(x^2, k=G.left, fill=NA, align="left")
    squaredSums.left <- sums.left^2
    var.tmp.left <- summedSquares.left[1:(n-G.left+1)] - 1/G.left*(squaredSums.left[1:(n-G.left+1)])
    var.left <- c(rep(NA,G.left-1), var.tmp.left) / G.left
    if (G.left == G.right) {
      summedSquares.right <- summedSquares.left
      squaredSums.right <- squaredSums.left
      var.tmp.right <- var.tmp.left
    } else {
      summedSquares.right <- rolling_sum(x^2, G.right) # zoo::rollsum(x^2, k=G.right, fill=NA, align="left")
      squaredSums.right <- sums.right^2
      var.tmp.right <- summedSquares.right[1:(n-G.right+1)] - 1/G.right*(squaredSums.right[1:(n-G.right+1)])
    }
    var.right <- c(var.tmp.right[2:(n-G.right+1)], rep(NA,G.right)) / G.right
    if (var.est.method == "mosum") {
      var <- (var.left + var.right) / 2
    } else if (var.est.method == "mosum.left") {
      # Note: This is Deprecated
      var <- var.left
    } else if (var.est.method == "mosum.right") {
      # Note: This is Deprecated
      var <- var.right
    } else if (var.est.method == "mosum.min") {
      var <- pmin(var.left, var.right)
    } else {
      stop("unknown variance estimation method")
    }
  }
  var.estimation <- var

  # CUSUM extension to boundary
  if (boundary.extension) {
    if (n > 2*G.left) {
      weights.left <- sqrt( (G.left+G.right) / (1:G.left) / ((G.left+G.right-1):(G.right))) 
      unscaledStatistic[1:G.left] <- cumsum(mean(x[1:(G.left+G.right)])-x[1:G.left]) * 
        weights.left
      var.estimation[1:G.left] <- var.estimation[G.left]
    }
    if (n > 2*G.right) {
      weights.right <- sqrt( (G.left+G.right) / (1:G.right) / ((G.left+G.right-1):(G.left))) 
      x.rev <- rev(x[(n-2*G.right+1):n])
      unscaledStatistic[(n-G.right+1):n] <- rev( cumsum(mean(x.rev[1:(2*G.right)])-x.rev[1:G.right])  * 
                                    weights.right)
      var.estimation[(n-G.right+1):n] <- var.estimation[n-G.right]
    }
  }
  res <- abs(unscaledStatistic) / sqrt(var.estimation)
  structure(
    list(x=x,
         G=G, 
         symmetric=symmetric,
         G.right=G.right.inputParam,
         var.est.method=var.est.method,
         var.custom=var.custom, 
         stat=res,
         rollsums=unscaledStatistic,
         var.estimation=var.estimation,
         boundary.extension=boundary.extension),
    class="mosum")
}

#' Plotting MOSUM objects
#'
#' Plotting method for objects of class "mosum".
#' @method plot mosum
#' @param x a \code{mosum} object
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}
#' @param critical.value.col a specification for the color of the
#' critical value, see \link[graphics]{par}
#' @param xlab graphical parameter
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' @importFrom stats plot.ts time
#' @export
plot.mosum <- function(x, alpha=0.05, critical.value.col="blue", 
                       xlab="Time", ...) {
  if (class(x$x)=="ts") {
    x_plot <- as.numeric(time(x$x))
  } else {
    x_plot <- seq_len(length(x$x))
  }
  plot(x=x_plot, y=x$stat, type='l', xlab=xlab, ...)
  G.left <- x$G
  if (x$symmetric) {
    G.right <- x$G
  } else {
    G.right <- x$G.right
  }
  abline(h=mosum.criticalValue(n=length(x$x), G.left=G.left,
                               G.right=G.right, alpha=alpha),
         col=critical.value.col, ...)
}
