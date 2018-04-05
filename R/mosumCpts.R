#' MOSUM changepoints
#' 
#' Extract the changepoints found by MOSUM method.
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G an integer value for the moving sum bandwidth
#' @param G.right if \code{!is.na(G.right)}, the asymmetric bandwidth (G,G.right)
#' will be used;
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, \code{mosum.criticalValue} is used with significance level \code{alpha}
#'  Alternatively it is possible to parse a user-defined numerical value with 
#'  \code{threshold.custom}.
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff threshold="critical.value"
#' @param threshold.custom a numeric value greater than 0 for the threshold of significance;
#' use iff threshold="custom"
#' @param criterion how to decide whether an exceeding point p
#' is a change point,
#' possible values are
#' \itemize{
#'    \item{"epsilon" p is the maximum of its local exceeding environment, 
#'    which has at least size epsilon*G}
#'    \item{"eta" there is no bigger exceeding in an eta*G environment of p}
#' }
#' @param epsilon a numeric value in (0,1] for the minimal size of exceeding
#' environments, relative to moving sum bandwidth (iff criterion=="epsilon")
#' @param eta a numeric value > 0 for the minimal mutual distance of 
#' changes, relative to moving sum bandwidth (iff criterion=="eta")
#' @param bootstrap flag indicating whether bootstrap replicates of estimated changepoints
#' should be computed
#' @param N_bootstrap number of bootstrap replicates to be generated (iff bootstrap)
#' @param alpha_CI numeric value in (0,1), such that the (1-alpha_CI)-confidence bootstrap intervals are computed (iff bootstrap)
#' @param ... further arguments to be parsed to call of \link[mosum]{mosum}
#' @return S3 \code{mosum.cpts} object, which contains the following fields:
#'    \item{m}{object of class \link[mosum]{mosum} representing the underlying MOSUM statistic}
#'    \item{threshold,alpha,threshold.custom}{input parameter}
#'    \item{epsilon,eta}{input parameter}
#'    \item{critical.value}{the critical value of the corresponding MOSUM test}
#'    \item{cpts}{list of computed changepoints}
#'    \item{bootstrap}{input parameter}
#'    \item{cpts_bootstrap}{bootstrap replicates and CIs, object of class \link[mosum]{cpts.bootstrap} (iff bootstrap)}
#' @examples x <- piecewiseStationary_timeSeries(lengths=rep(100,3), means=c(0,5,-2), sds=rep(1,3))
#' m.cpts <- mosum.cpts(x, G=40)
#' par(mfcol=c(2,1))
#' plot.ts(x)
#' plot(m.cpts)
#' @export
mosum.cpts <- function(x, G, G.right=NA, threshold = c("critical.value", "custom")[1],
                       alpha=0.05, threshold.custom=NULL, criterion="epsilon", 
                       epsilon=0.2, eta=1, 
                       bootstrap=F, N_bootstrap, alpha_CI=0.05, ...) {
  n <- length(x)
  
  m <- mosum(x, G, G.right, ...)
  
  # Consistency checks on input
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(criterion=="epsilon" || criterion=="eta")
  stopifnot(criterion!="epsilon" || epsilon >= 0)
  stopifnot(criterion!="eta" || eta >= 0)
  stopifnot(!bootstrap || N_bootstrap>0)
  
  G.left <- m$G
  if (m$symmetric) {
    G.right <- m$G
  } else {
    G.right <- m$G.right
  }
  G.min <- min(G.right, G.left)
  G.max <- max(G.right, G.left)
  K <- G.min / G.max
  changePoints <- numeric(0)
  
  if (threshold == "critical.value") {
    threshold_val <- mosum.criticalValue(n, G.left, G.right, alpha)
  } else if (threshold == "custom") {
    threshold_val <- threshold.custom
  } else {
    stop("Threshold must be set to either 'critical.value' or 'custom'")
  }
  
  # get exceeding TRUE/FALSE vector
  exceedings <- (m$stat > threshold_val)
  
  # adjust, in case of no boundary CUSUM extension
  if (!m$boundary.extension) {
    exceedings[n-G.right+1] <- F
  }
  
  if (criterion=="epsilon") {
    # get number of subsequent exceedings
    exceedingsCount <- (exceedings) * unlist(lapply(rle(exceedings)$lengths, seq_len))
    # get exceeding-intervals of fitting length
    minIntervalSize <- max(1, (G.min+G.max) / 2 * epsilon)
    intervalEndPoints <- which(diff(exceedingsCount) <= -minIntervalSize)
    intervalBeginPoints <- intervalEndPoints - exceedingsCount[intervalEndPoints] + 1
    if (!m$boundary.extension) {
      # manually adjust right border
      if (exceedings[n-G.right] && !((n-G.right) %in% intervalEndPoints)) {
        lastBeginPoint <- n - G.right - exceedingsCount[n-G.right] + 1
        stopifnot(exceedings[seq(lastBeginPoint,n-G.right)])
        stopifnot(!(lastBeginPoint %in% intervalBeginPoints))
        highestStatPoint <- which.max(m$stat[seq(lastBeginPoint,n-G.right)]) + lastBeginPoint - 1
        if (highestStatPoint-lastBeginPoint >= minIntervalSize/2) {
          # print(paste0("Found changepoint at the right border (G=(", G.left, ",", G.right, "))."))
          intervalEndPoints <- c(intervalEndPoints, n-G.right)
          intervalBeginPoints <- c(intervalBeginPoints, lastBeginPoint)
        }
      }
      # manually adjust left border
      if (exceedings[G.left] && !(G.left %in% intervalBeginPoints)) {
        firstEndPoint <- which(diff(exceedingsCount) < 0)[1]
        stopifnot(exceedings[seq(G.left,firstEndPoint)])
        stopifnot(!(firstEndPoint %in% intervalEndPoints))
        highestStatPoint <- which.max(m$stat[seq(G.left,firstEndPoint)]) + G.left - 1
        if (firstEndPoint - highestStatPoint >= minIntervalSize/2) {
          # print(paste0("Found changepoint at the left border (G=(", G.left, ",", G.right, "))."))
          intervalEndPoints <- c(firstEndPoint, intervalEndPoints)
          intervalBeginPoints <- c(G.left, intervalBeginPoints)
        }
      }
    }
    numChangePoints <- length(intervalBeginPoints)
    if (numChangePoints > 0) {
      for (i in 1:numChangePoints) {
        changePoint <- intervalBeginPoints[i] + which.max(m$stat[(intervalBeginPoints[i]):(intervalEndPoints[i])]) - 1
        changePoints <- c(changePoints, changePoint)
      }
    }
  } else { # (criterion=="eta")
    localMaxima <- (c((diff.default(m$stat) < 0), NA) & c(NA, diff.default(m$stat) > 0))
    # adjust, in case of no boundary CUSUM extension
    if (!m$boundary.extension) {
      localMaxima[n-G.right] <- T
    }
    p.candidates <- which(exceedings & localMaxima)
    changePoints <- eta_criterion_help(p.candidates, m$stat, eta, G.left, G.right)
  }
  ret <- structure(
    list(m=m,
         alpha=alpha,
         epsilon=epsilon,
         eta=eta,
         threshold=threshold,
         threshold.custom=threshold.custom,
         critical.value=mosum.criticalValue(n, G.left, G.right, alpha),
         cpts=changePoints,
         bootstrap=F), # note
    class="mosum.cpts")
  if (bootstrap) {
    ret$cpts_bootstrap <- cpts.bootstrap(ret, N_bootstrap, alpha_CI)
    ret$bootstrap <- T 
  } 
  ret
}

#' Plotting MOSUM changepoints
#' 
#' Plotting method for objects of class "mosum.cpts".
#' @method plot mosum.cpts
#' @param x a \code{mosum.cpts} object
#' @param critical.value.col a specification for the color of the
#' critical value, see \link[graphics]{par}
#' @param cpts.col a specification for the color of the
#' changepoints, see \link[graphics]{par}
#' @param xlab graphical parameter
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' and \link[graphics]{abline}
#' @export
plot.mosum.cpts <- function(x, critical.value.col="blue", cpts.col="red", xlab="Time", ...) {
  if (class(x$m$x)=="ts") {
    x_plot <- as.numeric(time(x$m$x))
  } else {
    x_plot <- seq_len(length(x$m$x))
  }
  plot(x=x_plot, y=x$m$stat, type='l')
  abline(h=x$critical.value,col=critical.value.col, xlab=xlab, ...)
  for (p in x$cpts) {
    abline(v=x_plot[p],col="red", ...)
  }
}