#' MOSUM procedure for change-point estimation
#' 
#' Computes the MOSUM detector, detects (multiple) change-points and estimates their locations.
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G an integer value for the moving sum bandwidth;
#' alternatively a number between \code{0} and \code{0.5} describing the moving sum bandwidth
#' relative to \code{length(x)}
#' @param G.right if \code{!is.na(G.right)}, the asymmetric bandwidth (G, G.right) will be used
#' @param var.est.method how the variance is estimated;
#' possible values are
#' \itemize{
#'    \item{\code{'custom'}}{a vector of \code{length(x)} is to be parsed by the user; use \code{var.custom} in this case to to so}
#'    \item{\code{'mosum'}}{both-sided MOSUM variance estimator}
#'    \item{\code{'mosum.min'}}{minimum of the sample variance estimates from the left and right summation windows}
#' }
#' @param var.custom a numeric vector (of the same length as \code{x}) containing
#' local estimates of the variance or long run variance; use iff \code{var.est.method=custom}
#' @param boundary.extension a logical value indicating whether the boundary
#' values should be filled-up with CUSUM values  
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, it is chosen from the asymptotic distribution at the given significance level \code{alpha}.
#' Alternatively it is possible to parse a user-defined numerical value with \code{threshold.custom}
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold='critical.value'}
#' @param threshold.custom a numeric value greater than 0 for the threshold of significance;
#' use iff \code{threshold='custom'}
#' @param criterion how to decide whether each point \code{p} at which MOSUM statistic 
#' exceeds the threshold is a change-point; possible values are
#' \itemize{
#'    \item{\code{'eta'}}{there is no bigger exceeding in an \code{eta*G} environment of \code{p}}
#'    \item{\code{'epsilon'}}{\code{p} is the maximum of its local exceeding environment, which has at least size \code{epsilon*G}}
#' }
#' @param eta a positive numeric value for the minimal mutual distance of 
#' changes, relative to moving sum bandwidth (iff \code{criterion='eta'})
#' @param epsilon a numeric value in (0,1] for the minimal size of exceeding
#' environments, relative to moving sum bandwidth (iff \code{criterion='epsilon'})
#' @param do.confint flag indicating whether to compute the confidence intervals for change-points
#' @param level use iff \code{do.confint=TRUE}; a numeric value (\code{0 <= level <= 1}) with which
#' \code{100(1-level)\%} confidence interval is generated
#' @param N_reps use iff \code{do.confint=TRUE}; number of bootstrap replicates to be generated
#' @return S3 \code{mosum.cpts} object, which contains the following fields:
#'    \item{x}{input data}
#'    \item{G.left,G.right}{left and right bandwidths}
#'    \item{var.est.method,var.custom,boundary.extension}{input parameters}
#'    \item{stat}{a series of MOSUM statistic values; the first \code{G} and last \code{G.right} values are \code{NA} iff \code{boundary.extension=F}}
#'    \item{rollsums}{a series of MOSUM detector values; equals \code{stat*sqrt(var.estimation)}}
#'    \item{var.estimation}{the local variance estimated according to \code{var.est.method}}
#'    \item{threshold,alpha,threshold.custom}{input parameters}
#'    \item{critical.value}{critical value of the corresponding MOSUM test}
#'    \item{criterion,eta,epsilon}{input parameters}
#'    \item{cpts}{a vector containing the estimated change-point locations}
#'    \item{do.confint}{input parameter}
#'    \item{ci}{object of class \code{cpts.ci} containing confidence intervals for change-points iff \code{do.confint=TRUE}}
#' @references A. Meier, C. Kirch and H. Cho (2018+)
#' mosum: A Package for Moving Sums in Change Point Analysis. \emph{Unpublished manuscript}.
#' @references B. Eichinger and C. Kirch (2018)
#' A MOSUM procedure for the estimation of multiple random change-points.
#' \emph{Bernoulli}, Volume 24, Number 1, pp. 526-564.
#' @examples 
#' x <- testData(lengths=rep(100, 3), means=c(0, 5, -2), sds=rep(1, 3))
#' m <- mosum(x, G=40)
#' plot(m)
#' summary(m)
#' @export
mosum <- function(x, G, G.right=NA, var.est.method='mosum', var.custom=NULL, boundary.extension=T,
                      threshold=c('critical.value', 'custom')[1], alpha=0.05, threshold.custom=NULL, 
                      criterion=c('eta', 'epsilon')[1], eta=0.4, epsilon=0.2,
                      do.confint=F, level=0.05, N_reps=1000) {
  
  # Consistency checks on input
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(criterion=='epsilon' || criterion=='eta')
  stopifnot(criterion!='epsilon' || epsilon >= 0)
  stopifnot(criterion!='eta' || eta >= 0)
  stopifnot(!do.confint || N_reps>0)

  n <- length(x)
  m <- mosum.stat(x, G, G.right, var.est.method, var.custom, boundary.extension)
  
  G.left <- m$G.left
  G.right <- m$G.right
  G.min <- min(G.right, G.left)
  G.max <- max(G.right, G.left)
  K <- G.min / G.max
  changePoints <- numeric(0)
  
  if (threshold == 'critical.value') {
    threshold_val <- mosumCriticalValue(n, G.left, G.right, alpha)
  } else if (threshold == 'custom') {
    threshold_val <- threshold.custom
  } else {
    stop('threshold must be either \'critical.value\' or \'custom\'')
  }
  
  # get exceeding TRUE/FALSE vector
  exceedings <- (m$stat > threshold_val)
  
  # adjust, in case of no boundary CUSUM extension
  if (!m$boundary.extension) {
    exceedings[n-G.right+1] <- F
  }
  
  if (criterion=='epsilon') {
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
          # print(paste0('Found change-point at the right border (G=(', G.left, ',', G.right, ')).'))
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
          # print(paste0('Found change-point at the left border (G=(', G.left, ',', G.right, ')).'))
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
  } else { # (criterion=='eta')
    localMaxima <- (c((diff.default(m$stat) < 0), NA) & c(NA, diff.default(m$stat) > 0))
    # adjust, in case of no boundary CUSUM extension
    if (!m$boundary.extension) {
      localMaxima[n-G.right] <- T
    }
    p.candidates <- which(exceedings & localMaxima)
    changePoints <- eta_criterion_help(p.candidates, m$stat, eta, G.left, G.right)
  }
  ret <- structure(
    list(x=x,
        G.left=G.left,
        G.right=G.right,
        var.est.method=m$var.est.method,
        var.custom=m$var.custom, 
        boundary.extension=m$boundary.extension,
        stat=m$stat,
        rollsums=m$rollsums,
        var.estimation=m$var.estimation,        
        threshold=threshold,
        alpha=alpha,
        threshold.custom=threshold.custom,
        critical.value=mosumCriticalValue(n, G.left, G.right, alpha),
        criterion=criterion,
        eta=eta,        
        epsilon=epsilon,
        cpts=changePoints,
        do.confint=F,
        ci=NA), # note
    class='mosum.cpts')
  if (do.confint) {
    ret$ci <- confint.mosum.cpts(ret, level=level, N_reps=N_reps)
    ret$do.confint <- T 
  } 
  ret
}

#' Plotting MOSUM change-points
#' 
#' Plotting method for objects of class \code{mosum.cpts}
#' @method plot mosum.cpts
#' @param x a \code{mosum.cpts} object
#' @param critical.value.col a specification for the color of the
#' critical value, see \link[graphics]{par}
#' @param cpts.col a specification for the color of the
#' change-points, see \link[graphics]{par}
#' @param xlab graphical parameter
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' and \link[graphics]{abline}
#' @details The input time series is plotted along with vertical lines indicating
#' the estimated change-point locations and a horizontal line indicating
#' the critical value used in the MOSUM test.
#' @export
plot.mosum.cpts <- function(x, critical.value.col='blue', cpts.col='red', xlab='Time', ...) {
  if (class(x$x)=='ts') {
    x_plot <- as.numeric(time(x$x))
  } else {
    x_plot <- seq_len(length(x$x))
  }
  plot(x=x_plot, y=x$stat, type='l', xlab=xlab, ...)
  abline(h=x$critical.value, col=critical.value.col)
  for (p in x$cpts) {
    abline(v=x_plot[p], col='red', ...)
  }
}

#' Summary of change-points estimated by MOSUM procedure
#' 
#' Summary method for objects of class \code{mosum.cpts}
#' @method summary mosum.cpts
#' @param object a \code{mosum.cpts} object
#' @param ... not in use
#' @details Provide information about each estimated change-point, 
#' including the bandwidths used for its estimation, associated p-value and (scaled) jump size.
#' @export
summary.mosum.cpts <- function(object, ...) { 
  n <- length(object$x)
  if(length(object$cpts) > 0){
    ans <- cbind(object$cpts, object$G.left, object$G.right, 
                 signif(mosum.pValue(object$stat[object$cpts], n, object$G.left, object$G.right), 3), 
                 round(object$stat[object$cpts]*sqrt(object$G.left+object$G.right)/sqrt(object$G.left*object$G.right), 3))
    colnames(ans) <- c('cpts', 'G.left', 'G.right', 'p.value', 'jump')
    if(object$do.confint) ans <- cbind(ans, object$ci$CI[, -1, drop=FALSE])
  } else{
    ans <- matrix(c(NA, object$G.left, object$G.right, NA, 0), nrow=1)
    colnames(ans) <- c('cpts', 'G.left', 'G.right', 'p.value', 'jump')
  }    
#  cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated at alpha = ', object$alpha, ' according to ', object$criterion, '-criterion', sep=''))
  if(object$criterion=='eta') cat(paste('\n with eta = ', object$eta, sep=''))
  if(object$criterion=='epsilon') cat(paste('\n with epsilon = ', object$epsilon, sep=''))
  cat(paste(' and ', object$var.est.method, ' variance estimate:', sep=''))
  cat('\n')
  cat('\n')
  print(ans, print.gap = 3)
}

#' Change-points estimated by MOSUM procedure
#' 
#' Print method for objects of class \code{mosum.cpts}
#' @method print mosum.cpts
#' @param x a \code{mosum.cpts} object
#' @param ... not in use
#' @export
print.mosum.cpts <- function(x, ...) {
  #cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated with bandwidths (', x$G.left, ', ', x$G.right, ')', 
            ' at alpha = ', x$alpha, sep='')) 
  cat(paste('\n according to ', x$criterion, '-criterion', sep=''))
  if(x$criterion=='eta') cat(paste(' with eta = ', x$eta, sep=''))
  if(x$criterion=='epsilon') cat(paste(' with epsilon = ', x$epsilon, sep=''))
  cat(paste(' and ', x$var.est.method, ' variance estimate:', sep=''))
  cat('\n')
  cat('\n')
  cat('  ')
  cat(x$cpts)
}