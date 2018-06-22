#' Multiscale MOSUM algorithm with localised pruning
#' 
#' Multiscale (asymmetric) MOSUM procedure with localised pruning using Schwarz criterion (SC).
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G a vector of bandwidths, given as either integers or numbers between 
#' \code{0} and \code{0.5} describing the moving sum bandwidths relative to \code{length(x)}
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, it is chosen from the asymptotic distribution at the given significance level \code{alpha}.
#' Alternatively, it is possible to parse a user-defined function with \code{threshold.function}
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}. Use iff \code{threshold='critical.value'}
#' @param threshold.function function object of form \code{function(G_l, G_r, length(x), alpha)}, to compute a
#' threshold of significance for different bandwidths (G_l, G_r); use iff \code{threshold='custom'}
#' @param criterion how to decide whether an exceeding point is a change-point; to be parsed to \link[mosum]{mosum}
#' @param eta,epsilon see \link[mosum]{mosum}
#' @param rule string for the choice of sorting criterion in merging step. 
#' Possible values are: 
#' \itemize{
#' \item{\code{'pval'}}{lowest p-value}
#' \item{\code{'jump'}}{largest (rescaled) jump size}
#' }
#' @param penalty string indicating which kind of penalty term is used; possible values are:
#' \itemize{
#' \item{\code{'log'}}{use \code{penalty = log(length(x))^pen.exp}}
#' \item{\code{'polynomial'}}{use \code{penalty = length(x)^pen.exp}}
#' }
#' @param pen.exp exponent for penalty term (see \code{penalty});
#' @param do.confint flag indicating whether confidence intervals for change-points should be computed
#' @param level use iff \code{do.confint=TRUE}; a numeric value (\code{0 <= level <= 1}) with which
#' \code{100(1-level)\%} confidence interval is generated
#' @param N_reps use iff \code{do.confint=TRUE}; number of bootstrap replicates to be generated
#' @param ... further arguments to be parsed to \link[mosum]{mosum} calls
#' @return S3 \code{multiscale.cpts} object, which contains the following fields:
#'    \item{x}{input data}
#'    \item{cpts}{estimated change-points}
#'    \item{cpts.info}{data frame containing information about estimated change-points}
#'    \item{sc}{SC value of the estimated change-point set}
#'    \item{pooled.cpts}{set of change-point candidates that have been considered by the algorithm}
#'    \item{G}{input parameter}
#'    \item{threshold,alpha,threshold.function}{input parameters}
#'    \item{criterion,eta,epsilon}{input parameters}
#'    \item{rule,penalty,pen.exp}{input parameters}
#'    \item{do.confint}{input parameter}
#'    \item{ci}{object of class \code{cpts.ci} containing confidence intervals for change-points iff \code{do.confint=TRUE}}
#' @details See Algorithm 2 in the first referenced paper for a comprehensive
#' description of the procedure and further details.
#' @references A. Meier, C. Kirch and H. Cho (2018+)
#' mosum: A Package for Moving Sums in Change Point Analysis. \emph{Unpublished manuscript}.
#' @references H. Cho and C. Kirch (2018+)
#' Multiple change-point detection via multiscale MOSUM procedure with localised pruning. \emph{Unpublished manuscript}.
#' @examples 
#' x <- testData(model='mix')
#' mlp <- multiscale.localPrune(x, G=c(8, 15, 30, 70), do.confint=TRUE)
#' print(mlp)
#' summary(mlp)
#' par(mfcol=c(2, 1), mar=c(2, 4, 2, 2))
#' plot.ts(x)
#' abline(v=mlp$cpts, col=2)
#' plot(mlp)
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export
multiscale.localPrune <- function(x, G=bandwidths.default(length(x)),
                            threshold=c('critical.value', 'custom')[1], alpha=.05, threshold.function = NULL,
                            criterion=c('eta', 'epsilon')[1], eta=0.4, epsilon=0.2,
                            rule=c('pval', 'jump')[1], penalty=c('log', 'polynomial')[1], pen.exp=1.01,
                            do.confint=F, level=0.05, N_reps=1000, ...) {
  
  THRESH_MANUAL_MERGING <- 24
  
  n <- length(x)
  
  if (class(G) == "integer" || class(G) == "numeric") {
    grid <- multiscale.grid(G)
  } else if(class(G) == 'multiscale.grid'){
    grid <- G
  } else stop('Expecting a vector of numbers')
  abs.bandwidth <- all(grid$grid>=1)
  
  stopifnot(is.element(rule, c('pval', 'jump', 'lr', 'rl')))
  stopifnot(is.element(criterion, c('eta', 'epsilon')))
  stopifnot((criterion=='eta' & eta <= 1 & eta > 0) || (criterion=='epsilon' & epsilon <= 1 & epsilon > 0))
  stopifnot(!do.confint || N_reps>0)
  
  if (penalty == 'log') {
    log_penalty <- T
  } else {
    if (penalty != 'polynomial') {
      stop('penalty has to set to log or polynomial')
    }
    log_penalty <- F
  }
  
  if (threshold != 'critical.value' && threshold != 'custom') {
    stop('threshold must be either \'critical.value\' or \'custom\'')
  }
  
  all.cpts <- matrix(NA, ncol=6, nrow=0)
  for (i in seq_len(nrow(grid$grid))) {
    G1 <- grid$grid[[i,1]]
    G2 <- grid$grid[[i,2]]
    if (threshold == 'critical.value') {
      m <- mosum(x, G=G1, G.right=G2, ..., 
                        threshold='critical.value', alpha=alpha, 
                        criterion=criterion, eta=eta, epsilon=epsilon)
    } else{
      threshold_val <- threshold.function(G1, G2, n, alpha)
      m <- mosum(x, G=G1, G.right=G2, ..., 
                         threshold='custom', threshold.custom=threshold_val, alpha=alpha,
                         criterion=criterion, eta=eta, epsilon=epsilon)
    }

    if(length(m$cpts)>0){
      if (!abs.bandwidth) {
        G1 <- floor(G1*n)
        G2 <- floor(G2*n)
      }
      all.cpts <- rbind(all.cpts, 
                        cbind(m$cpts, G1, G2, G1+G2,
                              mosum.pValue(m$stat[m$cpts], n, G1, G2), m$stat[m$cpts]*sqrt(G1+G2)/sqrt(G1*G2)))
    }
  }

  all.cpts <- all.cpts[sort(all.cpts[, 1], decreasing=FALSE, index.return=TRUE)$ix,,drop=F]
  all.cpts <- dup.merge(all.cpts, rule) # if there are duplicates, only select one according to 'rule'
  ac <- nrow(all.cpts)
  all.cpts <- cbind(all.cpts, matrix(0, nrow=ac, ncol=2))
  # the following for loop defines the 'neighbours' among the change-points according to their respective G-environment
  for(j in seq_len(ac)){
    k <- all.cpts[j, 1]
    if(j==1) all.cpts[j, 7] <- 0 else {
      ind <- ((j-1):1)[k-all.cpts[(j-1):1, 1] < all.cpts[j, 2]]
      if(length(ind) > 0 & sum(k-all.cpts[ind, 1] < all.cpts[ind, 3]) > 0) all.cpts[j, 7] <- max(which(k-all.cpts[ind, 1] < all.cpts[ind, 3]))
    }
    if(j==ac) all.cpts[j, 8] <- 0 else {
      ind <- ((j+1):ac)[all.cpts[(j+1):ac, 1]-k < all.cpts[j, 3]]
      if(length(ind) > 0 & sum(all.cpts[ind, 1]-k < all.cpts[ind, 2]) > 0) all.cpts[j, 8] <- max(which(all.cpts[ind, 1]-k < all.cpts[ind, 2]))
    }
  }
  
  cand_used <- rep(F, ac)
  all.unique.cpts <- c(0, all.cpts[, 1], n)
  auc <- length(all.unique.cpts)-2
  sums <- matrix(0, nrow=auc+1, ncol=4) # calculated for efficient computation of rss
  for(j in 1:(auc+1)){
    sums[j, 1:2] <- c(all.unique.cpts[j]+1, all.unique.cpts[j+1])
    sums[j, 3] <- sum(x[sums[j, 1]:sums[j, 2]])
    sums[j, 4] <- sum(x[sums[j, 1]:sums[j, 2]]^2)
  }
  min.cost <- sum(sums[, 4]-sums[, 3]^2/(sums[, 2]-sums[, 1]+1)) # min rss with all the candidates
  
  if(rule=='pval') rule.seq <- sort(all.cpts[, 5], decreasing=FALSE, index.return=TRUE)$ix
  if(rule=='jump') rule.seq <- sort(all.cpts[, 6], decreasing=!FALSE, index.return=TRUE)$ix
  if(rule=='lr') rule.seq <- pool
  if(rule=='rl') rule.seq <- rev(pool)
  
  current <- pool <- seq_len(ac); est.cpts.ind <- est.cpts <- integer(0)
  # current = C, pool = P, est.cpts.ind = B (index)
  while(length(pool)>0){
    # step 1
    j <- rule.seq[1]; adj <- 0
    # step 2
    left <- max(est.cpts.ind[est.cpts.ind<j]+1, j-all.cpts[j, 7])
    right <- min(est.cpts.ind[est.cpts.ind>j]-1, j+all.cpts[j, 8])
    cand.ind <- (left:right)[is.element(left:right, pool)]
    cand <- all.cpts[cand.ind, 1] # = D
    
    # step 3
    if(sum(current < left) > 0) li <- max(current[current < left]) else li <- 0
    if(sum(current > right) > 0) ri <- min(current[current > right]) else ri <- ac+1  
    
    ind_middl_tmp <- sums[(li+1):(ri-1), 2]
    ind_middl_tmp <- ind_middl_tmp[which(!cand_used[(li+1):(ri-1)])]
    ind_tmp <- c(sums[li+1,1]-1, ind_middl_tmp, sums[ri, 2])
    sub.sums <- extract_sub(ind_tmp, x)
    
    doExhaustiveSearch <- T
    # Too many candidates to do exhaustive search?
    if (length(cand) > THRESH_MANUAL_MERGING) {
      
      # Count neighbourhood size of remaining candidates
      cand_size <- rep(NA, length(rule.seq))
      cand_size[1] <- length(cand)
      for (i_tmp in seq(from=2, length.out=length(rule.seq)-1)) {
        jj <- rule.seq[i_tmp]
        left_jj <- max(est.cpts.ind[est.cpts.ind<jj]+1, jj-all.cpts[jj, 7])
        right_jj <- min(est.cpts.ind[est.cpts.ind>jj]-1, jj+all.cpts[jj, 8])
        cand.ind_jj <- (left_jj:right_jj)[is.element(left_jj:right_jj, pool)]
        cand_jj <- all.cpts[cand.ind_jj, 1] # = D
        cand_size[i_tmp] <- length(cand_jj)
      }
      
      if (any(cand_size <= THRESH_MANUAL_MERGING)) {
        # Proceed with next candidate, for which exhaustive search IS possible
        ind_star <- min(which(cand_size <= THRESH_MANUAL_MERGING))
        rule_tmp <- rule.seq[ind_star]; rule.seq[ind_star] <- rule.seq[1]; rule.seq[1] <- rule_tmp
        doExhaustiveSearch <- F
      } else {
        # No more exhaustive search possible at all
        # --> Do manual merging, until exhaustive search becomes possible
        while(length(cand) > THRESH_MANUAL_MERGING) {
          warn_msg <- paste0('Warning: ', length(cand), ' conflicting candidates, thinning manually')
          print(warn_msg) #; warning(warn_msg)
          k <- cand[which.min(diff(cand))]
          l <- which(sub.sums[, 2]==k)
          a <- sub.sums[l, ]; b <- sub.sums[l+1, ]
          # as change-points are merged, the minimum rss in the local environment needs to be updated
          adj <- adj + (a[2]-a[1]+1)*(b[2]-b[1]+1)/(b[2]-a[1]+1)*(a[3]/(a[2]-a[1]+1)-b[3]/(b[2]-b[1]+1))^2
          sub.sums[l+1, 1] <- a[1]; sub.sums[l+1, 3:4] <- sub.sums[l+1, 3:4]+a[3:4]
          sub.sums <- sub.sums[-l, ]
          cand <- setdiff(cand, k)
          k.ind <- which(all.cpts[, 1]==k)
          cand.ind <- setdiff(cand.ind, k.ind); pool <- setdiff(pool, k.ind); rule.seq <- setdiff(rule.seq, k.ind)
          cand_used[k.ind] <- T
        }
      }
    }
    
    if (doExhaustiveSearch) {
      # step 4
      # performs exhaustive search (Algorithm 2)
      out <- exhaust_bic(cand=cand, sub_sums=sub.sums, 
                         strength=pen.exp, log_penalty=log_penalty, 
                         n=n, auc=length(current), min_cost=min.cost)
      est.cpts <- c(est.cpts, out$est_cpts)
      est.cpts.ind <- c(est.cpts.ind, which(is.element(all.cpts[, 1], out$est_cpts)))
      
      # steps 5, 6
      # only those change-point candidates that have been judged in the 'fair' environment are removed
      rm.set <- c(j, est.cpts.ind,
                  cand.ind[(all.cpts[cand.ind, 1]-all.cpts[cand.ind, 2] >= sub.sums[1, 1] &
                              all.cpts[cand.ind, 1]+all.cpts[cand.ind, 3] <= sub.sums[nrow(sub.sums), 2]) |
                             (all.cpts[cand.ind, 1] >= suppressWarnings(min(out$est_cpts)) &
                                all.cpts[cand.ind, 1] <= suppressWarnings(max(out$est_cpts)))])
      pool <- setdiff(pool, rm.set)
      cand_used[rm.set] <- T
      rule.seq <- setdiff(rule.seq, rm.set)
      current <- which(is.element(setdiff(sums[, 2], n), c(all.cpts[pool, 1], est.cpts)))
      current_cands <- is.element(cand, all.cpts[current, 1])
      ind_star <- get_comb_ind(current_cands)
      min.cost <- min.cost + adj - out$bic[nrow(out$bic), 1] + out$bic[ind_star+1, 1]
    }
  }
  
  est.cpts <- sort(as.vector(est.cpts)); est.cpts.ind <- sort(as.vector(est.cpts.ind))
  est.cpts.info <- data.frame(cpts = all.cpts[est.cpts.ind, 1], 
                          G.left =  all.cpts[est.cpts.ind, 2], 
                          G.right =  all.cpts[est.cpts.ind, 3],
                          p.value = all.cpts[est.cpts.ind, 5],
                          jump = all.cpts[est.cpts.ind, 6])
  if (log_penalty) {
    penalty_term <- length(est.cpts)*log(n)^pen.exp
  } else {
    penalty_term <- length(est.cpts)*n^pen.exp
  }
  final.bic <- n/2*log(min.cost/n) + penalty_term
  if(!abs.bandwidth) G <- floor(n * sort(unique(c(grid$grid))))
  
  ret <- structure(list(x=x,
                        cpts=est.cpts, 
                        cpts.info=est.cpts.info,
                        sc=final.bic, 
                        pooled.cpts=all.cpts[, 1], 
                        G=G,
                        alpha=alpha,
                        threshold=threshold,
                        threshold.function=threshold.function,
                        criterion=criterion,
                        eta=eta,
                        epsilon=epsilon,
                        rule=rule,
                        penalty=penalty,
                        pen.exp=pen.exp,
                        do.confint=F,
                        ci=NA), # note 
                   class='multiscale.cpts')
  if (do.confint) {
    ret$ci <- confint.multiscale.cpts(ret, level=level, N_reps=N_reps)
    ret$do.confint <- T
  }
  ret
}

#' Remove duplicated from all.cpts data frame:
#' In case one change being added multiple times, choose the one
#' with smallest p-value (resp. highest jump)
#' @keywords internal
dup.merge <- function(all.cpts, rule='jump') {
  all.unique.cpts <- unique(all.cpts[, 1, drop=F])
  out <- matrix(NA, nrow=0, ncol=ncol(all.cpts))
  for(k in all.unique.cpts){
    ind <- which(all.cpts[, 1]==k)
    ind.min <- ind[all.cpts[ind, 4]==min(all.cpts[ind, 4])]
    if(length(ind.min) > 1 & rule=='pval') ind.min <- ind.min[which.min(all.cpts[ind.min, 5])] 
    if(length(ind.min) > 1 & rule=='jump') ind.min <- ind.min[which.max(all.cpts[ind.min, 6])] 
    out <- rbind(out, all.cpts[ind.min,])
  }
  out
}

#' Plotting multiscale.cpts change-points
#' 
#' Plotting method for objects of class 'multiscale.cpts'.
#' @method plot multiscale.cpts
#' @param x a \code{multiscale.cpts} object
#' @param shaded string indicating which to display as shaded area on the y-axis around change-point estimates.
#' Poissble values are \code{'bandwidth'} for the respective detection intervals and \code{'CI'} for (bootstrapped) confidence intervals
#' @param level,N_reps argument to be parsed to \link[mosum]{confint.multiscale.cpts}; use iff \code{shaded='CI'}
#' @param CI string indicating whether pointwise (\code{CI = 'pw'}) or uniform (\code{CI = 'unif'}) confidence intervals
#' shall be plotted; use iff \code{shaded='CI'}
#' @param ... not in use
#' @details The input time series is plotted along with the estimated change-points
#' and their detection intervals (if \code{shaded='bandwidth'}) or bootstrap confidence intervals
#' of the corresponding change-points (if \code{shaded='CI'}). The y-axis represents \code{1-p.value}
#' for each estimated change-point derived from the asymptotic distribution of MOSUM statistics.
#' @importFrom grDevices rainbow
#' @importFrom graphics abline axis plot points rect segments
#' @export
plot.multiscale.cpts <- function(x, shaded=c('CI', 'bandwidth')[1], 
                                 level=0.05, N_reps=10000, 
                                 CI = c('pw', 'unif')[1], ...) {
  if (shaded=='bandwidth') {
    main <- 'Change-point estimates and detection intervals'
  } else if (shaded=='CI') {
    if(CI=='pw') main <- 'Change-point estimates and pointwise confidence intervals'
    if(CI=='unif') main <- 'Change-point estimates and uniform confidence intervals'
    if(x$do.confint) b <- x$ci else b <- confint.multiscale.cpts(x, level=level, N_reps=N_reps)
  } else {
    stop('shaded argument has to be either \'CI\' or \'bandwidth\'.')
  }
  n <- length(x$x)
  cpts <- x$cpts.info
  cols <- rainbow(nrow(cpts),alpha=0.2)
  cols2 <- rainbow(nrow(cpts),alpha=1)
  # cols3 <- rainbow(nrow(cpts),alpha=0.2)
  newOrder <- 1:nrow(cpts) #sample(seq_len(nrow(cpts)), nrow(cpts))
  cols <- cols[newOrder]
  cols2 <- cols2[newOrder]
  
  y.min <- 1-1.1*max(cpts$p.value) # pvalue 
  plot(0,type='n',xlim=c(0,n), ylim=c(y.min,1),bty='n',axes=F,xlab='time',ylab='1-p.value',main=main)
  axis(side = 1)
  axis(side = 2)
  xx <- cpts$cpts # location
  if (shaded=='bandwidth') {
    xx_l <- xx-cpts$G.left+1
    xx_r <- xx+cpts$G.right
  }
  if (shaded == 'CI') {
    if(CI == 'pw'){
      xx_l <- b$CI[,2]
      xx_r <- b$CI[,3]
    } else if(CI == 'unif'){
      xx_l <- b$CI[,4]
      xx_r <- b$CI[,5]
    }
  }
  yy <- 1-cpts$p.value # pvalue
  points(x=xx, y=yy)
  rect(xleft=xx_l, xright=xx_r, ybottom=0, ytop=yy, col=cols, lty=0)
  # if (shaded=='CI' && include_unif_CI) {
  #   # also add uniform ones
  #   xx_ll <- b$CI[,4]
  #   xx_rr <- b$CI[,5]
  #   rect(xleft=xx_ll, xright=xx_rr, ybottom=0, ytop=yy, col=cols3, lty=0)
  # }
  segments(x0=xx,y0=0,y1=yy,col=cols2)
}

#' Summary of change-points estimated by multiscale MOSUM procedure
#' 
#' Summary method for objects of class \code{multiscale.cpts}
#' @method summary multiscale.cpts
#' @param object a \code{multiscale.cpts} object
#' @param ... not in use
#' @details Provide information about each estimated change-point, 
#' including the bandwidths used for its detection, associated p-value and (scaled) jump size;
#' if \code{object$do.confint=TRUE}, end points of the pointwise and uniform confidence intervals
#' are also provided.
#' @export
summary.multiscale.cpts <- function(object, ...) { 
  n <- length(object$x)
  ans <- object$cpts.info
  ans$p.value <- signif(ans$p.value, 3)
  ans$jump <- round(ans$jump, 3)
  if(object$do.confint) ans <- cbind(ans, object$ci$CI[, -1, drop=FALSE])
  #cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated at alpha = ', object$alpha, ' according to ', object$criterion, '-criterion', sep=''))
  if(object$criterion=='eta') cat(paste('\n with eta = ', object$eta, sep=''))
  if(object$criterion=='epsilon') cat(paste('\n with epsilon = ', object$epsilon, ':', sep=''))
  cat('\n')
  cat('\n')
  print(ans, print.gap = 3)
}

#' Change-points estimated by multiscale MOSUM procedure
#' 
#' Print method for objects of class \code{multiscale.cpts}
#' @method print multiscale.cpts
#' @param x a \code{multiscale.cpts} object
#' @param ... not in use
#' @export
print.multiscale.cpts <- function(x, ...) {
  #cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated with bandwidths\n'))
  cat('  ')
  cat(x$G)
  cat(paste('\nat alpha = ', x$alpha, ' according to ', x$criterion, '-criterion', sep=''))
  if(x$criterion=='eta') cat(paste(' with eta = ', x$eta, ':', sep=''))
  if(x$criterion=='epsilon') cat(paste(' with epsilon = ', x$epsilon, ':', sep=''))
  cat('\n')
  cat('\n')
  cat('  ')
  cat(x$cpts)
}