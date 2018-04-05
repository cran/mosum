#' Multiscale MOSUM algorithm (sBIC)
#' 
#' Multiscale MOSUM procedure with sBIC merging
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G an asymmetric bandwidth grid;
#' either an object of type \code{multiscale.grid} or an integer vector
#' of one-sided bandwidths to be passed as parameter to \link[mosum]{multiscale.grid}
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, \code{mosum.criticalValue} is used for each bandwidth (with significance
#' level \code{alpha}). Alternatively, it is possible to parse a user-defined function 
#' with \code{threshold.function};
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}. Use iff \code{threshold="critical.value"}
#' @param threshold.function function object of form \code{function(G_l,G_r)}, to compute a
#' threshold of significance for different bandwidths (G_l,G_r); use iff \code{threshold="custom"}
#' @param rule string for the choice of sorting criterion in merging step. 
#' Possible values are: 
#' \itemize{
#' \item{\code{pval} = lowest p-value}, 
#' \item{\code{peak} = largest (rescaled) peak.}
#' }
#' @param criterion Criterion how changepoint candidates are generated,
#' to be parsed to \link[mosum]{mosum.cpts}
#' @param eta,epsilon see \link[mosum]{mosum.cpts}
#' @param penalty string indicating which kind of penalty term is used; possible values are:
#' \itemize{
#' \item{\code{log}: in this case, penalty(n)=log(n)^pen.exp}
#' \item{\code{polynomial}: in this case, penalty=n^pen.exp}
#' }
#' @param pen.exp exponent for penalty term (see \code{penalty});
#' @param bootstrap flag indicating whether bootstrap replicates of estimated changepoints
#' should be computed
#' @param N_bootstrap number of bootstrap replicates to be generated (iff bootstrap)
#' @param alpha_CI numeric value in (0,1), such that the (1-alpha_CI)-confidence bootstrap intervals are computed (iff bootstrap)
#' @param ... further arguments to be parsed to \link[mosum]{mosum} calls
#' @return S3 \code{multiscale.cpts} object, which contains the following fields:
#'    \item{x}{the numeric input vector provided}
#'    \item{cpts}{estimated changepoints after merging}
#'    \item{cpts.info}{data frame containing information about estimated changepoints}
#'    \item{bic}{sBIC value of estimated changepoint set}
#'    \item{pooled.cpts}{set of changepoint candidates that have been considered during the algorithm}
#'    \item{G,rule,criterion}{input parameter}
#'    \item{threshold,threshold.function}{input parameter}
#'    \item{eta,epsilon,penalty,pen.exp}{input parameter}
#'    \item{bootstrap}{input parameter}
#'    \item{cpts_bootstrap}{bootstrap replicates and CIs, object of class \link[mosum]{cpts.bootstrap} (iff bootstrap)}
#' @details See Algorithm 2 in the first referenced paper for a comprehensive
#' description of the procedure and further details.
#' @references A. Meier, C. Kirch and H. Cho.
#' \emph{mosum: A Package for Moving Sums in Change Point Analysis.}
#' Unpublished manuscript, 2018+.
#' @references H. Cho and C. Kirch.
#' \emph{Multiple change-point detection via multiscale MOSUM procedure with localized pruning.}
#' Unpublished manuscript, 2018+.
#' @examples 
#' x <- piecewiseStationary_timeSeries(model="mix")
#' mcpts <- multiscale.cpts(x, G=c(8,15,30,70))
#' print(mcpts$cpts)
#' print(mcpts$cpts.info)
#' par(mfcol=c(2,1), mar=c(2,4,2,2))
#' plot.ts(x)
#' abline(v=mcpts$cpts, col=2)
#' plot(mcpts)
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export
multiscale.cpts <- function(x, G=bandwidths.default(length(x)),
                            threshold = c("critical.value", "custom")[1],
                            alpha=.05, threshold.function = NULL,
                            rule=c('pval', 'peak')[1], 
                            criterion=c('eta', 'epsilon')[1], eta=.4, epsilon=.2, 
                            penalty=c("log", "polynomial")[1], pen.exp=1.1,
                            bootstrap=F, N_bootstrap, 
                            alpha_CI=0.05, ...) {
  
  THRESH_MANUAL_MERGING <- 24
  
  n <- length(x)
  
  if (class(G) == "integer" || class(G) == "numeric") {
    grid <- multiscale.grid(G)
  } else {
    stopifnot(class(G) == "multiscale.grid")
    grid <- G
  }
  args <- list(...)
  abs.bandwidth <- all(grid$grid>=1)
  
  stopifnot(is.element(rule, c('pval', 'peak', 'lr', 'rl')))
  stopifnot(is.element(criterion, c('eta', 'epsilon')))
  stopifnot((criterion=='eta' & eta <= 1 & eta > 0) || (criterion=='epsilon' & epsilon <= 1 & epsilon > 0))
  stopifnot(!bootstrap || N_bootstrap>0)
  
  if (penalty == "log") {
    log_penalty <- T
  } else {
    if (penalty != "polynomial") {
      stop("penalty has to set to log or polynomial")
    }
    log_penalty <- F
  }
  
  if (threshold != "critical.value" && threshold != "custom") {
    stop("threshold must be either 'critical.value' or 'custom'")
  }
  
  all.cpts <- matrix(NA, ncol=6, nrow=0)
  all.points <- 1:n
  for (i in seq_len(nrow(grid$grid))) {
    G1 <- grid$grid[[i,1]]
    G2 <- grid$grid[[i,2]]
    if (threshold == "critical.value") {
      mcpt <- mosum.cpts(x, G=G1, G.right=G2, ..., 
                         threshold="critical.value", alpha=alpha, 
                         criterion=criterion, 
                         eta=eta, epsilon=epsilon)
    } else {
      threshold_val <- threshold.function(G1, G2)
      mcpt <- mosum.cpts(x, G=G1, G.right=G2, ..., 
                         threshold="custom", threshold.custom=threshold_val, 
                         criterion=criterion, 
                         eta=eta, epsilon=epsilon)
    }

    if(length(mcpt$cpts)>0){
      if (!abs.bandwidth) {
        G1 <- floor(G1*n)
        G2 <- floor(G2*n)
      }
      all.cpts <- rbind(all.cpts, cbind(mcpt$cpts, G1, G2, G1+G2, 
                                        mosum.pValue(mcpt$m$stat[mcpt$cpts], n, G1, G2), mcpt$m$stat[mcpt$cpts]*sqrt(G1+G2)/sqrt(G1*G2)))
      for(t in mcpt$cpts) all.points <- setdiff(all.points, (t-G1+1):(t+G2))
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
  if(rule=='peak') rule.seq <- sort(all.cpts[, 6], decreasing=!FALSE, index.return=TRUE)$ix
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
          warn_msg <- paste0("Warning: ", length(cand), " conflicting candidates, thinning manually")
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
      out <- exhaust_bic(cand=cand, 
                         sub_sums=sub.sums, 
                         strength=pen.exp, 
                         log_penalty=log_penalty, 
                         n=n, 
                         auc=length(current), 
                         min_cost=min.cost)

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
  est.cpts.info <- all.cpts[est.cpts.ind, 1:6, drop=FALSE]
  colnames(est.cpts.info) <- c('point', 'bandwidth_l', 'bandwidth_r', 'bandwidth', 'pValue', 'peak')
  
  if (log_penalty) {
    penalty_term <- length(est.cpts)*log(n)^pen.exp
  } else {
    penalty_term <- length(est.cpts)*n^pen.exp
  }
  final.bic <- n/2*log(min.cost/n) + penalty_term
  
  ret <- structure(list(x=x,
                        G=G,
                        threshold=threshold,
                        alpha=alpha,
                        threshold.function=threshold.function,
                        rule=rule,
                        criterion=criterion,
                        eta=eta,
                        epsilon=epsilon,
                        cpts=est.cpts, 
                        bic=final.bic, 
                        cpts.info=est.cpts.info,
                        pooled.cpts=all.cpts[, 1], 
                        penalty=penalty,
                        pen.exp=pen.exp,
                        bootstrap=F), # note 
                   class="multiscale.cpts")
  if (bootstrap) {
    ret$cpts_bootstrap <- cpts.bootstrap(ret, N_bootstrap, alpha_CI)
    ret$bootstrap <- T
  }
  ret
}

#' Remove duplicated from all.cpts data frame:
#' In case one change being added multiple times, choose the one
#' with smallest p-Value (resp. highest peak)
#' @keywords internal
dup.merge <- function(all.cpts, rule='peak') {
  all.unique.cpts <- unique(all.cpts[, 1, drop=F])
  out <- matrix(NA, nrow=0, ncol=ncol(all.cpts))
  for(k in all.unique.cpts){
    ind <- which(all.cpts[, 1]==k)
    ind.min <- ind[all.cpts[ind, 4]==min(all.cpts[ind, 4])]
    if(length(ind.min) > 1 & rule=='pval') ind.min <- ind.min[which.min(all.cpts[ind.min, 5])] 
    if(length(ind.min) > 1 & rule=='peak') ind.min <- ind.min[which.max(all.cpts[ind.min, 6])] 
    out <- rbind(out, all.cpts[ind.min,])
  }
  out
}

#' Plotting multiscale.cpts changepoints
#' 
#' Plotting method for objects of class "multiscale.cpts".
#' @method plot multiscale.cpts
#' @param x an \code{multiscale.cpts} object
#' @param shaded string indicating what to display as shaded area on the y-axis around change point estimates.
#' Poissble values are \code{"bandwidth"} for the respective detection bandwidth and \code{"CI"} for (bootstrapped) confidence intervals;
#' @param N_reps,alpha_CI argument to be parsed to \code{cpts.bootstrap}. Use iff \code{shaded="CI"}.
#' @param include_uniform_CI shall uniform confidence intervals also be plotted? Use iff \code{shaded="CI"}.
#' @param ... not in use
#' @importFrom grDevices rainbow
#' @importFrom graphics abline axis plot points rect segments
#' @export
plot.multiscale.cpts <- function(x, shaded=c("CI", "bandwidth")[1], 
                                 N_reps=10000, alpha_CI=0.05, 
                                 include_uniform_CI=F, ...) {
  if (shaded=="bandwidth") {
    main <- "Change point estimates and detection environments"
  } else if (shaded=="CI") {
    main <- "Change point estimates and confidence intervals"
    b <- cpts.bootstrap(x, N_reps, alpha_CI)
  } else {
    stop("shaded argument has to be either 'CI' or 'bandwidth'.")
  }
  n <- length(x$x)
  cpts <- x$cpts.info
  cols <- rainbow(nrow(cpts),alpha=0.2)
  cols2 <- rainbow(nrow(cpts),alpha=1)
  cols3 <- rainbow(nrow(cpts),alpha=0.2)
  newOrder <- 1:nrow(cpts) #sample(seq_len(nrow(cpts)), nrow(cpts))
  cols <- cols[newOrder]
  cols2 <- cols2[newOrder]
  
  y.min <- 1-1.1*max(cpts[,5]) # pvalue 
  plot(0,type='n',xlim=c(0,n), ylim=c(y.min,1),bty="n",axes=F,xlab="Time",ylab="1-pval",
       main=main)
  axis(side = 1)
  axis(side = 2)
  xx <- cpts[,1] # location
  if (shaded=="bandwidth") {
    xx_l <- xx-cpts[,2]+1
    xx_r <- xx+cpts[,3]-1
  }
  if (shaded == "CI") {
    xx_l <- b$CI[,2]
    xx_r <- b$CI[,3]
  }
  yy <- 1-cpts[,5] # pvalue
  points(x=xx, y=yy)
  rect(xleft=xx_l, xright=xx_r, ybottom=0, ytop=yy, col=cols, lty=0)
  if (shaded=="CI" && include_uniform_CI) {
    # also add uniform ones
    xx_ll <- b$CI[,4]
    xx_rr <- b$CI[,5]
    rect(xleft=xx_ll, xright=xx_rr, ybottom=0, ytop=yy, col=cols3, lty=0)
  }
  segments(x0=xx,y0=0,y1=yy,col=cols2)
}
