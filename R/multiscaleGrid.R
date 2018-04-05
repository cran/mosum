#' Default choice of multiple bandwidth set
#' 
#' Create bandwidths according to a default formula of the sample size
#' @param n integer representing the sample size
#' @param d.min integer for the minimal mutual distance of change points that can be expected
#' @param G.min integer for the minimal allowed bandwidth
#' @param G.max integer for the maximal allowed bandwidth
#' @details Returns an integer vector of bandwidths (G_1,...,G_m), 
#' with G_0 = G_1 = max(G.min, 2/3*d.min), and G_{j+1} = G_{j-1} + G_j (for j=1,...,m-1)
#' and m such that G_m <= G.max.
#' @return an integer vector of bandwidths
#' @export
bandwidths.default <- function(n, d.min=10, G.min=8, G.max=3*sqrt(n)) {
  G_0 <- G_1 <- max(G.min, round(2/3*d.min))
  G <- c(G_0, G_1)
  j <- 3
  G_j <- G[j-2] + G[j-1]
  while (G_j <= G.max) {
    G[j] <- G_j
    j <- j+1
    G_j <- G[j-2] + G[j-1]
  }
  G[-1]
}

#' Multiscale bandwidth grids
#' 
#' Create asymmetric bandwidth grids to be used with \code{multiscale.cpts}
#' @param bandwidths.left left parts of the bandwidths
#' @param bandwidths.right right parts of the bandwidths
#' @param method how the asymmetric bandwidths are created;
#' possible values are
#' \itemize{
#'    \item{"cartesian" create all bandwidths in the Cartesian product of 
#'           bandwidths.left and bandwidths.right}
#'    \item{"concatenate" join bandwidths.left and bandwidths.right element-wise} 
#' }
#' @param max.unbalance a numeric value for the maximal ratio between maximal and minimal bandwidth,
#' \code{1 <= max.unbalance <= Inf}; use iff method="cartesian"
#' @return S3 \code{multiscale.grid} object to be used in the \code{multiscale.grid} function
#' @examples multiscale.grid(c(10,15,30))
#' multiscale.grid((1:5)*5, max.unbalance=2)
#' @export
multiscale.grid <- function(bandwidths.left, bandwidths.right=bandwidths.left, 
                            method="cartesian", max.unbalance=4) {
  stopifnot(bandwidths.left > 0)
  stopifnot(bandwidths.right > 0)
  stopifnot(max.unbalance >= 1.0)
  H.left <- integer(0)
  H.right <- integer(0)
  if (method == "cartesian") {
    for (G.left in bandwidths.left) {
      for (G.right in bandwidths.right) {
        ratio <- max(G.left, G.right) / min(G.left, G.right)
        if (ratio <= max.unbalance) {
          H.left <- c(H.left, G.left)
          H.right <- c(H.right, G.right) 
        }
      }
    }
  } else {
    stopifnot(method == "concatenate")
    stopifnot(length(bandwidths.left)==length(bandwidths.left))
    H.left <- bandwidths.left
    H.right <- bandwidths.right 
  }
  stopifnot(length(H.left)==length(H.right))
  structure(list(grid=cbind(H.left, H.right),
                 max.unbalance=max.unbalance),
            class = "multiscale.grid"
  )
}
