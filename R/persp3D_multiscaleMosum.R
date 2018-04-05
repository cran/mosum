#' 3D Visualization of multiscale MOSUM statistics
#' 
#' 3D Visualization of multiscale MOSUM statistics
#' @param x a numeric input data vector
#' @param G a symmetric bandwidth grid;
#' either an object of type \code{multiscale.grid} or an integer vector
#' of bandwidths
#' @param threshold string indicating which threshold should be used for normalization.
#' By default, \code{mosum.criticalValue} is used for each bandwidth (with significance
#' level \code{alpha}). Alternatively, it is possible to parse a user-defined function 
#' with \code{threshold.function} (see details)
#' @param mosum.args a named list containing further arguments
#' to be parsed to the respective \code{mosum} function calls, see \link[mosum]{mosum}
#' @param alpha alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}. Use iff \code{threshold="critical.value"}
#' @param threshold.function function object of form \code{function(G)}, to compute a
#' threshold of significance for different bandwidths; use iff \code{threshold="custom"}
#' @param expand expansion factor applied to the z coordinates
#' @param theta azimuthal angle defining the viewing direction
#' @param phi colatitude angle defining the viewing direction
#' @param xlab,ylab,zlab,ticktype graphical parameters
#' @param clim,NAcol coloring parameters
#' @param ... Further arguments to be passed to function call of \link[plot3D]{persp3D}
#' @return see \link[plot3D]{persp3D}
#' @details The visualization is based on \link[plot3D]{persp3D}.
#' To make the MOSUM statistics of different bandwidths visually comparable, 
#' they are rescaled.
#' Rescaling is done either by their respective critical value to significance level \code{alpha}
#' (iff \code{threshold=="critical.value"}) or by a custom value as given by \code{threshold.function}
#' (iff \code{threshold=="custom"})
#' @examples
#' \dontrun{
#' 
#' #' # If you run the example be aware that this may take some time
#' print("example may take some time to run")
#' 
#' x <- piecewiseStationary_timeSeries(model="mix")
#' G <- 10:40
#' persp3D.multiscaleMosum(x, G, mosum.args=list(boundary.extension=F))
#' }
#' @importFrom plot3D persp3D
#' @export
persp3D.multiscaleMosum <- function(x, G, mosum.args=list(), 
                                    threshold = c("critical.value", "custom")[1],
                                    alpha=0.05, threshold.function = NULL,
                                    expand=.2, theta=120, phi=20, xlab="G", 
                                    ylab="t", zlab="T/c", ticktype="detailed", 
                                    clim=c(0,1+2/3), NAcol="#800000FF", ...) {
  
  stopifnot(is.null(mosum.args$G))
  stopifnot(is.null(mosum.args$x))
  stopifnot(is.null(mosum.args$G.right))
  
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

  # Collect all statistics for visualization
  m <- array(list(), nrow(grid$grid))
  for (i in seq_len(nrow(grid$grid))) {
    G <- grid$grid[[i,1]]
    argList <- mosum.args
    argList$x <- x
    argList$G <- G
    m[[i]] <- do.call(mosum, argList)
  }
  
  zz <- t(sapply(m, function(z) z$stat / mosum.criticalValue(n, z$G, z$G, alpha)))
  xx <- grid$grid[,1]
  yy <- seq_len(length(x)) 
  persp3D(x = xx, y = yy, z=zz, expand=expand, 
          theta=theta, phi=phi, xlab=xlab, ylab=ylab, 
          zlab=zlab, ticktype=ticktype, 
          clim=clim, NAcol=NAcol, ...)
}