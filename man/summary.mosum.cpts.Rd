% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mosum.R
\name{summary.mosum.cpts}
\alias{summary.mosum.cpts}
\title{Summary of change points estimated by MOSUM procedure}
\usage{
\method{summary}{mosum.cpts}(object, ...)
}
\arguments{
\item{object}{a \code{mosum.cpts} object}

\item{...}{not in use}
}
\description{
Summary method for objects of class \code{mosum.cpts}
}
\details{
Provide information about each estimated change point, 
including the bandwidths used for its estimation, associated p-value and (scaled) jump size;
if \code{object$do.confint=TRUE}, end points of the pointwise and uniform confidence intervals
are also provided.
}
\examples{
x <- testData(lengths = rep(100, 3), means = c(0, 5, -2), sds = rep(1, 3), seed = 1234)$x
m <- mosum(x, G = 40, do.confint = TRUE)
summary(m)
}
