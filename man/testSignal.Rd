% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_data.R
\name{testSignal}
\alias{testSignal}
\title{Piecewise constant test signal}
\usage{
testSignal(
  model = c("custom", "blocks", "fms", "mix", "stairs10", "teeth10")[1],
  lengths = NULL,
  means = NULL,
  sds = NULL
)
}
\arguments{
\item{model}{a string indicating from which model a realisation is to be generated;
possible values are "custom" (for user-specified model
using \code{lengths}, \code{means} and \code{sds}), and
"blocks", "fms", "mix", "stairs10", "teeth10" (for the referenced test signals)}

\item{lengths}{use iff \code{model = "custom"}; an integer vector for the lengths of the piecewise stationary segments}

\item{means}{use iff \code{model = "custom"}; a numeric vector for the means of the piecewise stationary segments}

\item{sds}{use iff \code{model = "custom"}; a numeric vector for the deviation scaling of the piecewise stationary segments.}
}
\value{
a list containing the following entries:
\itemize{
  \item{mu_t}{ mean vector of piecewise stationary model time series}
  \item{sigma_t}{ deviation scaling vector of piecewise stationary model time series}
}
}
\description{
Produce vectors of mean and dispersion values for generating piecewise stationary time series.
}
\details{
See Appendix B in the reference for details about the test signals.
}
\references{
P. Fryzlewicz (2014)
Wild Binary Segmentation for Multiple Change-Point Detection.
\emph{The Annals of Statistics}, Volume 42, Number 6, pp. 2243-2281.
}
\keyword{internal}
