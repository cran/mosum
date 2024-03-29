% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bandwidth.R
\name{bandwidths.default}
\alias{bandwidths.default}
\title{Default choice for the set of multiple bandwidths}
\usage{
bandwidths.default(n, d.min = 10, G.min = 10, G.max = min(n/2, n^(2/3)))
}
\arguments{
\item{n}{integer representing the sample size}

\item{d.min}{integer for the minimal mutual distance of change points that can be expected}

\item{G.min}{integer for the minimal allowed bandwidth}

\item{G.max}{integer for the maximal allowed bandwidth}
}
\value{
an integer vector of bandwidths
}
\description{
Create bandwidths according to a default function of the sample size
}
\details{
Returns an integer vector of bandwidths (G_1,...,G_m), 
with G_0 = G_1 = max(\code{G.min}, 2/3*\code{d.min}), G_{j+1} = G_{j-1} + G_j (for j = 1, ..., m-1)
and m satisfying G_m <= \code{G.max} while G_{m+1} > \code{G.max}.
}
\examples{
bandwidths.default(1000, 10, 10, 200)
}
\references{
A. Meier, C. Kirch and H. Cho (2021)
mosum: A Package for Moving Sums in Change-point Analysis.
\emph{Journal of Statistical Software}, Volume 97, Number 8, pp. 1-42.
<doi:10.18637/jss.v097.i08>.

H. Cho and C. Kirch (2022) Two-stage data segmentation permitting multiscale change points, heavy tails and dependence. \emph{Annals of the Institute of Statistical Mathematics}, Volume 74, Number 4, pp. 653-684.
}
