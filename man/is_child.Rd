% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{is_child}
\alias{is_child}
\title{Is index i_child a child of index i_parent?
ASSERT: i_child is of the form (i_parent XOR i_help),
        with i_help having exactly one non-zero bit}
\usage{
is_child(i_child, i_parent)
}
\description{
Is index i_child a child of index i_parent?
ASSERT: i_child is of the form (i_parent XOR i_help),
        with i_help having exactly one non-zero bit
}
\keyword{internal}
