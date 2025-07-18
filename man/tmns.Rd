% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efpca_utility.R
\name{tmns}
\alias{tmns}
\title{Trimmed mean}
\usage{
tmns(x, alpha)
}
\arguments{
\item{x}{numeric vector}

\item{alpha}{real number between 0 and 1}
}
\value{
the mean of the smallest (1-alpha)\% values in \code{x}
}
\description{
This function returns the mean of the smallest
(1-alpha)\% of its arguments.
}
