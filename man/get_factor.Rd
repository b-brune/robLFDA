% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efpca_utility.R
\name{get_factor}
\alias{get_factor}
\title{Calculate the outlyingness for a pair of observations given a covariance function
Get the standardized residuals}
\usage{
get_factor(x, t, sigma)
}
\arguments{
\item{x}{the two X-values X_ij and X_il}

\item{t}{the indices of the corresponding timepoints t_ij and t_il in the grid}
}
\value{
The standardized squared and squarerooted residual as described in above formulas
}
\description{
Calculate the outlyingness for a pair of observations given a covariance function
Get the standardized residuals
}
