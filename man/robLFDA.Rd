% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robLFDA.R
\name{robLFDA}
\alias{robLFDA}
\title{Required}
\usage{
robLFDA(
  data,
  design_df,
  t = seq(0, 1, 0.01),
  covariance_estimator = "MRCD",
  centered = TRUE,
  score_method = "non-parametric",
  number_of_components = "pve",
  k_scores_pve = TRUE,
  opt.h.cov = 0.1,
  pve = 0.9,
  cutoff_score_outliers = 2.5,
  mean_estimator = function(x) median(x, na.rm = TRUE),
  fixed_number_of_components = 10
)
}
\arguments{
\item{design_df}{... data.frame which includes the design of the experiment
with columns subject_id and t (measurement time)}

\item{t}{... vector, time grid that we are estimating the time domain on}

\item{covariance_estimator}{... the covariance estimator used for estiamtion of the eigenfunctions}

\item{centered}{... is the data centered or not?}

\item{score_method}{... either "parametric" (for mixed effects model linear in time)
or "non-parametric" (for a functional data analysis approach)}

\item{modeling_data}{... data.frame which contains the data in long format and has
columns subject_id, s (indexing in frequency domain),
timepoint (measurement time),
y (observed functions)}

\item{cv_scores}{... logical, cross validation for the score covariance smoothing?}
}
\value{
An object of class \code{robLFDAmod}
}
\description{
Required
}
