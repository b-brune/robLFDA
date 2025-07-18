% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fitting.R
\name{fit_time_dynamics_robust}
\alias{fit_time_dynamics_robust}
\title{Fits the time dynamics for the k components}
\usage{
fit_time_dynamics_robust(
  y,
  phi,
  design_df,
  parametric = TRUE,
  opt.h.cov = 0.2,
  cutoff_outliers = 2.5,
  mean_estimator = function(x) mean(x, na.rm = TRUE)
)
}
\description{
Fits the time dynamics for the k components
}
