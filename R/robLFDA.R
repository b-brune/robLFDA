#' Required
#' @param modeling_data ... data.frame which contains the data in long format and has
#'                      columns subject_id, s (indexing in frequency domain),
#'                      timepoint (measurement time),
#'                      y (observed functions)
#' @param design_df ... data.frame which includes the design of the experiment
#'                      with columns subject_id and t (measurement time)
#' @param t ... vector, time grid that we are estimating the time domain on
#' @param covariance_estimator ... the covariance estimator used for estiamtion of the eigenfunctions
#' @param centered ... is the data centered or not?
#' @param score_method ... either "parametric" (for mixed effects model linear in time)
#'                     or "non-parametric" (for a functional data analysis approach)
#' @param cv_scores ... logical, cross validation for the score covariance smoothing?
#'
#' @returns An object of class `robLFDAmod`
#'
#' @export

robLFDA <- function(
    data,
    design_df,
    t = seq(0, 1, 0.01),
    covariance_estimator="MRCD",
    centered=TRUE,
    score_method="non-parametric",
    number_of_components="pve",
    k_scores_pve = TRUE, # should the number of components be determined dynamically based on proportion of variance explained?
    opt.h.cov=0.1,
    pve=0.9,
    cutoff_score_outliers=2.5,
    mean_estimator = function(x) median(x, na.rm=TRUE),
    fixed_number_of_components=10
) {

  ## ------- STEP 1 -----------------------------------------------------------

  s = data$s |> unique() |> sort()

  if (isFALSE(centered)) {
    rlog::log_info("Estimating the mean function")

    mean_fn_rob = estimate_mean(data, method="robust") # Hier wird data frame format benötigt um das modellfitting sauber hinzubekommen

    # Add the fitted values and centered data to the modeling_data df # Hier auch
    data$fitted_rob = fitted(mean_fn_rob)
    data$y_centered_rob = with(data, y - fitted_rob)

  } else {
    mean_fn_rob = NULL
    data$fitted_rob = 0
    data$y_centered_rob = data$y
  }

  y_rob = t(matrix(data$y_centered_rob, nrow=length(s))) # Hier wird die Matrix erstellt mit der dann weitergerechnet wird


  ## ------- STEP 2 -----------------------------------------------------------
  rlog::log_info("Estimate the covariance matrix robustly")
  tmp = estimate_phi_y(
    y=y_rob,
    cov_est=covariance_estimator,
    number_of_components=number_of_components,
    pve=pve,
    k=fixed_number_of_components
  )
  phi_k_rob = tmp$phi_k
  cov_center = tmp$cov_center

  y_rob = t(apply(y_rob, 1, "-", cov_center))


  ## ------ STEP 3 ------------------------------------------------------------
  rlog::log_info("Model fitting for the score functions")
  scores = estimate_scores_y(y_rob, phi_k_rob)

  if (score_method == "parametric") {

    ######## DAS IST NOCH DIE ALTE VERSION
    rlog::log_info("Applying parametric smoothing for score functions.")
    time_models = fit_time_dynamics_robust(y_rob, phi_k_rob, design_df, parametric=TRUE)

  } else if (score_method == "non-parametric") {
    rlog::log_info("Applying non-parametric smoothing for score functions.")

    time_models = fit_time_dynamics_robust(
      y=y_rob,
      phi=phi_k_rob,
      design_df=design_df,
      parametric=FALSE,
      opt.h.cov=opt.h.cov,
      cutoff_outliers=cutoff_score_outliers,
      mean_estimator=mean_estimator
    )

  }

  ## --------------------------------------------------------------------------

  result = list(
    mean_model = mean_fn_rob,
    eigenfn = phi_k_rob,
    lambda_k = tmp$lambda_k,
    cov_center = cov_center,
    scores_raw = scores,
    scores_fitted_models = time_models,
    y = y_rob,
    ss = s,
    tt = t,
    design_df = design_df
  )

  class(result) = "robLFDAmod"

  return (
    result
  )

}

#' @export
predict.robLFDAmod <- function(
    object, ...
) {

  predict_grid(
    ss = object$ss,
    tt = object$tt,
    subject_ids = unique(object$design_df$subject_id),
    models=object$scores_fitted_models,
    design_df = object$design_df,
    y = object$y,
    phi=object$eigenfn,
    mean_fn_model=object$mean_model,
    cov_center=object$cov_center,
    ...
  )

}

