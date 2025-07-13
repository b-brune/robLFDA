#' @export
classicLFDA <- function(
    modeling_data,
    design_df,
    centered=TRUE,
    s = seq(0, 1, length.out=101),
    t = seq(0, 1, length.out=21),
    score_method="non-parametric",
    number_of_components="pve",
    k_scores_pve=FALSE,
    pve=0.9
) {



  ## ------- STEP 1 -----------------------------------------------------------

  if (isFALSE(centered)) {
    rlog::log_info("Estimating the mean function")
    mean_fn = estimate_mean(modeling_data, method = "tensor-product")

    # Add the fitted values and centered data to the modeling_data df
    modeling_data$fitted = fitted(mean_fn)
    modeling_data$y_centered = with(modeling_data, y - fitted)

  } else {
    mean_fn = NULL
    modeling_data$fitted = 0
    modeling_data$y_centered = modeling_data$y
  }

  ## ------- STEP 1 -----------------------------------------------------------


  y = t(matrix(modeling_data$y_centered, nrow=length(s)))


  ## ------- STEP 2 -----------------------------------------------------------
  rlog::log_info("Estimate the covariance function")
  tmp = estimate_phi_y(
    y,
    cov_est = "classic",
    number_of_components=number_of_components,
    pve=pve
  )
  phi_k = tmp$phi_k
  cov_center = tmp$cov_center

  ## ------ STEP 3 ------------------------------------------------------------
  rlog::log_info("Estimating the scores")
  scores = estimate_scores_y(y, phi_k)

  if (score_method == "parametric") {

    rlog::log_info("Applying parametric smoothing for score functions.")
    time_models = fit_time_dynamics_nonrobust(
      y=y,
      phi=phi_k,
      design_df=design_df,
      parametric=TRUE
    )

  } else if (score_method == "non-parametric") {
    rlog::log_info("Applying non-parametric smoothing for score functions.")

    time_models = fit_time_dynamics_nonrobust(
      y=y,
      phi=phi_k,
      design_df=design_df,
      parametric=FALSE
    )

  }


  ## --------------------------------------------------------------------------

  result = list(
    mean_model = mean_fn,
    eigenfn = phi_k,
    scores_raw = scores,
    scores_fitted_models = time_models,
    # predicted_values = pred,
    design_df = design_df,
    ss = s,
    tt = t
  )

  class(result) = "LFDAmod"

  return (
    result
  )

}

#' @export
predict.LFDAmod <- function(object, ...) {

  pred = predict_grid_nonrobust(
    ss = object$ss,
    tt = object$tt,
    subject_ids = unique(object$design_df$subject_id),
    models=object$scores_fitted_models,
    phi=object$eigenfn,
    mean_fn_model=object$mean_model
  )

  return (pred)
}
