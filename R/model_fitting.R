##
## Functions that help with model fitting
##

#' STEP 1: ESTIMATING THE MEAN FUNCTION (Preliminary)
#' Implements two possible ways: (Non-robust) tensor product splines or robust
#' centering using quantile smoothing.
#' Alternatively we can directly specify `y` as the centered variables (recommended)
#'
#' Input:
#' @param modeling_df modeling dataframe, with columns `y`, `s` (spatial domain) and `timepoijt` (time domain)
#' @param method, str, either `tensor-product` or `robust`
#'
#' @returns
#' The fitted model.
#'
estimate_mean = function(modeling_df, method = "tensor-product") {

  if (method == "tensor-product") {
    mod = mgcv::gam(y ~ te(s, timepoint, k = c(10, 5), bs = "cr"),
                    data = as.data.frame(modeling_data),
                    method = "REML")
  } else if (method == "robust") {

    quat.env <- new.env(parent = environment(quantreg::rqss))
    assign(x="data", value = modeling_data, envir=quat.env)

    expr <- substitute(
      quantreg::rqss(
        y ~ qss(cbind(s, timepoint), constraint="N", lambda=20),
        data = data,
        tau = 0.5)
    )

    mod <- eval(expr, envir = quat.env)


  } else {
    rlog::log_info(paste("Method ", method, "has not been implemented yet."))
  }

  return (mod)

}


#' STEP 2: Estimate phi for direct usage with modeling_df, requires centering of `y` beforehand, `centered_y` needs
#' to be handed over explicitly.
#' Input:
#' @param modeling_df modeling dataframe with column `s`
#' @param centered_y vector of centered y observations
#' @param number_of_components algorithm to select the number of components, either
#'                            "kneedle" or "pve"
#' @param pve proportion of variance explained if `number_of_components = "pve"`, does not work if `robest = "MRCD"`
#' @param cov_est the robust estimator used to estimate the covariance (either `classic` or `MRCD`),
#' defaults to `MRCD`
#'
#' @returns List with elements:
#' * `phi_k` matrix of estimated eigenfunctions
#' * `lambda_k` eigenvalues of the covariance function
#' * `cov_center` in case of MRCD estimation, the additional centering term
#'
#' @export
estimate_phi = function(modeling_df, centered_y, number_of_components = "pve", pve=0.95, k=2, cov_est="MRCD") {

  stopifnot({
      cov_est %in% c("classic", "MRCD", "MRCT", "PASS")
      number_of_components %in% c("kneedle", "pve", "fixed")
      pve < 1
    })

  # Turn y into a matrix:
  s = unique(modeling_df$s)
  y = t(matrix(centered_y, nrow=length(s)))

  if (cov_est == "classic") {
    cov_matrix = 1 / (ncol(y) - 1) * crossprod(y)

    cov_center = NULL

  } else if (cov_est == "MRCD") {

    mrcd_obj = rrcov::CovMrcd(y, alpha=0.75, target="identity")
    cov_center = mrcd_obj@center

    # Obtain the original covariance matrix
    cov_matrix = (mrcd_obj@cov - mrcd_obj@rho * mrcd_obj@target) / (1 - mrcd_obj@rho)

  } else if (cov_est == "MRCT") {

    res = mrct(data = y)

    cov_center = colMeans(y[res$optimal.subset, ])
    cov_matrix = res$k * cov(y[res$optimal.subset, ])

  } else if (cov_est == "PASS") {
    rr = performPASSeigenAnalysis(fdata = y)
    phi_k = rr$efun
    lambda_k = rr$itr.ratio
    k = which(lambda_k > pve)[1] + 1
    phi_k = phi_k[1:k, ]

    return(list(phi_k = phi_k, lambda_k = lambda_k, cov_center = NULL))
  }

  ev = RSpectra::eigs_sym(cov_matrix, k=min(ncol(y), 100))

  ev$values = ev$values[ev$values >= 0]
  ev$vectors = ev$vectors[, ev$values >= 0]

  if (number_of_components == "pve") {
    k = min(which(cumsum(ev$values) / sum(ev$values) > pve))
    explained = round((cumsum(ev$values) / sum(ev$values))[k] * 100, 2)
  } else if (number_of_components == "kneedle") {
    k = kneedle::kneedle(1:length(ev$values), ev$values, concave=TRUE)[1] - 1
  } else if (number_of_components == "fixed") {
    k = k
  }

  rlog::log_info(
    glue::glue("Selected {k} FPCs explaining {round((cumsum(ev$values) / sum(ev$values) * 100)[k], 2)}% of the variance using `{number_of_components}`.")
  )

  phi_k = matrix(ev$vectors[, 1:k] * sqrt(ncol(y)), ncol=k)

  colnames(phi_k) = paste("phi_", 1:k, sep="")

  return (list(phi_k = phi_k, lambda_k = ev$values / sqrt(ncol(y)), cov_center = cov_center))
}


#' STEP 2: ESTIMATING THE FPCs from the centered data, requires the `y` observations in matrix form (columnwise)
#' @param y matrix of observations
#' @param number_of_components algorithm to select the number of components, either
#'                            "kneedle" or "pve"
#' @param pve proportion of variance explained if `number_of_components = "pve"`, does not work if `robest = "MRCD"`
#' @param cov_est the robust estimator used to estimate the covariance (either `classic` or `MRCD`),
#' defaults to `MRCD`
#'
#' @returns @returns List with elements:
#' * `phi_k` matrix of estimated eigenfunctions
#' * `lambda_k` eigenvalues of the covariance function
#' * `cov_center` in case of MRCD estimation, the additional centering term
#'
#' @export
estimate_phi_y = function(y, number_of_components = "pve", pve=0.95, k=2, cov_est="MRCD") {

  stopifnot({
    cov_est %in% c("classic", "MRCD", "MRCT", "PASS")
    number_of_components %in% c("kneedle", "pve", "fixed")
    pve < 1
  })

  # Turn y into a matrix

  if (cov_est == "classic") {
    cov_matrix = 1 / (ncol(y) - 1) * crossprod(y)

    cov_center = NULL

  } else if (cov_est == "MRCD") {

    mrcd_obj = rrcov::CovMrcd(y, alpha=0.75, target="identity")
    # cov_matrix = mrcd_obj@cov
    cov_center = mrcd_obj@center

    # Obtain the original covariance matrix
    cov_matrix = (mrcd_obj@cov - mrcd_obj@rho * mrcd_obj@target) / (1 - mrcd_obj@rho)

  } else if (cov_est == "MRCT") {

    res = mrct(data = y)

    cov_center = colMeans(y[res$optimal.subset, ])
    cov_matrix = res$k * cov(y[res$optimal.subset, ])

  } else if (cov_est == "PASS") {
    rr = performPASSeigenAnalysis(fdata = y)
    phi_k = rr$efun
    lambda_k = rr$itr.ratio
    k = which(lambda_k > pve)[1] + 1
    phi_k = phi_k[1:k, ]

    return(list(phi_k = phi_k, lambda_k = lambda_k, cov_center = NULL))
  }


  ev = RSpectra::eigs_sym(cov_matrix, k=min(ncol(y), 100))

  ev$values = ev$values[ev$values >= 0]
  ev$vectors = ev$vectors[, ev$values >= 0]

  if (number_of_components == "pve") {
    k = min(which(cumsum(ev$values) / sum(ev$values) > pve))
    explained = round((cumsum(ev$values) / sum(ev$values))[k] * 100, 2)
  } else if (number_of_components == "kneedle") {
    k = kneedle::kneedle(1:length(ev$values), ev$values, concave=TRUE)[1] - 1
  } else if (number_of_components == "fixed") {
    k = k
  }

  rlog::log_info(
    glue::glue("Selected {k} FPCs explaining {round((cumsum(ev$values) / sum(ev$values) * 100)[k], 2)}% of the variance using `{number_of_components}`.")
  )

  phi_k = matrix(ev$vectors[, 1:k] * sqrt(ncol(y)), ncol=k)

  colnames(phi_k) = paste("phi_", 1:k, sep="")

  return (list(phi_k = phi_k, lambda_k = ev$values / sqrt(ncol(y)), cov_center = cov_center))
}


#' STEP 3: ESTIMATE THE SCORES
#' Estimate phi for direct usage with modeling_df, requires centering of `y` beforehand, `centered_y` needs
#' to be handed over explicitly.
#' @param modeling_df modeling dataframe with column `s`
#' @param centered_y vector of centered y observations
#' @param phi the estimated eigenfunctions as returned by `estimate_phi` or `estimate_phi_y`
#' @param cov_center the additional centering term in case of MRCD covariance estimation
#'
#' @returns A matrix of scores with `ncol(phi)` columns
#' @export

estimate_scores = function(modeling_df, centered_y, phi, cov_center = NULL) {

  s = unique(modeling_df$s)
  y = t(matrix(centered_y, nrow=length(s)))

  if (!is.null(cov_center)) {
    y =  t(apply(y, 1, "-", cov_center))
  }

  k = ncol(phi)

  res = matrix(NA, nrow=nrow(y), ncol=k)

  for (i in 1:k) {
    res[, i] = apply(y, 1, \(x) mean(x * phi[, i]))
  }

  colnames(res) = paste("xi", 1:k, sep="")

  return (res)
}

#' STEP 3: ESTIMATE THE SCORES
#' Estimate the scores using y and the eigenfunctions
#' @param y the matrix of observations
#' @param phi the estimated eigenfunctions of the same length, i.e. dimensions `ncol(y) x k`
#' @param cov_center the additional centering term in case of MRCD covariance estimation
#'
#' @returns A matrix of scores with `ncol(phi)` columns
#' @export
estimate_scores_y = function(y, phi, cov_center = NULL) {

  k = ncol(phi)

  res = matrix(NA, nrow=nrow(y), ncol=k)

  for (i in 1:k) {
    res[, i] = apply(y, 1, \(x) mean(x * phi[, i]))
  }

  colnames(res) = paste("xi", 1:k, sep="")

  return (res)
}


#' STEP 3: MODELING THE SCORES

#' Transforms the data into the necessary structure
#' @param design_df design_df for the model
#' @param estimated_scores the estimated scores as returned by `estimate_scores` or `estimate_scores_y`
#'
#' @returns The object `X`, list of length 2 with observations in `X$x` and corresponding timepoints in `pp`
data_for_time_dynamic_fitting <- function(design_df, estimated_scores) {

  design_df$estimated_scores = estimated_scores

  tmp = design_df |>
    group_by(subject_id) |>
    group_split()

  X = list(x = list(), pp = list())
  X$x = lapply(tmp, dplyr::pull, "estimated_scores")
  X$pp = lapply(tmp, dplyr::pull, "t")
  names(X$x) <- names(X$pp) <- unique(design_df$subject_id)

  return (X)
}

#' STEP 3: MODELING THE SCORES
#'
#' Fits the time dynamics (either robustly or non-robustly); auxiliary function
#' Transforms the data into the necessary structure
#' @param y matrix of centered observations
#' @param phi the estimated eigenfunctions
#' @param design_df the design df for the problem
#' @param parametric logical, parametric or non-parametric smoothing
#' @param opt.h.cov bandwidth for non-parametric smoothing
#' @param cutoff_outliers cutoff c for flagging of outlying curves
#' @param mean_estimator function specifying how the "measure of outlyingness" should be averaged, defaults to median()
#'
#' @returns The object `X`, list of length 2 with observations in `X$x` and corresponding timepoints in `pp`

fit_time_dynamics_robust <- function(
    y,
    phi,
    design_df,
    parametric=TRUE,
    opt.h.cov=0.2,
    cutoff_outliers=2.5,
    mean_estimator=function(x) median(x, na.rm=TRUE)
) {

  estimated_scores = estimate_scores_y(y, phi)
  k = ncol(phi)

  if (isTRUE(parametric)) {

    score_df = dplyr::bind_cols(design_df, estimated_scores)

    models = lapply(1:k, \(i) robustlmm::rlmer(stats::formula(paste("xi", i, " ~ 1 + (1 + t | subject_id)", sep="")), score_df))

  } else if (isFALSE(parametric)) {

    models = lapply(1:k, \(i) {

      X = data_for_time_dynamic_fitting(design_df, estimated_scores[, i])

      efpca(
        X=X,
        opt.h.cov = opt.h.cov,
        rho.param=1e-3,
        ncov=floor(length(unique(design_df$t)) / 2), #length(unique(score_df$t)),
        max.kappa=1e3,
        prediction_grid = seq(0, 1, length.out = 101),
        cutoff_outliers = cutoff_outliers,
        mean_estimator = mean_estimator###
      )
    })

  }
  return (models)


}

#' Fit time dynamics non-robust (auxiliary function)
#' @param y matrix of centered observations
#' @param phi the estimated eigenfunctions
#' @param design_df the design df for the problem
#' @param parametric logical, parametric or non-parametric smoothing
fit_time_dynamics_nonrobust = function(
    y,
    phi,
    design_df, # df with subject_id and measurement timepoints
    parametric = TRUE
) {

  score_df = design_df |>
    dplyr::bind_cols(estimate_scores_y(y, phi))

  k = ncol(phi)


  if (isTRUE(parametric)) {
    models = lapply(1:k, \(x) lme4::lmer(stats::formula(paste("xi", x, " ~ 0 + (1 + t | subject_id)", sep="")), score_df))
  } else if (isFALSE(parametric)) {

    models = lapply(1:k, \(x) {
      xi_string = paste("xi", x, sep="")
      ydata = score_df |>
        rename(.id = subject_id, .index = t, .value = all_of(xi_string)) |>
        dplyr::select(.id, .index, .value)

      r = refund::fpca.sc(ydata=ydata, pve=0.9)
      return(r)
    })

  }

  return (models)

}

#' PREDICT TRAJECTORIES FOR THE WHOLE GRID

#' Predict the whole grid of information from the fitted model
#' @param ss Spatial domain indices
#' @param tt Complete time domain indices
#' @param subject_ids Which subjects to predict
#' @param models the models for the time dynamics (as returned by \code{fit_time_dynamics()})
#' @param design_df Design dataframe
#' @param y matrix of observations
#' @param phi the estimated eigenfunctions
#' @param mean_fn_model the model for the estimated mean function, optional, can also be performed by hand
#' @param cov_center the centering term returned in case of MRCD
#'
#' @returns The predicted y-values for all combinations of ss and tt and subject_id
#'          in a data frame.
#' @export

predict_grid <- function(
    ss,
    tt,
    subject_ids,
    models,
    design_df = NULL,
    y = NULL,
    phi,
    mean_fn_model,
    cov_center=NULL
) {


  result = data.frame(
    subject_id = rep(subject_ids, each = length(ss) * length(tt)),
    s = rep(ss, length(tt) * length(subject_ids)),
    t = rep(rep(tt, each = length(ss)), length(subject_ids))) |>
    mutate(ramp = sapply(t, \(x) min(x / 0.2, 1)))

  if (!is.null(cov_center)) {
    cov_center_df = data.frame(s = ss, cov_center = cov_center)
    result = result |> dplyr::left_join(cov_center_df, by = "s")
  }

  k = length(models)

  if (is(models[[1]], "rlmerMod") | is(models[[1]], "lm") | is(models[[1]], "lmrob")) {

    score_matrix = sapply(1:k, \(x) predict(models[[x]], newdata=result))
    colnames(score_matrix) = paste("xi", 1:k, sep="")

    result = result |>
      dplyr::bind_cols(score_matrix)

  } else if (is(models[[1]], "efpca")) {

    if (!is.null(cov_center)) y =  t(apply(y, 1, "-", cov_center))

    estimated_scores = estimate_scores_y(y, phi)

    pp = vector("list", length=k)

    for (i in 1:k) {
      # browser()
      X = data_for_time_dynamic_fitting(design_df, estimated_scores[, i])

      pp[[i]] = fitted.efpca(models[[i]], X=X)
    }

    score_matrix = as.data.frame(sapply(1:k, \(x) c(t(pp[[x]]$pred))))

    colnames(score_matrix) = paste("xi", 1:k, sep="")
    score_matrix$subject_id = rep(subject_ids, each = length(tt))
    score_matrix$t = rep(tt, length(subject_ids))

    result = result |>
      dplyr::left_join(score_matrix, by = c("subject_id", "t"))

  } else {
    rlog::log_error("Prediction not possible with given time dynamics model.")
  }

  phifn = apply(phi, 2, \(x) rep(x, length(tt) * length(subject_ids)))
  colnames(phifn) = paste("phi", 1:k, sep="")


  # Add the mean function if necessary
  if (!is.null(mean_fn_model)) {

    tmp = expand.grid(s=unique(result$s), timepoint=unique(result$t))
    tmp$mean_fn = predict(mean_fn_model, newdata=tmp)

    result = result |>
      dplyr::bind_cols(phifn) |>
      left_join(tmp, by = c("s" = "s", "t" = "timepoint"))

    ff = paste("mean_fn + ", paste("phi", 1:k, " * xi", 1:ncol(phi), collapse=" + ", sep=""), collapse=" + ")

  } else {

    result = result |>
      dplyr::bind_cols(phifn) |>
      dplyr::mutate(mean_fn = 0)

    ff = paste("phi", 1:k, " * xi", 1:ncol(phi), collapse=" + ", sep="")
  }

  ## TODO: Das ist noch nicht so super schön gelöst...
  if (!is.null(cov_center)) {
    ff = paste(ff, "+", "cov_center")
  }

  result = result |>
    dplyr::mutate(yhat = eval(parse(text = ff)))

  return(result)
}


#' Non-robust prediction
#' #' Predict the whole grid of information from the fitted model
#' @param ss Spatial domain indices
#' @param tt Complete time domain indices
#' @param subject_ids Which subjects to predict
#' @param models the models for the time dynamics (as returned by \code{fit_time_dynamics()})
#' @param phi the estimated eigenfunctions
#' @param mean_fn_model the model for the estimated mean function, optional, can also be performed by hand
#'
predict_grid_nonrobust = function(ss, tt, subject_ids, models, phi, mean_fn_model) {

  result = data.frame(
    subject_id = rep(subject_ids, each = length(ss) * length(tt)),
    s = rep(ss, length(tt) * length(subject_ids)),
    t = rep(rep(tt, each = length(ss)), length(subject_ids))
  ) |>
    mutate(ramp = sapply(t, \(x) min(x / 0.2, 1)))

  k = length(models)

  if (is(models[[1]], "lmerMod") | is(models[[1]], "lm")) {
    score_matrix = sapply(1:k, \(x) predict(models[[x]], newdata=result))
    colnames(score_matrix) = paste("xi", 1:k, sep="")

    result = result |>
      dplyr::bind_cols(score_matrix)
  }
  else if (is(models[[1]], "fpca")) {

    score_matrix = as.data.frame(sapply(1:k, \(x) c(t(models[[x]]$Yhat))))
    colnames(score_matrix) = paste("xi", 1:k, sep="")
    score_matrix$subject_id = rep(subject_ids, each = length(tt))
    score_matrix$t = rep(tt, length(subject_ids))

    result = result |>
      dplyr::left_join(score_matrix, by = c("subject_id", "t"))

  }

  phifn = apply(phi, 2, \(x) rep(x, length(tt) * length(subject_ids)))
  colnames(phifn) = paste("phi", 1:k, sep="")

  # Add the mean function if necessary
  if (!is.null(mean_fn_model)) {

    result = result |>
      dplyr::bind_cols(phifn) |>
      dplyr::mutate(mean_fn = predict(mean_fn_model, newdata=data.frame(s=result$s, timepoint=result$t) ))

    ff = paste("mean_fn + ", paste("phi", 1:k, " * xi", 1:k, collapse=" + ", sep=""), collapse=" + ")

  } else {

    result = result |>
      dplyr::bind_cols(phifn) |>
      dplyr::mutate(mean_fn = 0)

    ff = paste("phi", 1:k, " * xi", 1:k, collapse=" + ", sep="")
  }


  result = result |>
    dplyr::mutate(yhat = eval(parse(text = ff)))

  return(result)
}

