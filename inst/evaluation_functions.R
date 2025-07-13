### Evaluation functions

postprocess_rmse_df = function(rmse_df) {
  
  if ("outlying_curve" %in% colnames(rmse_df)) {
    rmse_list = list(
      
      rmse_ocFALSE_osFALSE_obsFALSE = rmse_df |>
        filter(!outlying_curve, !outlying_subject, !observed) |>
        pull(err),
      
      rmse_ocFALSE_osFALSE_obsTRUE = rmse_df |>
        filter(!outlying_curve, !outlying_subject, observed) |>
        pull(err),
      
      rmse_ocFALSE_osTRUE_obsFALSE = rmse_df |>
        filter(!outlying_curve, outlying_subject, !observed) |>
        pull(err),
      
      rmse_ocFALSE_osTRUE_obsTRUE = rmse_df |>
        filter(!outlying_curve, outlying_subject, observed) |>
        pull(err),
      
      rmse_ocTRUE_osTRUE_obsTRUE = rmse_df |>
        filter(outlying_curve, outlying_subject, observed) |>
        pull(err)
    )
    
  } else {
    rmse_list = list(
      
      rmse_osFALSE_obsFALSE = rmse_df |>
        filter(!outlying_subject, !observed) |>
        pull(err),
      
      rmse_osFALSE_obsTRUE = rmse_df |>
        filter(!outlying_subject, observed) |>
        pull(err),
      
      rmse_osTRUE_obsFALSE = rmse_df |>
        filter(outlying_subject, !observed) |>
        pull(err),
      
      rmse_osTRUE_obsTRUE = rmse_df |>
        filter(outlying_subject, observed) |>
        pull(err)
    )
    
  }
  
  return(rmse_list)
}

get_rmses <- function(model, instance) {
  #' Calulates the rmses for the fitted model
  #' in and out of sample

  uncontaminated = create_data_matrices(
    score_matrix_list = instance$raw_data$score_matrix, 
    basis_functions = instance$raw_data$basis_functions, 
    ss = instance$raw_data$ss, 
    tt = instance$raw_data$tt, 
    mean_function= function(s, t) 0
  )
  
  
  y_raw = do.call("cbind", lapply(uncontaminated, c))
  
  # Get predicted values:
  pred = predict(model)
  outlying_subjects = instance$outlying_subjects
  
  if (instance$whole_subject) {
    if (instance$outlier_type == "amplitude") {
      y_raw[, outlying_subjects] = y_raw[, outlying_subjects] * instance$outlier_size
    }
    if (instance$outlier_type == "shift") {
      y_raw[, outlying_subjects] = y_raw[, outlying_subjects] + instance$outlier_size
    }
  }
  
  if (isFALSE(instance$whole_subject)) {
    results_table = suppressMessages(
      pred |> 
        left_join(instance$modeling_df, by=c("subject_id", "s", "t"="timepoint")) |>
        left_join(instance$design_df, by = c("subject_id", "t")) |>
        replace_na(list(outlying_curve=FALSE, observed=FALSE)) |> ### STIMMT DAS HIER???
        mutate(outlying_subject = subject_id %in% outlying_subjects) |>
        mutate(y_raw = c(y_raw)) |>
        mutate(y_corr = coalesce(y, y_raw)) |>
        group_by(subject_id, observed, outlying_curve, outlying_subject, t) |>
        summarize(rmse = mean((y_corr - yhat)^2)) |>
        ungroup() |>
        group_by(outlying_curve, outlying_subject, observed) |>
        summarize(err = sqrt(mean(rmse))) |>
        ungroup()
    )
  } else {
    results_table = suppressMessages(
      pred |> 
        left_join(instance$modeling_df, by=c("subject_id", "s", "t"="timepoint")) |>
        left_join(instance$design_df, by = c("subject_id", "t")) |>
        replace_na(list(observed=FALSE)) |> ### STIMMT DAS HIER???
        mutate(outlying_subject = subject_id %in% outlying_subjects) |>
        mutate(y_raw = c(y_raw)) |>
        mutate(y_corr = coalesce(y, y_raw)) |>
        group_by(subject_id, observed, outlying_subject, t) |>
        summarize(rmse = mean((y_corr - yhat)^2)) |>
        ungroup() |>
        group_by(outlying_subject, observed) |>
        summarize(err = sqrt(mean(rmse))) |>
        ungroup()
    )
  }
  
  return(postprocess_rmse_df(results_table))
}


rmse <- function(a, b) {
  if (length(a) != length(b)) {
    stop("`a` and `b` need to be the same length.")
  }
  
  return(
    sqrt(mean((a - b)^2))
  )
}




#' Calculate the angle between two functions
#' @param phi1 vector/function 1
#' @param phi2 vector/function 2
#' 
#' @return The angle between \code{phi1} and \code{phi2}
#' @details The angle is calculated as 
#' \eq{acos(<phi1, phi2> / (||phi1|| ||phi2||))}
#' 
#' @export

angle <- function(phi1, phi2) {
  
  if (length(phi1) != length(phi2)) {
    stop("`phi1` and `phi2` need to be the same length.")
  }
  
  phi1_norm = sqrt(sum(phi1^2))
  phi2_norm = sqrt(sum(phi2^2))
  dotprod = sum(phi1 * phi2)
  
  angle = acos(min(sum(phi1 * phi2)  / (phi1_norm * phi2_norm), 1))
  
  return (angle)
}


get_angles <- function(model, instance) {
  phi = model$eigenfn
  k = ncol(phi)
  
  phi_true = instance$raw_data$basis_functions[, 1:k]
  
  sapply(1:k, \(i) min(angle(phi[, i], phi_true[, i]), angle(-phi[, i], phi_true[, i]) ))
}

###
### Outlier detection analysis
###

#' Calculates a table of outlying curves (i.e. those recognized using
#' the spike detection algorithm)
get_outlying_curve_confusion_matrix <- function(model, instance, method = "or", cutoff = 2) {
  
  outlyingness_matrix <- matrix(unlist(lapply(model$scores_fitted_models, "[[", "outlyingness")), ncol=length(model$scores_fitted_models))

  if (method == "or") {
    flag = rowSums(outlyingness_matrix > cutoff) >= 1
  } else if (method == "and") {
    flag = rowSums(outlyingness_matrix > cutoff) == ncol(outlyingness_matrix)
  } else if (method == "combined") {
    flag = rowMeans(outlyingness_matrix) > 2
  }
    
  conf_mat = table(
    true = instance$design_df$outlying_curve, 
    flagged = flag
    )
  
  return(conf_mat)
}


#' Calculates a table of outlying curves (i.e. those recognized using
#' the spike detection algorithm)
detect_outlying_scorefunctions <- function(model, alpha=0.05) {
  
  k = length(model$scores_fitted_models)
  score_models = model$scores_fitted_models
  estimated_scores = robLFDA:::estimate_scores_y(model$y, model$eigenfn)
  lambda_k = model$lambda_k
  norms_score = matrix(nrow=length(unique(model$design_df$subject_id)), ncol=k)
  norms_orth = matrix(nrow=length(unique(model$design_df$subject_id)), ncol=k)
  
  cutoffs_score_bonferroni = vector("numeric", length = k)
  cutoffs_score_separate = vector("numeric", length = k) 
  cutoffs_orth_bonferroni = vector("numeric", length = k)
  cutoffs_orth_separate = vector("numeric", length = k) 
  
  for (i in 1:k) {
    
    X = robLFDA:::data_for_time_dynamic_fitting(model$design_df, estimated_scores[, i])  
    pred = fitted(score_models[[i]], X=X, pve=0.95)
    
    eigen_tmp = eigen(score_models[[i]]$cov.fun) 
    
    ev = eigen_tmp$values/ nrow(eigen_tmp$vectors) # normalization of eigenvalues
    
    ev_score = ev[1:ncol(pred$xis)]
    ev_orth = ev[ev > 1e-6][(ncol(pred$xis)+1):length(ev[ev > 1e-6])]
    
    norms_score[, i] = sapply(1:nrow(pred$pred), \(i) mean(pred$pred[i, ]^2)) / lambda_k[i]
    
    norms_orth[, i] = sapply(1:nrow(pred$pred), \(j) mean((X$x[[j]] - pred$pred[j, score_models[[i]]$grid %in% X$pp[[j]]])^2))
    
    # Monte Carlo simulation of the quantiles
    qs_score = quantile(replicate(10000, sum(ev_score / lambda_k[i] * rchisq(length(ev_score), df=1))), c(1 - alpha/k, 1-alpha))
    qs_orth = quantile(replicate(10000, sum(ev_orth * rchisq(length(ev_orth), df=1))), c(1 - alpha/k, 1-alpha))
    
    cutoffs_score_bonferroni[i] = qs_score[1]
    cutoffs_score_separate[i] = qs_score[2]
    
    cutoffs_orth_bonferroni[i] = qs_orth[1]
    cutoffs_orth_separate[i] = qs_orth[2]
  }
  
  
  calc_flags <- function(norms, cutoffs_separate, cutoffs_bonferroni) {
    flag_individual = sapply(1:k, \(i) norms[, i] > cutoffs_separate[i])
    
    list(
      # bonferroni method
      bonferroni = rowSums(norms) > sum(cutoffs_bonferroni),
      # flagged at least one individual
      min1 = rowSums(flag_individual) >= 1,
      # flagged both individuals separately
      both = rowSums(flag_individual) == 2 
    )
  }
  
  score_results = calc_flags(norms_score, cutoffs_score_separate, cutoffs_score_bonferroni)
  orth_results = calc_flags(norms_orth, cutoffs_orth_separate, cutoffs_orth_bonferroni)
  
  
  return(list(score = score_results, orth = orth_results))#bonferroni = flag_bonferroni, min1 = flag_min1, both=flag_both))
}


get_confusion_matrices <- function(model, instance, alpha = 0.05, which = "score") {
  flags = detect_outlying_scorefunctions(
    model = model, # instance = instance, 
    alpha=alpha)
  
  if (which == "score") {
    outl = (1:length(flags$score$bonferroni)) %in% instance$outlying_subjects
    flags$bonferroni = flags$score$bonferroni
    flags$min1 = flags$score$min1
    flags$both = flags$score$both
  } else if (which == "orth") {
    outl = (1:length(flags$orth$bonferroni)) %in% instance$outlying_subjects
    flags$bonferroni = flags$orth$bonferroni
    flags$min1 = flags$orth$min1
    flags$both = flags$orth$both
  }
  
  return(
    list(
      bonferroni = table(outl, flags$bonferroni),
      min1 = table(outl, flags$min1),
      both = table(outl, flags$both)
    ))
}

## Calculate different diagnostic measures from confusion matrices:
TPR = function(confmat) {
  if (is.null(confmat)) return(NA)
  if (length(confmat) != 4) return (0)
  
  return(confmat["TRUE", "TRUE"] / (confmat["TRUE", "FALSE"] + confmat["TRUE", "TRUE"]))
}

FPR = function(confmat) {
  
  if (is.null(confmat)) return(NA)
  if (length(confmat) != 4) return (0)
  
  return(confmat["FALSE", "TRUE"] / (confmat["FALSE", "TRUE"] + confmat["FALSE", "FALSE"]))
}

F1_score = function(confmat) {
  if (is.null(confmat)) return(NA)
  if (length(confmat) != 4) return (0)
  
  precision <- diag(confmat) / colSums(confmat)
  recall <- diag(confmat) / rowSums(confmat)
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  #Assuming that F1 is zero when it's not possible compute it
  f1[is.na(f1)] <- 0
  
  #Binary F1 or Multi-class macro-averaged F1
  return(f1["TRUE"])
}

precision = function(confmat) {
  if (is.null(confmat)) return(NA)
  if (length(confmat) != 4) return (0)
  
  precision <- confmat["TRUE", "TRUE"] / (confmat["TRUE", "TRUE"] + confmat["FALSE", "TRUE"])
  
  return(precision)
}

recall = function(confmat) {
  if (is.null(confmat)) return(NA)
  if (length(confmat) != 4) return (0)
  
  recall <- confmat["TRUE", "TRUE"] / (confmat["TRUE", "TRUE"] + confmat["TRUE", "FALSE"])
  
  return(recall)
}


####
#### Accuracy of the reconstructed score functions
####

score_function_rmse <- function(model, instance) {
  
  k = length(model$scores_fitted_models)
  estimated_scores = robLFDA:::estimate_scores_y(model$y, model$eigenfn)
  
  
  mses = vector("list", k)
  rel_mses = vector("list", k)
  
  for (mod in 1:k) {
    X = robLFDA:::data_for_time_dynamic_fitting(model$design_df, estimated_scores[, mod])  
    
    pp = fitted(model$scores_fitted_models[[mod]], X)
    
    pred = t(pp$pred)
    true = sapply(instance$raw_data$score_matrix, function(x) x[, mod])
    
    mses[[mod]] = vector("numeric", length(X$x))
    rel_mses[[mod]] = vector("numeric", length(X$x))
    
    for (i in 1:length(X$x)) {
      mses[[mod]][i] = min(mean((pred[, i] - true[, i])^2), mean((pred[, i] + true[, i])^2))
      rel_mses[[mod]][i] = mses[[mod]][i] / mean(true[, i]^2)
    } 
  }  
  
  mean_reconstruction_errors_outlying_subjects = unlist(lapply(mses, \(x) sqrt(mean(x[instance$outlying_subjects]))))
  mean_reconstruction_errors_non_outlying_subjects = unlist(lapply(mses, \(x) sqrt(mean(x[-instance$outlying_subjects]))))
  
  relative_mean_reconstruction_errors_outlying_subjects = unlist(lapply(rel_mses, \(x) sqrt(mean(x[instance$outlying_subjects]))))
  relative_mean_reconstruction_errors_non_outlying_subjects = unlist(lapply(rel_mses, \(x) sqrt(mean(x[-instance$outlying_subjects]))))
  
  return(
    list(
      outlying = mean_reconstruction_errors_outlying_subjects, 
      inlying = mean_reconstruction_errors_non_outlying_subjects,
      rel_outlying = relative_mean_reconstruction_errors_outlying_subjects,
      rel_inlying = relative_mean_reconstruction_errors_non_outlying_subjects
    )
  )
}




score_function_rmse_classic <- function(model, instance) {
  k = length(model$scores_fitted_models)
  
  mses = vector("list", k)
  for (mod in 1:k) {

    pred = t(model$scores_fitted_models[[mod]]$Yhat)
    true = sapply(instance$raw_data$score_matrix, function(x) x[, mod])
    
    mses[[mod]] = vector("numeric", ncol(pred))
    
    for (i in 1:ncol(pred)) {
      mses[[mod]][i] = min(mean((pred[, i] - true[, i])^2), mean((pred[, i] + true[, i])^2))
    }
  }  
  
  mean_reconstruction_errors_outlying_subjects = unlist(lapply(mses, \(x) sqrt(mean(x[instance$outlying_subjects]))))
  mean_reconstruction_errors_non_outlying_subjects = unlist(lapply(mses, \(x) sqrt(mean(x[-instance$outlying_subjects]))))
  
  return(
    list(
      outlying = mean_reconstruction_errors_outlying_subjects, 
      inlying = mean_reconstruction_errors_non_outlying_subjects
    )
  )
}


