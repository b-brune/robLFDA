library(batchtools)
library(data.table)
library(tidyverse)
library(patchwork)

# Create registry:
# reg <- makeExperimentRegistry(
#   file.dir = "robLFDA_experiments",
#   seed = 1,
#   source=c("DATA_GENERATION.R", "evaluation_functions.R"),
#   packages=c("tidyverse", "robLFDA")
# )

# Load registry (if already created)
reg <- loadRegistry(file.dir="robLFDA_experiments", writeable=TRUE)

# Specify number of nodes for paralellization:
# reg$cluster.functions <- makeClusterFunctionsSocket(ncpus = 80, fs.latency = 65)


###
### Dataset generation
###
create_dataset <- function(
    data,
    job,
    n,
    whole_subject = TRUE,
    outlier_type = "none",
    outlier_size = 3,
    outlying_subject_proportion=0.2,
    outliers_per_individual = 0.05,
    score_fn_type="fourier",
    sample_sizes = list(n_min=5, n_max=15),
    ...) {
  
  raw_data = create_n_samples(
    n=n, 
    tt=seq(0, 1, 0.01), 
    ss=seq(0, 1, 0.01), 
    sd_nonmeaningful=0.05,
    n_min=sample_sizes$n_min, n_max=sample_sizes$n_max,
    mean_function=function(s, t) 0,
    score_fn_type=score_fn_type
  )
  
  contaminated_data = contaminate_n_samples(
    raw_data, 
    whole_subject=whole_subject, 
    outlier_type=outlier_type, 
    outlier_size=outlier_size, 
    outlying_subject_proportion = outlying_subject_proportion,
    outliers_per_individual = outliers_per_individual
  )
  
  modeling_dataframe = reformat_for_modeling(contaminated_data) 

  outlying_subjects = (1:n)[contaminated_data$outlying_subjects]
  
  design_df = modeling_dataframe |> 
    dplyr::select(subject_id, timepoint) |>
    distinct() |>
    rename(t = timepoint) |>
    mutate(outlying_subject = subject_id %in% outlying_subjects)
  
  if (is.null(contaminated_data$outlying_samples)) {
    contaminated_data$outlying_samples = FALSE
  }
  
  design_df = design_df |> mutate(outlying_curve = contaminated_data$outlying_samples, observed = TRUE)
  
  return(
    list(
      design_df = design_df,
      modeling_df = modeling_dataframe,
      outlying_subjects = outlying_subjects,
      outlying_samples = contaminated_data$outlying_samples,
      raw_data = raw_data,
      whole_subject = whole_subject,
      outlier_size = outlier_size,
      outlier_type = outlier_type
    )
  )
}

addProblem(name="centered_dataset", fun=create_dataset)


###
### Algorithm: Model (classic and robust)
###

wrapper.robLFDA <- function(
    data,
    job,
    instance, 
    opt.h.cov = 0.1,
    cutoff_score_outliers = 2.5,
    mean_estimator = function(x) mean(x, na.rm = TRUE),
    ...
) {
  
  # Fit the model from the instance:
  res = robLFDA:::robLFDA(
    instance$modeling_df,
    instance$design_df,
    t=instance$raw_data$tt,
    centered=TRUE,
    covariance_estimator="MRCD",
    score_method="non-parametric",
    pve=0.95,
    number_of_components="pve",
    opt.h.cov=opt.h.cov,
    k_scores_pve=TRUE,
    cutoff_score_outliers=cutoff_score_outliers,
    mean_estimator = mean_estimator
  )
  
  # Calculate the reconstruction errors:
  rmses = get_rmses(res, instance)
  
  # Calculate angles of phi vs. \hat{phi}
  angles = get_angles(res, instance)
  
  # Calculate confusion matrices for cellwise outliers:
  confmat_outlying_curves = get_outlying_curve_confusion_matrix(
    res, instance, method = "or",
    cutoff = cutoff_score_outliers)
  confmat_outlying_curves_and = get_outlying_curve_confusion_matrix(
    res, instance, method = "and",
    cutoff = cutoff_score_outliers)
  confmat_outlying_curves_comb = get_outlying_curve_confusion_matrix(
    res, instance, method = "combined",
    cutoff = cutoff_score_outliers)
  
  # Calculate confusion matrices for score- and orthogonal distance:
  confmat_outlying_subjects = get_confusion_matrices(res, instance = instance, which="score")
  confmat_outlying_subjects_orth = get_confusion_matrices(res, instance = instance, which="orth")
  
  # Calculate reconstruction errors for score functions
  score_fn_rmses = score_function_rmse(res, instance)
  
  return (
    list(
      rmses = rmses, 
      angles= angles, 
      score_fn_rmses = score_fn_rmses, 
      confmat_outlying_curves = confmat_outlying_curves, 
      confmat_outlying_curves_and = confmat_outlying_curves_and,
      confmat_outlying_curves_comb = confmat_outlying_curves_comb,
      confmat_outlying_subjects = confmat_outlying_subjects)
  )
}

addAlgorithm(name="robLFDA", fun=wrapper.robLFDA)


wrapper.classicLFDA <- function(
    data,
    job,
    instance, 
    ...
) {
  
  result = robLFDA:::classicLFDA(
    instance$modeling_df,
    instance$design_df,
    t=instance$raw_data$tt,
    centered=TRUE,
    score_method="non-parametric",
    pve=0.95,
    number_of_components="pve"
  )
  
  ### Calculate the metrics we are interested in:
  # Reconstruction errors
  rmses = get_rmses(result, instance)
  # Angles
  angles = get_angles(result, instance)
  # Reconstruction errors for the score functions
  score_fn_rmses = score_function_rmse_classic(result, instance)
  
  return (list(rmses = rmses, angles = angles, score_fn_rmses = score_fn_rmses))
}

addAlgorithm(name="classicLFDA", fun=wrapper.classicLFDA)

#########
######### 

## 
##  Settings reported in the paper (and a few more):
##

dataset_settings <- list(
  centered_dataset = data.table::CJ(
    n=c(50, 100),
    sample_sizes=list(list(n_min=5, n_max=15), list(n_min=15, n_max=30)),
    outlier_size=c(1, 3, 5, 10, 20),
    outlying_subject_proportion = c(0.05, 0.1, 0.2),
    outliers_per_individual = c(0.05, 0.1, 0.2),
    score_fn_type=c("fourier"),
    outlier_type=c("amplitude", "shape_score"),
    whole_subject=c(TRUE, FALSE),
    sorted=FALSE)[!(outlier_type == "shape_score" & outlier_size > 1)][!(whole_subject & outliers_per_individual > 0.05)]
)


algo_settings <- list(
  robLFDA = CJ(
    opt.h.cov=c(0.05, 0.1, 0.2), 
    cutoff_score_outliers=2, 
    mean_estimator = list(function(x) median(x, na.rm = TRUE)), 
    sorted = FALSE), 
  classicLFDA = data.table(dummy = 0)
)

addExperiments(dataset_settings, algo_settings, repls=100)

# Test if the job works:
testJob(1)

ids = getJobTable()[, .(job.id)]
ids[, chunk := chunk(job.id, n.chunks=80)]

# Submit jobs:
# submitJobs(ids)


