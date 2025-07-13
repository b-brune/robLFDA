library(data.table)
library(tidyverse)

# After connecting to the registry, do

jobs = getJobTable()
res = reduceResultsDataTable()
res = ijoin(jobs, res)


# Then feed res to the following function and save the results.
preprocess_raw_batchtools_results <- function(res) {

  for (colname in names(res$prob.pars[[1]])) {
    print(glue::glue("Unwrapping {colname}"))
    if (colname == "sample_sizes") {
      res <- res |>
        mutate(!!colname := lapply(res$prob.pars, "[[", colname))
    } else if (colname == "grid_s") {
      res <- res |>
        mutate(!!colname := unlist(lapply(lapply(res$prob.pars, "[[", colname), function(x) if (length(x) >= 150) "fine" else "rough")))
    }
    else {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$prob.pars, "[[", colname)))
    }
  }

  res <- res |> select(-prob.pars)

  for (colname in names(res$result[[1]])) {
    res <- res |>
      mutate(!!colname := lapply(res$result, "[[", colname))
  }

  res <- res |> select(-result)

  for (colname in unique(unlist(lapply(res$algo.pars, names)))) {
    print(glue::glue("Unwrapping {colname}"))

    if (colname == "opt.h.cov") {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$algo.pars, function(x) if(!is.null(x$opt.h.cov)) x$opt.h.cov else NA)))
    } else if (colname == "cutoff_score_outliers") {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$algo.pars, function(x) if(!is.null(x$cutoff_score_outliers)) x$cutoff_score_outliers else NA)))
    } else if (colname == "covariance_estimator") {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$algo.pars, function(x) if(!is.null(x$covariance_estimator)) x$covariance_estimator else NA)))
    } else if (colname == "mean_estimator") {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$algo.pars, function(x) if(!is.null(x$mean_estimator)) x$mean_estimator else NA)))
    } else if (colname == "dummy") {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$algo.pars, function(x) if(!is.null(x$dummy)) x$dummy else NA)))
    }
    else {
      res <- res |>
        mutate(!!colname := unlist(lapply(res$algo.pars, "[[", colname)))
    }
  }


  res$opt.h.cov = unlist(lapply(res$opt.h.cov, \(x) ifelse(is.null(x), NA, x)))
  res$cutoff_score_outliers = unlist(lapply(res$cutoff_score_outliers, \(x) ifelse(is.null(x), NA, x)))



  res <- res |> select(-algo.pars)

  res$number_of_curves <- unlist(lapply(res$sample_sizes, \(x) ifelse(x$n_min == 5, "5-15", "15-30")))

  check_if_fn_is_median <- function(fn) {
    test_vec = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 100)
    if (!is.function(fn)) {
      return("mean / none")
    } else if (fn(test_vec) == 1) {
      return("median")
    } else {
      return("trimmed_mean")
    }
  }


  res$mean_estimator <- unlist(lapply(res$mean_estimator, check_if_fn_is_median))

  res = bind_cols(res, batchtools::unwrap(res[, "rmses"]))

  # tmp1 <- bind_cols(res[1:60000, ], unwrap(res[1:60000, "rmses"]))
  # tmp2 <- bind_cols(res[60001:nrow(res), ], unwrap(res[60001:nrow(res), "rmses"]))
  # res <- bind_rows(tmp1, tmp2)
  #
  res$rmse_ocTRUE_osTRUE_obsTRUE = unlist(lapply(res$rmse_ocTRUE_osTRUE_obsTRUE, \(x) ifelse(is.null(x), NA, x)))

  res <- res |> select(-rmses)
  res = bind_cols(res, batchtools::unwrap(res |> select(angles)))

  res <- res |> select(-angles)

  res <- bind_cols(res, batchtools::unwrap(batchtools::unwrap(res[, "score_fn_rmses"])))

  if (!any(res$whole_subject)) {
    res <- bind_cols(res, batchtools::unwrap(res[, "confmat_outlying_subjects"]))
  }

  res <- res |> select(-c(submitted, started, done, error, mem.used, batch.id, log.file, job.hash,
                          job.name, repl, time.queued, time.running, resources, tags))

  return (res)
}
