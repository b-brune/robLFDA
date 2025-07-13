# Data generation

# Use simulation setting from Park and Staicu?

# Y_{ij}(s) = \mu(s, T_ij) + \xi_1(T_ij)\phi_1(s) + \xi_2(s)(T_ij)\phi_2(s)

mu <- function(s, t) {
  # 0
  1 + 2 * s + 3 * t + 4 * s * t
}


xi_np = function(t, mean1=0, mean2=0, sd1=sqrt(3), sd2=sqrt(1.5)) {
  
  xi = rnorm(n = 1, mean = mean1, sd = sd1) * sqrt(2) * cos(2 * pi * t) + 
    rnorm(n = 1, mean = mean2, sd = sd2) * sqrt(2) * sin(2 * pi * t) 
  
  return (xi)
}

xi_linear = function(t, mean1=0, mean2=0, sd1=sqrt(3), sd2=sqrt(1.5)) {
  
  xi = rnorm(n = 1, mean = mean1, sd = sd1) * sqrt(3) * (2*t - 1) + rnorm(n = 1, mean = mean2, sd = sd2) * 1 
  
  return (xi)
}



xi_rem = function(t, mean=0, sd=1) {
  
  bik = rnorm(2, mean = mean, sd = sd)
  xi = bik[1] + bik[2] * t
  
  return (xi)
}

xi_exp = function(t, rho=0.9, lam=4.5) {
  
  corr.mat = rho^(abs(outer(t, t, "-")))
  xi = MASS::mvrnorm(1, rep(0, length(t)), Sigma = lam * corr.mat)
  
  return (xi)
}



create_data_matrices <- function (
    score_matrix_list, 
    basis_functions,
    mean_function,
    tt, 
    ss
    ) {
  #' Create the data matrices
  #' @param score_matrix_list List of matrices of scores geneated by xi-functions or
  #'                          white noise
  #' @param basis_functions basis functions returned by fda::create_basis_functions
  #' @param mean_function mean function evaluated in tt and ss
  #' @param tt time resolution
  #' @param ss frequency resolution
  #'
  #' @returns list of data matrices
  
  mu_eval = sapply(tt, \(x) mean_function(ss, x)) 
  
  result <- lapply(1:length(score_matrix_list), 
                   \(i) mu_eval + Reduce("+", 
                               lapply(1:ncol(score_matrix_list[[i]]), \(k) t(kronecker(t(basis_functions[, k]), score_matrix_list[[i]][, k])))
                   )
  )
  
  # Each basis function product can be written as (t(phi_k(.)) (x) xi_{ik}(.)) ((x) = kronecker product)
  
  return (result)
}


sample_subject <- function(
    data_matrix,
    tt,
    n_min=20, 
    n_max=40) {
  #' Samples random measurements from the whole data
  #' @param data_matrix as returned by create_data_matrices
  #' @param n_min minimum number of samples
  #' @param n_max maximum number of samples
  #'
  #' @returns list of length two with elements data (matrix) and timepoints (sampling 
  #'          timepoints corresponding to the columns)

  n_samples = runif(1, n_min, n_max)
  timepoints = sort(sample(tt, n_samples))
  flag = tt %in% timepoints
  
  return ( 
    list(data = data_matrix[, flag],
         timepoints = tt[flag]
    )
  )
}


sample_data <- function(
    tt, 
    data_matrix_list,
    ...
) {
  #' Actually applies the sample_subject function to the data_matrix_list
  #' @param tt numeric vector of timepoints
  #' @param data_matrix_list as returned by create_data_matrices
  #' @param ... additional parameters handed over to sample_subject function
  #' @returns list of sampled data and corresponding timepoints for each subject
  
  return(
    lapply(data_matrix_list, \(x) sample_subject(x, tt=tt, ...))  
  )
} 


create_n_samples <- function(n, tt, ss, mean_function, reduce=TRUE, n_basis = 10, n_meaningful_basis=2, sd_nonmeaningful=0.1, score_fn_type="fourier", ...) {
  #' Creates n samples
  #' @param n number of samples
  #' @param tt resolution in time domain
  #' @param ss resolution in frequency domain
  #' @param mean_function mean function evaluable in ss and tt
  #' @param reduce boolean, should the data be sampled down? defaults to `TRUE`
  #' @param n_basis number of basis functions
  #' @param n_meaningful_basis number of meaningful basis functions
  #' @param sd_nonmeaningful sd for the non-meaningful score functions (unstructured, N(0, sd_nonmeaningful^2) rv's)
  #' @param ... additional parameters handed over to sample_subject function
  #'
  #' @returns list of length three with reduced data, basis_functions and score_matrix 
  #'          (basis functions and score matrix can be used to recreate all the raw data)
  #' @example
  #' ss = seq(0, 1, 0.01)
  #' tt = seq(0, 1, 0.01)
  #' data_object = create_n_samples(10, tt, ss, mu)
  
  
  # Create basis functions (fourier or B-splines or polynomial)
  basis = fda::create.fourier.basis(rangeval=c(0, 1), nbasis=n_basis + 1)
  basis_vectors = fda::eval.basis(ss, basis)[, -1]
  
  create_score_functions = function(number, tt, score_fn_type) {
    sds1 = sqrt(3) / 1:number
    sds2 = sqrt(1.5) / 1:number
    
    n = matrix(NA, ncol = number, nrow=length(tt))
    for (i in 1:number) {
      if (score_fn_type == "fourier") {
        n[, i] = xi_np(tt, sd1=sds1[i], sd2 = sds2[i])
      } else if (score_fn_type == "linear") {
        n[, i] = xi_linear(tt, sd1=sds1[i], sd2 = sds2[i])
      }
    }
    
    return(n)
  }
  
  # Create score matrices:
  score_matrices = replicate(n, 
                             cbind(create_score_functions(n_meaningful_basis, tt=tt, score_fn_type=score_fn_type), 
                                   replicate(n_basis - n_meaningful_basis, rnorm(length(tt), mean=0, sd=sd_nonmeaningful))), 
                             simplify=FALSE)
  
  # Combine them
  data_matrices = create_data_matrices(score_matrices, basis_vectors, ss=ss, tt=tt, mean_function=mean_function)
  
  # Sample random trajectories for model fitting
  if (reduce) {
    reduced_data = sample_data(tt, data_matrices, ...)
  }
  # Tracking of basis vectors and score matrices enables reconstruction of all trajectories without having to store a 
  # huge amount of data!
  
  return (list(data = reduced_data, basis_functions = basis_vectors, score_matrix = score_matrices, ss=ss, tt=tt))
}




contaminate_n_samples <- function (
    data_object, 
    whole_subject = FALSE, 
    outlier_type,
    outlier_size,
    outlying_subject_proportion = 0.1,
    outliers_per_individual = 0.3,
    ...) {
  #' Contaminates the samples with different types of outliers
  #' @param data_object as returned from create_n_samples
  #' @param whole_subject should the whole subject be contaminated or random subset?
  #' @param outlier_type which outlier type? currently: shift, amplitude, shape_eigen
  #' (different eigenfunction), shape_score (different score function)
  #' @param outlier_size size of the outliers (for amplitude and shift)
  #' @param outlying_subject_proportion number of outlying subjects
  #' @param outliers_per_individual outlying samples per individual (ignored if
  #' whole_subject = TRUE)
  #' 
  #' @returns list of data, outlying_subjects and outlying samples if necessary
  
  
  # Now the contamination should happen
  # Grundlegend Unterscheidung zwischen whole subject und single outlying curves
  
  if (outlying_subject_proportion == 0) {
    return (c(data_object, list(outlying_subjects = rep(FALSE, length(data_object$data)), outlying_samples = NULL)))
    
  }
  
  # Sample the outlying subjects:
  outlying_subjects = sample(c(TRUE, FALSE), length(data_object$data), replace=TRUE, prob = c(outlying_subject_proportion, 1 - outlying_subject_proportion))
  
  if (!any(outlying_subjects)) {
    outlying_subjects[1] = TRUE
  }
  
  outlying_samples = c()

    
  # Single observations
  if (isTRUE(whole_subject)) {
    
    if (outlier_type == "shift") {
      
      for (o in seq_along(outlying_subjects)) {
        if (outlying_subjects[o]) {
          data_object$data[[o]]$data = data_object$data[[o]]$data + outlier_size 
          outlying_samples = c(outlying_samples, rep(TRUE, ncol(data_object$data[[o]]$data)))
        } else {
          outlying_samples = c(outlying_samples, rep(FALSE, ncol(data_object$data[[o]]$data)))
        }
      }
      
    } else if (outlier_type == "amplitude") {
      
      for (o in seq_along(outlying_subjects)) {
        if (outlying_subjects[o]) {
          data_object$data[[o]]$data = data_object$data[[o]]$data * outlier_size 
          outlying_samples = c(outlying_samples, rep(TRUE, ncol(data_object$data[[o]]$data)))
        } else {
          outlying_samples = c(outlying_samples, rep(FALSE, ncol(data_object$data[[o]]$data)))
        }
      }
      
    } else if (outlier_type == "shape_score") {
      
      # Create a second set of data to use for replacement?
      outlying_score_matrices = replicate(length(data_object$data), 
                                          cbind(
                                            replicate(2, rnorm(length(data_object$tt), mean=0, sd = 3)),
                                            replicate(8, rnorm(length(data_object$tt), mean=0, sd=0.1))), 
                                          simplify=FALSE)
      
      contaminated_data_matrices = create_data_matrices(outlying_score_matrices, 
                                                        basis_functions=data_object$basis_functions, 
                                                        mean_function = function(s, t) { return(0) },
                                                        ss=data_object$ss,
                                                        tt = data_object$tt)
      
      for (o in seq_along(outlying_subjects)) {
        
        if (outlying_subjects[o]) {
          data_object$data[[o]]$data = contaminated_data_matrices[[o]][, data_object$tt %in% data_object$data[[o]]$timepoints]  
          outlying_samples = c(outlying_samples, rep(TRUE, ncol(data_object$data[[o]]$data)))
        } else {
          outlying_samples = c(outlying_samples, rep(FALSE, ncol(data_object$data[[o]]$data)))
        }
        
      }
      
      
    } else if (outlier_type == "shape_eigen") {
      # Not yet implemented!
      
      
    }
  } 
  
  # Whole subjects
  if (isFALSE(whole_subject)) {
    outlying_samples = c()
    
    if (outlier_type == "shift") {
      for (i in seq_along(outlying_subjects)) {
        if (outlying_subjects[i]) {
          nobs = ncol(data_object$data[[i]]$data)
          outl = sample(c(TRUE, FALSE), size=nobs, prob=c(outliers_per_individual, 1 - outliers_per_individual), replace=TRUE)
          if (!any(outl)) {
            outl[sample(nobs, 1)] = TRUE
          }
          data_object$data[[i]]$data[, outl] =data_object$data[[i]]$data[, outl] + outlier_size
          outlying_samples = c(outlying_samples, outl)
        }else {
          outlying_samples = c(outlying_samples, rep(FALSE, ncol(data_object$data[[i]]$data)))
        }
      }
      
    } else if (outlier_type == "amplitude") {
      
      for (i in seq_along(outlying_subjects)) {
        
        if (outlying_subjects[i]) {
          nobs = ncol(data_object$data[[i]]$data)
          outl = sample(c(TRUE, FALSE), size=nobs, prob=c(outliers_per_individual, 1 - outliers_per_individual), replace=TRUE)
          data_object$data[[i]]$data[, outl] = data_object$data[[i]]$data[, outl] * outlier_size
          outlying_samples = c(outlying_samples, outl)
        } else {
          outlying_samples = c(outlying_samples, rep(FALSE, ncol(data_object$data[[i]]$data)))
        }
      }
      
    } else if (outlier_type == "shape_score") {
      
           # Create a second set of data to use for replacement?
      outlying_score_matrices = replicate(length(data_object$data), 
                                          cbind(
                                            replicate(2, rnorm(length(data_object$tt), mean=0, sd = 3)),
                                            replicate(8, rnorm(length(data_object$tt), mean=0, sd=0.1))), 
                                          simplify=FALSE)
      
      contaminated_data_matrices = create_data_matrices(
        outlying_score_matrices, 
        basis_functions=data_object$basis_functions, 
        mean_function = function(s, t) { return(0) },
        ss=data_object$ss,
        tt = data_object$tt
        )
      
      for (o in seq_along(outlying_subjects)) {
        if (outlying_subjects[o]) {
          nobs = ncol(data_object$data[[o]]$data)
          outl = sample(c(TRUE, FALSE), size=nobs, prob=c(outliers_per_individual, 1 - outliers_per_individual), replace=TRUE)
          data_object$data[[o]]$data[, outl] = contaminated_data_matrices[[o]][, data_object$tt %in% data_object$data[[o]]$timepoints][, outl] 
          outlying_samples = c(outlying_samples, outl)
        } else {
          outlying_samples = c(outlying_samples, rep(FALSE, ncol(data_object$data[[o]]$data)))
        }
      }
      
    }
     
  }
  
  return (c(data_object, list(outlying_subjects = outlying_subjects, outlying_samples = outlying_samples)))
  
}


reformat_for_modeling = function(reduced_data) {

  create_matrix_from_data = function(red_data, i) {

    cbind(subject_id = i, s = rep(s, length(red_data$timepoints)), timepoint = rep(red_data$timepoints, each = length(s)), y = c(red_data$data))

  }

  s = reduced_data$ss
  n = length(reduced_data$data)
  
  reduced_list = vector("list", length=n)
  for (i in 1:n) {
    reduced_list[[i]] = create_matrix_from_data(reduced_data$data[[i]], i)
  }
  
  
  return (as.data.frame(do.call("rbind", reduced_list)))

}
