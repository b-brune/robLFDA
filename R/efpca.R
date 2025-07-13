#' Elliptical-FPCA main function
#'
#' @param X number of values
#'
#' @export
efpca <- function(
    X,
    opt.h.cov,
    cutoff_outliers=2.5,
    rho.param=1e-3,
    k = 3,
    s = 3,
    ncov=50,
    max.kappa=1e3,
    prediction_grid = NULL,
    mean_estimator = function(x) mean(x, na.rm=TRUE)) {

  # X is a list of 2 named elements "x", "pp"
  # which contain the lists of observations and times for each item
  # X$x[[i]] and X$pp[[i]] are the data for the i-th individual

  ## matrixx() creates all pairs of observations and corresponding timepoints
  ma <- matrixx(X)

  # Compute the estimated cov function
  #
  rlog::log_info("Calculating the estimated covariance function.")
  cov.fun2 <- cov.fun.hat2(X=X, h=opt.h.cov, ma=ma, ncov=ncov, trace=FALSE)

  # Perform the smoothing step
  yy <- as.vector(cov.fun2$G)
  xx <- cov.fun2$grid

  smoothing_data = dplyr::tibble(cov_fun = yy, grid1 = xx[, 1], grid2 = xx[, 2])

  smoothed_cov_model = mgcv::gam(cov_fun ~ s(grid1, grid2), data=smoothing_data, family="gaussian")


  # Predict from the smoothed covariance function and symmetrize it.

  # Get predictions on the right grid
  if (is.null(prediction_grid)) {
    tt = unique(cov.fun2$grid[, 1])
  } else {
    tt = prediction_grid
  }

  tmp <- predict(smoothed_cov_model, newdata=tidyr::expand_grid(grid1 = tt, grid2 = tt), na.action=na.pass)
  cov.fun2$G <- matrix(tmp, length(tt), length(tt))
  cov.fun <- ( cov.fun2$G + t(cov.fun2$G) ) / 2 # This is the covariance function


  # predicted scores, fitted values

  la1 <- eigen(cov.fun)$values[1]
  rho.param <- la1 / (max.kappa - 1)

  # Flagging the outliers:

  outlyingness = lapply(seq_along(X$x), function(i) {
    calculate_outlyingness(x=X$x[[i]], t_true = X$pp[[i]], grid = tt, sigma=cov.fun, mean_estimator=mean_estimator)
  })

  flagged = lapply(outlyingness, ">", cutoff_outliers)

  result = list(
      cov.fun=cov.fun,
      rho.param = rho.param,
      flagged = flagged,
      outlyingness = outlyingness,
      grid=tt
    )

  class(result) = "efpca"

  return(
    result
  )
}


#' @export
fitted.efpca <- function(object, X, max_components = 20, pve=0.9, ...) {

  result = predict_score_functions(
    X = X,
    cov.fun = object$cov.fun,
    tt = object$grid,
    s = max_components,
    rho = object$rho.param,
    pve = pve,
    flagged_observations=object$flagged
  )

  return (result)
}





#' Estimates the covariance matrix
#' @param estimate_both_directions added by BB, do we really need to perform the smoothing step twice or
#'                                 can we just assume the matrix is symmetric :)
#' @export
cov.fun.hat2 <- function(X, h, ma, ncov=50, trace=FALSE) {
  # this function uses the diagonal
  if(trace) print("Computing cov function")

  if(missing(ma)) ma <- matrixx(X=X) # calculate the centered pairs

  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  betahat <- rep(0, np)
  sigmahat <- rep(0, ncov)

  ## Calculate the diagonal
  for(j in 1:ncov) {
    sigmahat[j] <- gtthat(X=X, t0=tt[j], h=h) # solves equation (7) for each t_0
  }

  ## Calculate the off-diagonal elements
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    betahat[j] <- gtshat(X=X, t0=t0, s0=s0, h=h, matx=ma, eps=1e-6) # estimate the off-diagonal covriances

  }
  G <- betahat <- matrix(betahat, ncov, ncov)
  for(i in 1:ncov) {
    for(j in 1:ncov) {
      G[i,j] <- betahat[i,j] * sigmahat[i]^2
    }
  }

  G <- ( G + t(G) ) / 2

  if(trace) print("Done computing cov function")

  return(list(G=G, grid=pps))
}


#' This calculates gamma(t_0, t_0) as in equation (7) in the paper
#' @export
gtthat <- function(X, t0, h, b=.5, cc=1.54764, initial.sc, max.it=300, eps=1e-10) {
  i <- 0
  err <- 1 + eps
  t <- unlist(X$pp)
  x <- unlist(X$x)

  if (missing(initial.sc)) { sc <- stats::mad(x) } else { sc <- initial.sc }


  while( ( (i <- i+1) < max.it ) && (err > eps) ) {

    kerns <- k.epan((t - t0)/h)
    sc2 <- sqrt(sc^2 * sum(kerns * rho(x / sc, cc)) / (b * sum(kerns)))
    err <- abs(sc2 / sc - 1)
    sc <- sc2

  }

  return(sc)
}

#' @export
gtshat <- function(X, t0, s0, h, matx, cc=3.443689, eps=1e-6){ # 3.443689 4.685065 1.345

  n <- length(X$x)
  err <- 1+eps
  j <- 0
  if (missing(matx)) matx <- matrixx(X=X)

  M <- matx$m
  MT <- matx$mt

  w2 <- k.epan((MT[, 1] - t0) / h)
  w3 <- k.epan((MT[, 2] - s0) / h)
  we <- w2*w3
  M <- M[ we > 0, ]
  we <- we[ we > 0]

  if (length(we) == 0) return(NA)

  B <- stats::median( M[,2] / M[,1 ])
  sigma.hat <- stats::mad( M[,2] - B * M[,1])

  while ( ( (j <- j+1) < 1000 ) && (err > eps) ){

    w1 <- tukey.weight( ( (M[,2]-B*M[,1]) / sigma.hat ) / cc )
    w <- w1 * we
    B2 <- sum(w * M[, 1] * M[, 2]) / sum(w * (M[ ,1]^2))
    err <- abs(B2/B - 1)

    if( is.na(err) || is.nan(err) ) return(NA)

    B <- B2
  }

  return(B)
}




#' Predict smoothed trajectories from the estimated covariance model
#' @param X list of lists of length two with measurements and timepoints
#' @param cov.fun smoothed covariance matrix (length(tt) x length(tt) matrix)
#' @param tt timepoints to predict on
#' @param s maximum number of components
#' @param rho number of measurements
#' @param pve pve to decide on `k`, the number of components
predict_score_functions <- function(
    X,
    cov.fun,
    tt,
    flagged_observations,
    s = 20,
    rho = 0,
    pve = 0.9) {


  # compute eigenfunctions of the smoothed covariance surface
  eg <- eigen(cov.fun)

  rlog::log_info("Selecting number of components based on pve")

  # kick out negative eigenvalues
  ev = eg$values[eg$values > 1e-12]

  ef = eg$vectors[, eg$values > 1e-12]

  prop = cumsum(ev) / sum(ev)
  k = which(prop == min(prop[prop >= pve]))


  s1 <- max(k, min(s, length(ev))) # number of components to use for score calculation / covariance matrix of X

  # Use a maximum of s1 = max(k,s) eigenvalues / eigenfunctions
  lam <- ev[1:s1]
  ef <- ef[, 1:s1]

  # Standardize eigenfunctions
  # L2.normal.mesh calculates the norm and the line after that one normalizes the
  # eigenfunctions so that they have norm 1

  # Here we create the mesh on which we evaluate the covariance grid after smoothing;
  # this should contain all potential timepoints we're interested in (tt)

  eval_mesh = seq(min(tt), max(tt), length.out=ncol(cov.fun))

  normas <- apply(ef, 2, L2.norma.mesh, mesh=eval_mesh) # rep(1, max(k,s))
  efn <- scale(ef, center=FALSE, scale=normas)

  # make them into functions
  ff <- vector('list', s1)

  for(i in 1:s1) {
    # Calls approxfun to interpolate the grid and get an interpolation function on
    # eval_mesh
    ff[[i]] <- stats::approxfun(eval_mesh, efn[, i], method='linear')
  }


  # Compute predicted scores
  n <- length(X$x) # number of subjects
  xis <- matrix(NA, n, k) # matrix of eigenfunction values (n x k) -> one value per eigenfunction and subject
  pr <- matrix(NA, n, length(tt)) # predictions -> interpolated trajectories

  for (i in 1:n) {
    # Filter out the clean observations:
    # If too many observations (i.e. > 50%) were flagged, we simply keep everything
    okay = (is.na(flagged_observations[[i]]) | !flagged_observations[[i]])

    if (sum(okay) < (length(X$x[[i]]) / 2)) {
      okay = rep(TRUE, length(X$x[[i]]))
      flagged_observations[[i]] = !okay
    }

    xx = X$x[[i]][okay]
    ti = X$pp[[i]][okay]

    xic <- xx # centered X_i

    phis <- matrix(0, length(ti), s1) # matrix to get the measurements
    phis_full <- matrix(0, length(tt), s1)

    ## Evaluate the functions resulting from the interpolation
    for(j in 1:s1) {
      phis[,j] <- ff[[j]](ti)
      phis_full[,j] <- ff[[j]](tt)
    }


    # Sigma_Y
    siy <- phis[, 1:s1, drop=FALSE] %*% diag( lam[1:s1] ) %*% t(phis[, 1:s1, drop=FALSE])
    # Sigma_Y^{-1} (X - mu)

    rhs <- as.vector( solve(siy + rho * diag(length(ti)), xic ) )

    # scores  = \phi' \lambdas Sigma_Y^{-1} (X - mu)
    if (k > 1) {
      xis[i, ] <- t( phis[, 1:k, drop=FALSE] %*% diag( lam[1:k] ) ) %*% rhs
    } else {
      xis[i, ] <- t( phis[, 1, drop=FALSE] * lam[1] ) %*% rhs
    }
    # \hat{X - mu} = \sum phi_j \xi_j
    pr[i, ] <- phis_full[, 1:k, drop=FALSE] %*% as.vector( xis[i,] )
  }


  return(
    list(
      xis=xis,
      pred=pr,
      flagged_observations = ifelse(is.na(unlist(flagged_observations)), FALSE, unlist(flagged_observations))
    )
  )
}
