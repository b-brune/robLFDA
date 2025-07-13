library(fda)
library(robustbase)
library(SpatialNP)
library(MASS)
library(kde1d)
library(mvtnorm)

mrct.rescale <- function(data,
                         data.cov,
                         alpha,
                         beta = 0,
                         H,
                         scaling.iterations = 5,
                         scaling.tolerance = 10^(-4),
                         operator = "id",
                         proj.dim = 3,
                         grid.diff
){

  k0 <- iter <- 1
  k1 <- Inf
  error <- (k0-k1)^2

  eigen.cov <- eigen(data.cov)

  data.sqrtcov <- eigen.cov$vectors %*% diag(abs(eigen.cov$values)^(1/2)) %*% t(eigen.cov$vectors)

  mrctoperator <- mrct.operator(data=data,
                                data.cov=data.cov,
                                data.sqrtcov=data.sqrtcov,
                                grid.diff=grid.diff,
                                k0=k0,
                                alpha=alpha,
                                beta=beta,
                                operator=operator,
                                proj.dim=proj.dim)
  std.data <- mrctoperator$std.data
  std.data.cov <- std.data %*% t(std.data)
  mhd <- diag(std.data.cov)

  while(iter <= scaling.iterations & error > scaling.tolerance){
    #print(iter)
    mrctoperator <- mrct.operator(data=data,
                                  data.cov=data.cov,
                                  data.sqrtcov=data.sqrtcov,
                                  grid.diff=grid.diff,
                                  k0=k0,
                                  alpha=alpha,
                                  beta=beta,
                                  operator=operator,
                                  proj.dim=proj.dim)
    std.data <- mrctoperator$std.data
    std.data.cov <- std.data %*% t(std.data)
    mhd <- diag(std.data.cov)

    if(operator == "id"){
      #browser()
      realizations <- mrct.distmhd(data.cov=data.cov,
                                   alpha=alpha,
                                   k=k0,
                                   iter = dim(data)[1])

      k1 <- stats::median(mhd)/stats::median(realizations)

      #reg <- MASS::rlm(sort(mhd)[1:H] ~ sort(realizations)[1:H])
      #k1 <- reg$coefficients[2]
      #print(paste0("Scaling parameter based in rlm: ", round(reg$coefficients[2],3)))
    }else if(operator == "idfone" | operator == "idproj" | operator == "projsq"){
      realizations <- mrct.cutoff(data.cov = data.cov,
                                  L = mrctoperator$L,
                                  L2 = mrctoperator$L2,
                                  alpha = alpha,
                                  beta = beta,
                                  k = k0,
                                  iter = dim(data)[1]
      )

      k1 <- stats::median(mhd) / stats::median(realizations)
      #k1 <- stats::quantile(mhd,0.75) / stats::quantile(realizations,0.75)

      # reg <- MASS::rlm(sort(mhd)[1:H] ~ sort(realizations)[1:H])
      # k1 <- reg$coefficients[2]
      #print(paste0("Scaling parameter based in rlm: ", round(reg$coefficients[2],3)))
    }else{
      realizations <- mrct.cutoff(data.cov = data.cov,
                                  L = mrctoperator$L,
                                  alpha = alpha,
                                  k = k0,
                                  iter = dim(data)[1]
      )

      k1 <- stats::median(mhd) / stats::median(realizations)

      # reg <- MASS::rlm(sort(mhd)[1:H] ~ sort(realizations)[1:H])
      # k1 <- reg$coefficients[2]
      #print(paste0("Scaling parameter based in rlm: ", round(reg$coefficients[2],3)))

    }

    error <- (k1-k0)^2
    iter <- iter + 1
    k0 <- k1
  }

  return(list("aMHD" = mhd,
              "std.data" = std.data,
              "hsubset" = order(mhd)[1:H],
              "scalingparameter" = ifelse(scaling.iterations == 0,1,k1))
  )
}

mrct.cutoff <- function(data.cov,
                        L,
                        L2 = matrix(0,nrow=dim(L)[1],ncol=dim(L)[2]),
                        alpha,
                        beta = 0,
                        k = 1,
                        iter = 2000,
                        seed = 123){

  set.seed(seed)

  eigen.cov <- eigen(data.cov)
  lambda <- eigen.cov$values
  psi <- eigen.cov$vectors
  #reg.inv.cov <- solve(data.cov + alpha / k * t(L) %*% L)
  # new for linear combination of two operators
  reg.inv.cov <- solve(data.cov + alpha / k * t(L) %*% L + beta / k * t(L2) %*% L2 )

  realizations <- c()
  for(i in 1:iter){
    eta <- rnorm(length(lambda))

    scaled.eigen <- rep(lambda * eta, rep.int(nrow(psi),ncol(psi))) * psi # in columns
    realizations[i] <- sum((reg.inv.cov %*% rowSums(scaled.eigen))^2)
  }

  return(realizations)
}

mrct.distmhd <- function(data.cov,
                         alpha,
                         k = 1,
                         iter = 2000,
                         seed =123){
  set.seed(seed)

  p <- dim(data.cov)[1]
  Q <- data.cov %*% solve(data.cov + alpha/k*diag(rep(1,p)))
  eigv.Q <- eigen(Q,symmetric = T)$values

  realizations <- c()
  for(i in 1:iter){
    beta <- stats::rchisq(p,1)
    realizations[i] <- sum(eigv.Q^2 * beta)
  }
  return(realizations)
}

mrct.kernel <- function(d,
                        sigma = 1,
                        l = 1,
                        method = "gaussian"){

  if(method == "gaussian"){
    return(sigma^2*exp(-d^2/(2*l^2)))
  }else if(method == "quadratic"){
    return(sigma*exp(-d^2/l))
  }else if(method == "linear"){
    return(sigma*exp(-abs(d)/l))
  }

}

mrct.operator <- function(data,
                          data.cov,
                          data.sqrtcov,
                          grid.diff,
                          k0,
                          alpha,
                          beta,
                          operator,
                          proj.dim = 3){

  if(operator == "id"){
    std.data <- t(solve(data.cov + alpha/k0*diag(rep(1,ncol(data)))) %*%  data.sqrtcov %*% t(data))
    P <- data.sqrtcov %*% solve(data.cov + alpha/k0*diag(rep(1,ncol(data)))) %*%  data.sqrtcov
    div <- diag(1,nrow=dim(data.cov)[1])
  }else if (operator == "fone"){

    # first differential operator
    # - derivative at endpoints are approximated by backward and forward differences
    # - der. at inner points of intervall are approximated by central differences
    # - results in a matrix R^{p \times p}
    # p <- dim(data)[2]
    # div <- matrix(0,ncol=p,nrow=p)
    # div[(1:(p - 1)) * (p + 1)] <- 1
    # div <- div - t(div)
    # div[1,1] <- div[p,(p-1)] <- -1
    # div[p,p] <- div[1,2] <- 1
    # div <- div / grid.length

    # only use backward differences (R^{(p-1) \times p})
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]
    #browser()
    std.data <- t(solve(data.cov + alpha/k0* t(div)%*%div) %*%  data.sqrtcov %*% t(data))
    P <- data.sqrtcov %*% solve(data.cov + alpha/k0* t(div)%*%div) %*%  data.sqrtcov
  }else if(operator == "ftwo"){
    # second derivative
    # - derivative at endpoints are approximated by backward and forward differences
    # - der. at inner points of intervall are approximated by central differences
    # - results in a matrix R^{p \times p}
    # p <- dim(data)[2]
    # div <- matrix(0,ncol=p,nrow=p)
    # div[(1:(p - 1)) * (p + 1)] <- 1
    # div <- div + t(div)
    # diag(div) <- -2
    # div[1,1:3] <- div[p,(p-2):p] <- c(1,-2,1)
    # div <- div / grid.length^2

    # only using central differences (R^{(p-2) \times p})
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div + t(div)
    diag(div) <- -2
    div <- div[-c(1,p),] / grid.diff^2

    std.data <- t(solve(data.cov + alpha/k0* t(div)%*%div) %*%  data.sqrtcov %*% t(data))
    P <- data.sqrtcov %*% solve(data.cov + alpha/k0* t(div)%*%div) %*%  data.sqrtcov
  }else if(operator == "fthree"){
    # third derivative (tbc.)

  }else if(operator =="idfone"){
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    std.data <- t(solve(data.cov + alpha/k0*diag(rep(1,ncol(data))) + beta/k0 * t(div) %*% div) %*%  data.sqrtcov %*% t(data))
    P <- data.sqrtcov %*% solve(data.cov + alpha/k0*diag(rep(1,ncol(data))) + beta/k0 * t(div) %*% div) %*%  data.sqrtcov
  }else if(operator == "idproj"){
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    eigen.data <- eigen(data.cov)
    W <- eigen.data$vectors[,1:5]
    div <- div %*% (diag(rep(1,ncol(data))) - W %*% t(W))

    std.data <- t(solve(data.cov + alpha/k0*diag(rep(1,ncol(data))) + beta/k0 * t(div) %*% div) %*%  data.sqrtcov %*% t(data))
    P <- data.sqrtcov %*% solve(data.cov + alpha/k0*diag(rep(1,ncol(data))) + beta/k0 * t(div) %*% div) %*%  data.sqrtcov
    # std.data <- t(solve(data.cov + alpha/k0* t(div)%*%div) %*%  data.sqrtcov %*% t(data))
    # P <- data.sqrtcov %*% solve(data.cov + alpha/k0* t(div)%*%div) %*%  data.sqrtcov
  }else if(operator == "projsq"){
    #browser()
    p <- dim(data)[2]

    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    eigen.data <- eigen(data.cov)
    W <- eigen.data$vectors[,1:proj.dim]
    div1 <- div %*% (diag(rep(1,ncol(data))) - W %*% t(W))
    div2 <- diag(rep(1,ncol(data))) - W %*% t(W)

    std.data <- t(solve(data.cov + alpha/k0* t(div1) %*% div1 + beta/k0 * t(div2) %*% div2) %*%  data.sqrtcov %*% t(data))
    P <- data.sqrtcov %*% solve(data.cov + alpha/k0*t(div1) %*% div1 + beta/k0 * t(div2) %*% div2) %*%  data.sqrtcov
  }

  if(operator == "idfone" | operator == "idproj"){
    return(list("std.data" = std.data,
                "P" = P,
                "L" = diag(rep(1,ncol(data))),
                "L2" = div)
    )
  }else if(operator == "projsq"){
    return(list("std.data" = std.data,
                "P" = P,
                "L" = div1,
                "L2" = div2)
    )
  }else{
    return(list("std.data" = std.data,
                "P" = P,
                "L" = div)
    )
  }

}

mrct.operatordata <- function(data,
                              data.cov,
                              data.sqrtcov,
                              grid.diff,
                              k0,
                              alpha,
                              beta,
                              operator,
                              proj.dim = 3){

  if(operator == "id"){
    operator.data <- diag(rep(1,ncol(data))) %*% t(data)
  }else if (operator == "fone"){

    # only use backward differences
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    operator.data <- div %*% t(data)
  }else if(operator == "ftwo"){
    # p <- dim(data)[2]
    # div <- matrix(0,ncol=p,nrow=p)
    # div[(1:(p - 1)) * (p + 1)] <- 1
    # div <- div + t(div)
    # diag(div) <- -2
    # div[1,1:3] <- div[p,(p-2):p] <- c(1,-2,1)
    # div <- div / grid.length^2

    # only using central differences (R^{(p-2) \times p})
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div + t(div)
    diag(div) <- -2
    div <- div[-c(1,p),] / grid.diff^2

    operator.data <- div %*% t(data)
  }else if(operator =="idfone"){
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    div <- alpha*diag(rep(1,ncol(data)))[-p,] + beta *div

    operator.data <- div %*% t(data)
  }else if(operator =="idproj"){
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    eigen.cov <- eigen(data.cov)
    W <- eigen.cov$vectors[,1:5]
    div <- div %*% (diag(rep(1,ncol(data))) - W %*% t(W))

    # add Identity to regularization
    div <- alpha*diag(rep(1,ncol(data)))[-p,] + beta *div

    operator.data <- div %*% t(data)
  }else if(operator =="projsq"){
    p <- dim(data)[2]
    div <- matrix(0,ncol=p,nrow=p)
    diag(div) <- -1
    div[(1:(p - 1)) * (p + 1)] <- 1
    div <- div / grid.diff
    div <- div[-p,]

    eigen.cov <- eigen(data.cov)
    W <- eigen.cov$vectors[,1:proj.dim]
    div1 <- div %*% (diag(rep(1,ncol(data))) - W %*% t(W))
    div2 <- diag(rep(1,ncol(data))) - W %*% t(W)

    div <- alpha*div1 + beta*div2[-p,]

    operator.data <- div %*% t(data)
  }

  return(t(operator.data))
}

#' @export
mrct <- function(data,
                 h = 0.75,
                 alpha = 0.01,
                 beta = 0.01,
                 initializations = 5,
                 subset.iteration = 10,
                 seed = 123,
                 scaling.iterations = 10,
                 scaling.tolerance = 10^(-4),
                 criterion = "sum",
                 operator = "id",
                 proj.dim = 3,
                 grid.diff = diff(seq(0,1,length.out=dim(data)[2]))[1],
                 fixed.initialization = F
){

  if(h < 0 | h > 1){
    stop("The value of `h` cannot be negative or above 1!")
  }
  if(h < 0.5){
    stop("`h` = ",h," was selected. For values below 0.5 the results might break down if too many outliers are present. Please consider a value between 0.5 and 1.")
  }

  N <- nrow(data)
  p <- ncol(data)

  # In this setting only one initialization is considered
  if(fixed.initialization){
    data.med <- data - matrix(rep(robustbase::colMedians(data),N), ncol = p, byrow = T)
    data.med <- robustbase::rowMedians(abs(data.med))

    initializations <- 1
  }

  # Allocation of data frames
  objval <- c()
  hsubsets <- matrix(NA,
                     ncol=floor(h*N),
                     nrow=initializations)
  scalingparameters <- c()

  #####################################################################################
  ### Part 1: Start with >=1 initial subsets and iterate until convergence for each ###
  #####################################################################################

  for(i in 1:initializations){

    set.seed(seed+i)
    if(fixed.initialization){
      subset.initial <- order(data.med)[1:floor(h*N)]
    }else{
      subset.initial <- sample(1:N,floor(h*N))
    }

    data.centered <- data - matrix(rep(robustbase::colMedians(data),N), ncol = p, byrow = T)
    data.cov <- t(data.centered[subset.initial,]) %*% data.centered[subset.initial,] / length(subset.initial)

    k <- 1

    subset.old <- subset.initial
    subset.new <- Inf

    while(k <= subset.iteration
          & setequal(subset.old,subset.new) == F
    ){

      if(k >= 2){subset.old <- subset.new}

      # Calculate new h-subset
      #browser()
      #print(k)
      tmp <- mrct.rescale(data = data.centered,
                          data.cov = data.cov,
                          H = floor(N*h),
                          alpha = alpha,
                          beta = beta,
                          scaling.iterations = scaling.iterations,
                          operator = operator,
                          grid.diff = grid.diff,
                          proj.dim = proj.dim
      )

      subset.new <- tmp$hsubset

      data.centered <- data - matrix(rep(colMeans(data[subset.new,]),N),byrow=T,ncol=p)
      data.cov <- stats::cov(data.centered[subset.new,])

      k <- k + 1
    }

    mhd <- tmp$aMHD / tmp$scalingparameter

    if(criterion == "sum"){
      objval[i] <- sum(diag(data.cov))
    }else if(criterion == "cluster"){
      res <- stats::kmeans(mhd,centers = 2)
      objval[i] <- res$tot.withinss/res$betweenss
    }

    hsubsets[i,] <- subset.new
    scalingparameters[i] <- tmp$scalingparameter
  }

  ###################################
  ### Part 2: Select final subset ###
  ###################################
  # Take subset with smallest obj value
  subset.optimal <- hsubsets[order(objval)[1],]
  scalingparameter.optimal <- scalingparameters[order(objval)[1]]

  # "Extract" optimal covariance and scaled a-MHD
  data.centered <- data - matrix(rep(colMeans(data[subset.optimal,]),N),byrow=T,ncol=p)
  data.cov <- stats::cov(data.centered[subset.optimal,])
  eigen.cov <- eigen(data.cov)
  data.sqrtcov <- eigen.cov$vectors %*% diag(abs(eigen.cov$values)^(1/2)) %*% t(eigen.cov$vectors)

  mrctoperator <- mrct.operator(data=data.centered,
                                data.cov=data.cov,
                                data.sqrtcov=data.sqrtcov,
                                grid.diff=grid.diff,
                                k0=scalingparameter.optimal,
                                alpha=alpha,
                                operator=operator,
                                beta = beta,
                                proj.dim=proj.dim)

  if(operator == "id"){
    std.data <- mrctoperator$std.data
    mhd <- diag(std.data %*% t(std.data)) / scalingparameter.optimal
    dist <- mrct.distmhd(data.cov=data.cov,
                         alpha = alpha,
                         k = scalingparameter.optimal)

    quant <- stats::quantile(dist,0.975) # quantile
  }else if(operator == "idfone" | operator == "idproj" | operator == "projsq"){
    std.data <- mrctoperator$std.data
    mhd <- diag(std.data %*% t(std.data)) / scalingparameter.optimal
    L <- mrctoperator$L

    dist <- mrct.cutoff(data.cov = data.cov,
                        L = L,
                        L2 = mrctoperator$L2,
                        alpha = alpha,
                        beta = beta,
                        k = scalingparameter.optimal)

    quant <- stats::quantile(dist,0.95) # quantile
  }else{
    std.data <- mrctoperator$std.data
    mhd <- diag(std.data %*% t(std.data)) / scalingparameter.optimal
    L <- mrctoperator$L

    dist <- mrct.cutoff(data.cov = data.cov,
                        L = L,
                        alpha = alpha,
                        k = scalingparameter.optimal)

    quant <- stats::quantile(dist,0.95) # quantile
  }

  ##############################
  ### Part 3: Weighting step ###
  ##############################
  non.outliers <- as.numeric(which(mhd <= quant))

  if(length(non.outliers) == 0){
    data.centered.weighted <- data - matrix(rep(colMeans(data),N),byrow=T,ncol=p)
    data.cov.weighted <- stats::cov(data.centered.weighted)
  }else if(length(non.outliers) == 1){
    data.centered.weighted <- data - matrix(rep(data[non.outliers,],N),byrow=T,ncol=p)
    #print(1)
    data.cov.weighted <- stats::sd(data.centered.weighted[non.outliers,])^2

    theoretical.outliers <- as.numeric(which(mhd >  quant))

    output <- list("theoretical" = theoretical.outliers,
                   "theoretical.w" = theoretical.outliers,
                   "aMHD" = mhd,
                   "aMHD.w" = mhd,
                   "quant" = quant,
                   "quant.w" = quant,
                   "k" = scalingparameter.optimal,
                   "k.w" = scalingparameter.optimal,
                   "optimal.subset" = subset.optimal,
                   "subsets" = hsubsets,
                   "objval" = objval,
                   "std.data" = std.data,
                   "centered.data" = data.centered,
                   "operator.data" = mrct.operatordata(data=data.centered,
                                                       data.cov=data.cov,
                                                       data.sqrtcov=data.sqrtcov,
                                                       grid.diff=grid.diff,
                                                       k0=scalingparameter.optimal,
                                                       alpha=alpha,
                                                       operator=operator,
                                                       beta = beta,
                                                       proj.dim=proj.dim)
    )

    return(output)
  }else{
    data.centered.weighted <- data - matrix(rep(colMeans(data[non.outliers,]),N),byrow=T,ncol=p)
    data.cov.weighted <- stats::cov(data.centered.weighted[non.outliers,])
  }

  temp <- mrct.rescale(data = data.centered.weighted,
                       data.cov = data.cov.weighted,
                       alpha = alpha,
                       H = floor(h*N),
                       scaling.tolerance = scaling.tolerance,
                       scaling.iterations = scaling.iterations,
                       operator = operator,
                       grid.diff = grid.diff,
                       beta = beta,
                       proj.dim=proj.dim
  )

  mhd.weighted <- temp$aMHD/temp$scalingparameter
  if(operator == "id"){
    dist.weighted <- mrct.distmhd(data.cov = data.cov.weighted,
                                  alpha = alpha,
                                  k = temp$scalingparameter)

    quant.weighted <- stats::quantile(dist.weighted,0.975) # quantile
  }else if(operator == "idfone" | operator == "idproj" | operator == "projsq"){

    dist <- mrct.cutoff(data.cov = data.cov.weighted,
                        L = mrctoperator$L,
                        L2 = mrctoperator$L2,
                        alpha = alpha,
                        beta = beta,
                        k = temp$scalingparameter)

    quant.weighted <- stats::quantile(dist,0.95) # quantile
  }else{
    L <- mrctoperator$L

    dist <- mrct.cutoff(data.cov = data.cov.weighted,
                        L = L,
                        alpha = alpha,
                        k = temp$scalingparameter)

    quant.weighted <- stats::quantile(dist,0.99) # quantile
  }

  ###########################################################
  ### Part 4: Determine outliers (with/without weighting) ###
  ###########################################################
  theoretical.outliers.weighted <- as.numeric(which(mhd.weighted >  quant.weighted))
  theoretical.outliers <- as.numeric(which(mhd >  quant))

  output <- list("theoretical" = theoretical.outliers,
                 "theoretical.w" = theoretical.outliers.weighted,
                 "aMHD" = mhd,
                 "aMHD.w" = mhd.weighted,
                 "quant" = quant,
                 "quant.w" = quant.weighted,
                 "k" = scalingparameter.optimal,
                 "k.w" = temp$scalingparameter,
                 "optimal.subset" = subset.optimal,
                 "subsets" = hsubsets,
                 "objval" = objval,
                 "std.data" = std.data,
                 "centered.data" = data.centered,
                 "operator.data" = mrct.operatordata(data=data.centered,
                                                     data.cov=data.cov,
                                                     data.sqrtcov=data.sqrtcov,
                                                     grid.diff=grid.diff,
                                                     k0=scalingparameter.optimal,
                                                     alpha=alpha,
                                                     operator=operator,
                                                     beta = beta,
                                                     proj.dim=proj.dim)
  )

  # calculate k*covariance of the optimal subset!!!

  return(output)
}
