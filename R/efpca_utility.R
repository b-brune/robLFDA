#' Epanechnikov kernel
k.epan <- function(x) {
  #' Evaluates the Epanechnikov kernel for a given x
  #' @param x the value
  a <- 0.75 * (1 - x^2)
  k.epan <- a * (abs(x) <= 1)
  return(k.epan)
}


## ************************************************************************** ##
## SPIKE DETECTION
## ************************************************************************** ##

#' Calculate the outlyingness for a pair of observations given a covariance function
get_factor <- function(x, t, sigma) {
  #' Get the standardized residuals
  #' @param x the two X-values X_ij and X_il
  #' @param t the indices of the corresponding timepoints t_ij and t_il in the grid
  #'
  #' @returns The standardized sqrt residual

  if (t[1] == t[2]) return(0)

  beta = sigma[t[1], t[2]] / sigma[t[2], t[2]]
  diff = x[1] - beta * x[2]
  std = sigma[t[1], t[1]] - (sigma[t[1], t[2]]^2 / sigma[t[2], t[2]])

  if (std < 1e-4 | is.na(std)) return (NA)

  return(sqrt(diff^2 / std))
}


#' Calculate the outlyingness for the observations from one individual
calculate_outlyingness <- function(
    x,
    t_true,
    grid,
    sigma,
    mean_estimator=function(x) mean(x, na.rm=TRUE) # estimates the average for mean calculation
) {
  #' Calculates the outlyingness matrix for one set of observations
  #'
  #' @param x observation vector
  #' @param t_true timestamps
  #' @param grid the grid we're extrapolating to
  #' @param sigma covariance matrix


  t = grid %in% t_true
  sigma_subset = sigma[t, t]

  res = matrix(0, nrow=length(x), ncol=length(x))

  for (i in 1:length(x)) {
    for (j in 1:length(x)) {
      index = c(i, j)

      xx = x[index]
      tt_true = t_true[index]

      res[i, j] = get_factor(xx, index, sigma_subset)
    }
  }

  mean_stand_residuals = apply(res, 1, mean_estimator)

  return(mean_stand_residuals)
}


#########
#########

integral <- function(efe,mesh){
  # we need that the first and last element of mesh
  # are the extreme points of the interval where we integrate
  m <- length(mesh)
  t0 <- min(mesh)
  tm1 <- max(mesh)
  primero<- (mesh[1]-t0)*efe[1]
  ultimo<-(tm1-mesh[m])*efe[m]
  sumo<-primero+ultimo
  menos1<-m-1
  for (i in 1:menos1)
  {
    a1<-(efe[i]+efe[i+1])/2
    amplitud<- mesh[i+1]-mesh[i]
    sumo<-sumo+amplitud*a1
  }
  return(sumo)
}



#################################
# Calcula el producto interno entre dos datos (como listas)  sobre una grilla de puntos MESH en (0,1) ordenados
##########################

#' @export
L2.dot.product.mesh <- function(dato1, dato2, mesh)
{
  return(integral(dato1*dato2,mesh))
}



#################################
# Calcula la norma de un dato funcional sobre una grilla de puntos MESH en (0,1) ordenados
# (OJO ES LA NORMA no el cuadrado de la norma!!!)
##########################

#' @export
L2.norma.mesh <- function(dato, mesh)
{
  return(sqrt(L2.dot.product.mesh(dato,dato,mesh)))
}


#' Operator squared norm distance
#'
#' This function computes the squared "operator
#' norm distance" between its arguments. This is
#' closely related to the squared Frobenius norm
#' of the difference of the matrices.
#'
#' @param gammax the matrix representation of the first operator
#' @param gamma0 the matrix representation of the second operator
#'
#' @return the squared Frobenius norm of the difference of
#' the matrices, divided by the number of
#' elements.
#'
#' @export
norma <- function(gammax, gamma0){
  norm <- NA
  kx <- dim(gammax)[2]
  k0 <- dim(gamma0)[2]
  if(kx!=k0) { print("cuidado!!") }
  if(kx==k0) {
    gamaresta<-gammax-gamma0
    #  gamarestavec<-c(lowerTriangle(gamaresta,diag=TRUE),upperTriangle(gamaresta,diag=FALSE))
    # norm<-sum(gamarestavec*gamarestavec) /(kx*kx)
    norm<- sum(gamaresta^2)/(kx*kx)
  }
  return( norm )
}


#' All possible subsets
#'
#' This function returns all possible subsets of a given set.
#'
#' This function returns all possible subsets of a given set.
#'
#' @param n the size of the set of which subsets will be computed
#' @param k the size of the subsets
#' @param set an optional vector to be taken as the set
#'
#' @return A matrix with subsets listed in its rows
#'
#' @export
subsets <- function(n, k, set = 1:n) {
  if(k <= 0) NULL else if(k >= n) set
  else rbind(cbind(set[1], Recall(n - 1, k - 1, set[-1])), Recall(n - 1, k, set[-1]))
}



#' Trimmed mean of squares
#'
#' This function returns the mean of the smallest
#' (1-alpha)% squares of its arguments.
#'
#' @param x numeric vector
#' @param alpha real number between 0 and 1
#'
#' @return the mean of the smallest (1-alpha)% squared values in \code{x}
#'
#' @export
tm <- function(x, alpha) {
  n <- length(x)
  return( mean( (sort(x^2, na.last=NA))[1:(n - floor(alpha*n))], na.rm=TRUE ) )
}



#' Trimmed mean
#'
#' This function returns the mean of the smallest
#' (1-alpha)% of its arguments.
#'
#' @param x numeric vector
#' @param alpha real number between 0 and 1
#'
#' @return the mean of the smallest (1-alpha)% values in \code{x}
#'
#' @export
tmns <- function(x, alpha) {

  n <- length(x)
  return( mean( (sort(x, na.last=NA))[1:(n - floor(alpha*n))], na.rm=TRUE ) )
}


#' @export
matrixx <- function(X) {
  # build all possible pairs Y_{ij}, Y_{il}, j != l
  # and the corresponding times t_{ij}, t_{i,l}, j != l
  # They are returned in "m" and "mt" below

  n <- length(X$x)
  M <- MT <- NULL
  for (i in 1:n){

    comb <- subsets(length(X$x[[i]]),2)

    if (class(comb)[1] != 'matrix') comb <- matrix(comb, byrow=TRUE, ncol=2)

    Maux <- matrix(X$x[[i]][comb], ncol=2)
    MTaux <- matrix(X$pp[[i]][comb], ncol=2)
    M <- rbind(M, Maux, cbind(Maux[, 2], Maux[, 1]))
    MT <- rbind(MT, MTaux, cbind(MTaux[, 2], MTaux[, 1]))
  }
  return(list(m = M, mt = MT))
}

#' @export
rho <- function(a, cc=3) {
  tmp <- 1-(1-(a/cc)^2)^3
  tmp[abs(a/cc)>1] <- 1
  return(tmp)
}

#' @export
tukey.weight <- function(a) {
  tmp <- 6*(1-a^2)^2
  tmp[ abs(a) > 1 ] <- 0
  return(tmp)
}



