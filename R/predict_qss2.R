#' @export

predict.qss2 <- function (object, newdata, ...)
{
  x <- object[, 1]
  y <- object[, 2]
  z <- object[, 3]
  tri.area <- function(v) {
    0.5 * ((v[2, 1] - v[1, 1]) * (v[3, 2] - v[1, 2]) - (v[3,
                                                          1] - v[1, 1]) * (v[2, 2] - v[1, 2]))
  }
  barycentric <- function(v) {
    b <- rep(0, 3)
    Area <- tri.area(v[1:3, ])
    b[1] <- tri.area(v[c(4, 2, 3), ])/Area
    b[2] <- tri.area(v[c(1, 4, 3), ])/Area
    b[3] <- tri.area(v[c(1, 2, 4), ])/Area
    if (any(b < 0) || any(b > 1))
      stop("barycentric snafu")
    b
  }
  if (is.list(newdata)) {
    fnames <- (dimnames(object)[[2]])[1:2]
    if (all(!is.na(match(fnames, names(newdata))))) {
      newx <- newdata[[fnames[1]]]
      newy <- newdata[[fnames[2]]]
    }
    else (stop("qss object and newdata frame names conflict"))
  }
  else if (is.matrix(newdata))
    if (ncol(newdata) == 2) {
      newx <- newdata[, 1]
      newy <- newdata[, 2]
    }
  else (stop("newdata matrix must have 2 columns"))
  trinew <- interp::tri.mesh(x, y)
  tri <- interp::triSht2tri(trinew)
  if (!all(interp::in.convex.hull(trinew, newx, newy, strict=FALSE)))
    stop("some newdata points outside convex hull")
  p <- length(x)
  m <- length(newx)
  V <- matrix(0, m, 3)
  B <- matrix(0, m, 3)
  for (i in 1:m) {
    Tmp <- interp::tri.find(trinew, newx[i], newy[i])
    V[i, ] <- c(Tmp[[1]], Tmp[[2]], Tmp[[3]])
    v <- rbind(cbind(x[V[i, ]], y[V[i, ]]), c(newx[i], newy[i]))
    B[i, ] <- barycentric(v)
  }
  ra <- c(t(B))
  ja <- as.integer(c(t(V)))
  ia <- as.integer(3 * (0:m) + 1)
  D <- new("matrix.csr", ra = ra, ja = ja, ia = ia, dimension = c(m,
                                                                  p))
  list(x = newx, y = newy, z = c(D %*% z), D = D[, -1])
}


#' @export
predict.rqss <- function (object, newdata, interval = "none", level = 0.95, ...)
{
  ff <- object$fake.formula
  Terms <- delete.response(terms(object$formula, "qss"))
  Names <- all.vars(parse(text = ff))
  if (any(!(Names %in% names(newdata))))
    stop("newdata doesn't include some model variables")
  nd <- eval(model.frame(ff, data = newdata), parent.frame())
  qssterms <- attr(Terms, "specials")$qss
  if (length(qssterms)) {
    tmp <- quantreg::untangle.specials(Terms, "qss")
    dropv <- tmp$terms
    m <- length(dropv)
    if (length(dropv))
      PLTerms <- Terms[-dropv]
    attr(PLTerms, "specials") <- tmp$vars
  }
  else {
    PLTerms <- Terms
    m <- 0
  }
  if (requireNamespace("MatrixModels") && requireNamespace("Matrix"))
    X <- as(MatrixModels::model.Matrix(PLTerms, data = nd,
                                       contrasts = contrasts, sparse = TRUE), "matrix.csr")
  else X <- model.matrix(PLTerms, data = nd)
  p <- ncol(X)
  y <- X %*% object$coef[1:p]
  X <- SparseM:::as.matrix.csr(X)
  if (m > 0) {
    for (i in 1:m) {
      qss <- object$qss[[i]]
      names <- all.vars(Terms[dropv[i]])
      names <- names[names %in% Names]
      dimnames(qss$xyz)[[2]] <- c(names, "zfit")
      newd <- nd[names]
      if (ncol(qss$xyz) == 3) {
        g <- predict.qss2(qss$xyz, newdata = newd, ...)
        y <- y + g$z
        if (interval == "confidence")
          X <- cbind(X, g$D)
      }
      else if (ncol(qss$xyz) == 2) {
        g <- quantreg:::predict.qss(qss, newdata = newd, ...)
        y <- y + g$y
        if (interval == "confidence")
          X <- cbind(X, g$D)
      }
      else stop("invalid fitted qss object")
    }
  }
  if (interval == "confidence") {
    v <- sqrt(diag(X %*% summary(object, cov = TRUE)$V %*%
                     t(X)))
    calpha <- qnorm(1 - (1 - level)/2)
    y <- cbind(y, y - v * calpha, y + v * calpha)
    dimnames(y)[[2]] <- c("yhat", "ylower", "yupper")
  }
  y
}
