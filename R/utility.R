#' Plots a matrix using geom_tile()
#'
#' @param cc the matrix that should be plotted
#' @param index the index that should be used for the x- and y-axis of the plot
#'
#' @returns A ggplot

plot_matrix = function(cc, index=NULL) {
  if (!is.null(index)) {
    colnames(cc) <- rownames(cc) <- index
    labs = "t"
  } else {
    labs = ""
  }

  cc2 = reshape2::melt(cc)

  require(ggplot2)

  return(
    cc2 |>
      ggplot2::ggplot(ggplot2::aes(Var1, Var2, fill=value)) +
      ggplot2::geom_tile() +
      ggplot2::xlab(labs) +
      ggplot2::ylab(labs) +
      ggplot2::scale_fill_gradient2(low = "blue", mid="white", high="red", midpoint=0)
  )
}


#' Calculate the RMSE
#' @param a true values
#' @param b estimated values of same length as `a`
#'
#' @returns The RMSE between \code{a} and \code{b}
#' @export

rmse <- function(a, b) {
  if (length(a) != length(b)) {
    stop("`a` and `b` need to be the same length.")
  }

  return(
    sqrt(mean((a - b)^2))
  )
}

#' Calculate the mean absolute deviation
#' @param a true values
#' @param b estimated values of same length as `a`
#'
#' @returns The mean absolute deviation between \code{a} and \code{b}
#'
#' @export

mean_abs_dev <- function(a, b) {

  if (length(a) != length(b)) {
    stop("`a` and `b` need to be the same length.")
  }

  return(
    mean(abs(a - b))
  )
}


