#' Shift vector
#'
#' Move k elements from the start of a vector to its end.
#'
#' @noRd
shift_vector <- function(v, k) {
  if (k == 0) {
    return(v)
  }
  c(v[-(1:k)], v[1:k])
}

#' Rearrange vector
#'
#' Move elements from the start of a vector to its end, so that the minimal
#' element will be first.
#'
#' @examples
#' v <- c(5, 3, 2, 1, 4)
#' rearranged <- rearrange_vector(v)
#' all(rearranged == c(1, 4, 5, 3, 2)) # TRUE
#'
#' @noRd
rearrange_vector <- function(v) {
  shift_vector(v, which.min(v) - 1)
}

#' Is the numeric value representing the whole number
#'
#' This code is copied from [base::integer()] example.
#'
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (!is.numeric(x)) {
    return(rep(FALSE, length(x)))
  }
  abs(x - round(x)) < tol
}

#' Is matrix symmetric
#'
#' We did not use the `matrixcalc::is.positive.semi.definite` function, because
#' here we have no checks(because they were done before) and the `tol`
#' argument there is taken relative, but we want the absolute.
#' Also, the default of the `tolerance` argument here is `1e-06`,
#' which is the same as in the `MASS::mvrnorm` function.
#'
#' Watch out that this function does NOT checks weather
#' the `matrix_of_interest` is indeed a matrix.
#'
#' @noRd
is.positive.semi.definite.matrix <- function(matrix_of_interest, tolerance = 1e-06) {
  eigenvalues <- eigen(
    matrix_of_interest,
    symmetric = TRUE,
    only.values = TRUE
  )[["values"]]

  return(all(eigenvalues >= -tolerance * abs(eigenvalues[1]))) # 1st is the biggest eigenvalue
}

wrong_argument_abort <- function(i, x = "") {
  rlang::abort(c("There was a problem identified with provided argument",
    "i" = i,
    "x" = x
  ))
}

#' Used primarily in tests
#'
#' Maybe later in vignette?
#'
#' @noRd
to_perm <- function(v) permutations::as.cycle(permutations::as.word(v))
