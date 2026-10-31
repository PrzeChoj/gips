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

#' Check positive semidefiniteness of a previously validated matrix
#'
#' The tolerance is relative to the largest absolute matrix entry.
#' This helper does not check the matrix type, dimensions, finiteness or symmetry.
#'
#' @noRd
is.positive.semi.definite.matrix <- function(
    matrix_of_interest, tolerance = sqrt(.Machine$double.eps)) {
  eigenvalues <- eigen(
    matrix_of_interest,
    symmetric = TRUE,
    only.values = TRUE
  )[["values"]]

  return(all(eigenvalues >= -tolerance * max(abs(matrix_of_interest))))
}

#' Check positive definiteness of a previously validated symmetric matrix
#'
#' Cholesky tests factorization without imposing a minimum relative eigenvalue.
#' Matrix type, dimensions, finiteness and symmetry must be checked by callers.
#'
#' @noRd
is.positive.definite.matrix <- function(matrix_of_interest) {
  tryCatch(
    {
      chol(matrix_of_interest, pivot = FALSE)
      TRUE
    },
    error = function(e) FALSE
  )
}

is_scalar_identity_matrix <- function(matrix_of_interest) {
  if (nrow(matrix_of_interest) != ncol(matrix_of_interest)) {
    return(FALSE)
  }

  if (nrow(matrix_of_interest) <= 1) {
    return(TRUE)
  }

  # Only exact scalar matrices may use the simplified posterior calculation.
  if (matrix_of_interest[1, 2] != 0) {
    return(FALSE)
  }

  # fast diagonal check
  diagonals <- diag(matrix_of_interest)
  d <- diagonals[1]
  if (any(diagonals != d)) {
    return(FALSE)
  }

  all(matrix_of_interest == diag(d, nrow(matrix_of_interest)))
}

wrong_argument_abort <- function(i, x = "") {
  rlang::abort(c("There was a problem identified with provided argument",
    "i" = i,
    "x" = x
  ))
}

#' Used primarily in tests
#'
#' @noRd
to_perm <- function(v) permutations::as.cycle(permutations::as.word(v))

change_log_probabilities_unnorm_to_probabilities <- function(log_probabilities_unnorm) {
  log_probabilities_unnorm <- log_probabilities_unnorm - max(log_probabilities_unnorm)
  probabilities_unnorm <- exp(log_probabilities_unnorm)
  probabilities <- probabilities_unnorm / sum(probabilities_unnorm)

  probabilities
}

pretty_plot_matrix <- function(S, title = "") {
  plot(gips(S, 1, was_mean_estimated = FALSE), type = "heatmap") +
    ggplot2::labs(title = title, x = "", y = "", fill = "value")
}

pretty_plot_block_matrix <- function(S, perm, title = "") {
  # S is the original cov estimator - not diagonalized!
  plot(gips(S, 1, perm = perm, was_mean_estimated = FALSE), type = "block_heatmap") +
    ggplot2::labs(title = title, x = "", y = "", fill = "value")
}

get_block_ends <- function(structure_constants) {
  cumsum(structure_constants[["r"]] * structure_constants[["d"]])
}

OEIS_A051625 <- c(1, 2, 5, 17, 67, 362, 2039, 14170, 109694, 976412, 8921002, 101134244, 1104940280, 13914013024, 191754490412, 2824047042632, 41304021782824, 708492417746000, 11629404776897384, 222093818836736752, 4351196253952132832, 88481681599705382144)
OEIS_A000142 <- c(1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000, 51090942171709440000, 1124000727777607680000)
