#' Prepare orthogonal matrix
#'
#' Calculate the orthogonal matrix `U_Gamma` for decomposition in
#' [Theorem 1 from references](https://arxiv.org/abs/2004.03503).
#'
#' Let \eqn{X} be a matrix invariant under the permutation `perm`.
#' Let the permutations cyclic group for it be called Gamma:
#' \eqn{\Gamma = <perm> = \{perm, perm^2, ...\}}.
#'
#' Then \eqn{U_\Gamma} is an orthogonal matrix that block-diagonalizes \eqn{X}.
#'
#' To be more precise, the matrix `t(U_Gamma) %*% X %*% U_Gamma` has a
#' block-diagonal structure, which is ensured by
#' [Theorem 1 from references](https://arxiv.org/abs/2004.03503).
#'
#' The formula for `U_Gamma` can be found in
#' [Theorem 6 from references](https://arxiv.org/abs/2004.03503).
#'
#' A nice example is demonstrated in the **Block Decomposition - \[1\], Theorem 1**
#' section of `vignette("Theory", package="gips")` or its
#' [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'
#' @param perm An object of a `gips_perm` or anything
#'     a `gips_perm()` can handle. It can also be of a `gips` class,
#'     but it will be interpreted as the underlying `gips_perm`.
#' @param perm_size Size of a permutation.
#'     Required if `perm` is neither `gips_perm` nor `gips`.
#' @param basis A matrix with basis vectors in COLUMNS; identity by default.
#' @returns A square matrix of size `perm_size` by `perm_size` with
#'     columns from vector elements \eqn{v_k^{(c)}} according to
#'     [Theorem 6 from references](https://arxiv.org/abs/2004.03503).
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [arXiv link](https://arxiv.org/abs/2004.03503);
#' \doi{10.1214/22-AOS2174}
#'
#' @seealso
#' * [project_matrix()] - A function used in examples
#'     to show the properties of `prepare_orthogonal_matrix()`.
#' * **Block Decomposition - \[1\], Theorem 1** section of
#'     `vignette("Theory", package = "gips")` or its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html) -
#'     A place to learn more about the math behind the `gips` package
#'     and see more examples of `prepare_orthogonal_matrix()`.
#'
#' @examples
#' gperm <- gips_perm("(1,2,3)(4,5)", 5)
#' U_Gamma <- prepare_orthogonal_matrix(gperm)
#'
#' number_of_observations <- 10
#' X <- matrix(rnorm(5 * number_of_observations), number_of_observations, 5)
#' S <- cov(X)
#' X <- project_matrix(S, perm = gperm) # this matrix in invariant under gperm
#'
#' block_decomposition <- t(U_Gamma) %*% X %*% U_Gamma
#' round(block_decomposition, 5) # the non-zeros only on diagonal and [1,2] and [2,1]
#' @export
prepare_orthogonal_matrix <- function(perm, perm_size = NULL, basis = NULL) {
  if (inherits(perm, "gips")) {
    validate_gips(perm)
    perm <- perm[[1]]
  }
  if (!inherits(perm, "gips_perm")) {
    perm <- gips_perm(perm, perm_size)
  }
  if (is.null(basis)) {
    return(prepare_orthogonal_matrix_default_fast_(perm))
  }
  v_object <- lapply(perm, function(subcycle) {
    get_v_matrix_for_subcycle(subcycle, basis)
  })
  arrange_v_object(v_object)
}

#' Prepare orthogonal matrix for the standard basis
#'
#' @param perm An object of class `gips_perm`.
#' @returns A square matrix of size `attr(perm, "size")`.
#' @noRd
prepare_orthogonal_matrix_default_fast_ <- function(perm) {
  p <- attr(perm, "size")
  if (p == 0) {
    return(matrix(numeric(), nrow = 0, ncol = 0))
  }

  v_matrix <- matrix(0, nrow = p, ncol = p)
  feature_matrix <- matrix(NA_real_, nrow = p, ncol = 3)
  column_start <- 1L

  for (cycle_index in seq_along(perm)) {
    subcycle <- perm[[cycle_index]]
    cycle_length <- length(subcycle)
    local_indices <- seq_len(cycle_length)
    k_s <- local_indices - 1
    cycle_columns <- column_start + local_indices - 1

    v_matrix[subcycle, cycle_columns[1]] <- 1 / sqrt(cycle_length)

    max_beta <- floor((cycle_length - 1) / 2)
    if (max_beta >= 1) {
      betas <- seq_len(max_beta)
      trygonometric_argument <- 2 * pi * outer(k_s, betas) / cycle_length
      v_matrix[subcycle, cycle_columns[2 * betas]] <-
        cos(trygonometric_argument) * sqrt(2 / cycle_length)
      v_matrix[subcycle, cycle_columns[2 * betas + 1]] <-
        sin(trygonometric_argument) * sqrt(2 / cycle_length)
    }

    if (cycle_length %% 2 == 0) {
      v_matrix[subcycle, cycle_columns[cycle_length]] <-
        cos(pi * k_s) / sqrt(cycle_length)
    }

    feature_matrix[cycle_columns, ] <- cbind(
      floor(local_indices / 2) / cycle_length,
      cycle_index,
      local_indices %% 2
    )
    column_start <- column_start + cycle_length
  }

  sorting_indices <- order(
    feature_matrix[, 1], feature_matrix[, 2],
    feature_matrix[, 3]
  )
  v_matrix[, sorting_indices]
}

#' Get V matrix
#'
#' @param subcycle An integer vector interpreted as cycle of a permutation.
#' @param basis An ortogonal matrix.
#' @returns A matrix p x length(subcycle).
#' Essentially a object v_k^c for one c value from paper, defined right before
#' Theorem 6.
#' @noRd
get_v_matrix_for_subcycle <- function(subcycle, basis) {
  cycle_length <- length(subcycle)
  k_s <- 1:cycle_length - 1
  chosen_basis_columns <- basis[, subcycle, drop = FALSE] # matrix p x curr_cycle_length

  first_element <- apply(chosen_basis_columns, 1, sum) /
    sqrt(cycle_length)
  v_matrix <- matrix(nrow = nrow(basis), ncol = cycle_length)
  v_matrix[, 1] <- first_element

  max_beta <- floor((cycle_length - 1) / 2)
  if (max_beta >= 1) {
    betas <- 1:max_beta
    trygonometric_argument <- 2 * pi * outer(k_s, betas) / cycle_length
    even_elements <-
      chosen_basis_columns %*% cos(trygonometric_argument) *
      sqrt(2 / cycle_length)
    odd_elements <-
      chosen_basis_columns %*% sin(trygonometric_argument) *
      sqrt(2 / cycle_length)

    v_matrix[, 2 * betas] <- even_elements
    v_matrix[, 2 * betas + 1] <- odd_elements
  }
  if (cycle_length %% 2 == 0) {
    last_element <- chosen_basis_columns %*% cos(pi * k_s) /
      sqrt(cycle_length)
    v_matrix[, cycle_length] <- last_element
  }
  v_matrix
}

#' Arrange V object
#'
#' @returns A matrix p x p with columns from V object elements, sorted according to
#'     Theorem 6.
#' @noRd
arrange_v_object <- function(v_object) {
  # we put v_k^v earlier than v_k'^c' if
  # i) [k/2]/p_c < [k'/2]/p_c' or
  # ii) [k/2]/p_c == [k'/2]/p_c' and c<c' or
  # iii) [k/2]/p_c == [k'/2]/p_c' and c == c' and k is even and k' is odd

  v_matrix <- do.call(cbind, v_object)

  features_list <- lapply(1:length(v_object), function(i) {
    cycle_length <- ncol(v_object[[i]])
    f_1 <- floor(1:cycle_length / 2) / cycle_length
    f_2 <- rep(i, cycle_length)
    f_3 <- 1:cycle_length %% 2
    matrix(c(f_1, f_2, f_3), ncol = 3)
  })

  feature_matrix <- do.call(rbind, features_list)
  sorting_indices <- order(
    feature_matrix[, 1], feature_matrix[, 2],
    feature_matrix[, 3]
  )
  v_matrix[, sorting_indices]
}
