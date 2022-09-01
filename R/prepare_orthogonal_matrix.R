#' Prepare orthogonal matrix
#'
#' Calculate orthogonal matrix U_Gamma for decomposition in
#' [Theorem 1 from references](https://doi.org/10.1214/22-AOS2174).
#'
#' Given X - invariant under the permutation `perm`. Call Gamma the
#' permutations cyclic group \eqn{<perm> = \{perm, perm^2, ...\}}.
#'
#' Then, U_Gamma is such an orthogonal matrix that X is *"pretty"* in it.
#'
#' To be more precise, the matrix `t(U_Gamma) %*% X %*% U_Gamma` has
#' a lot of zeros (see examples).
#'
#' To calculate it, [Theorem 6 from references](https://doi.org/10.1214/22-AOS2174) is used.
#'
#' @param perm Object of a `gips_perm` or a `permutations::cycle` class.
#' @param perm_size Size of permutation.
#'     Required if `perm` is of a `permutations::cycle` class.
#' @param basis A matrix with basis vectors in COLUMNS. Identity by default.
#' @returns A matrix `perm_size` x `perm_size` with columns from
#'     vector elements \eqn{\{v_k^{(c)}\}} according to
#'     [Theorem 6 from references](https://doi.org/10.1214/22-AOS2174).
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kolodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [DOI: 10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)
#'
#' @seealso
#' * [project_matrix()] - A function used in examples
#'     to show the properties of `prepare_orthogonal_matrix()`.
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
#' t(U_Gamma) %*% X %*% U_Gamma # the non-zeros only on diagonal and [1,2] and [2,1]
#' @export
prepare_orthogonal_matrix <- function(perm, perm_size = NULL, basis = NULL) {
  if (!inherits(perm, "gips_perm")) {
    perm <- gips_perm(perm, perm_size)
  }
  if (is.null(basis)) {
    basis <- diag(nrow = attr(perm, "size"))
  }
  v_object <- lapply(perm, function(subcycle) {
    get_v_matrix_for_subcycle(subcycle, basis)
  })
  arrange_v_object(v_object)
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
