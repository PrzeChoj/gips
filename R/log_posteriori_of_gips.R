#' A log of a posteriori that the covariance matrix is invariant under permutation
#'
#' More precisely, it is the logarithm of an unnormalized
#' posterior probability.
#' It is the goal function for optimization algorithms in [find_MAP()] function.
#' The `perm_proposal` that maximizes this function is
#' the Maximum A Posteriori (MAP) Estimator.
#'
#' It is calculated using
#' [formulas (33) and (27) from references](https://doi.org/10.1214/22-AOS2174).
#'
#' If `Inf` or `NaN` is reached, it produces a warning.
#'
#' @export
#'
#' @param g An object of a `gips_perm` class.
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kolodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [DOI: 10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)
#'
#' @seealso
#' * [calculate_gamma_function()] - The function that calculates the value
#'     needed for `log_posteriori_of_gips()`.
#' * [find_MAP()] - The functions that tries
#'     to optimize the `log_posteriori_of_gips` function.
#'
#' @returns Returns a value of
#'     the logarithm of an unnormalized A Posteriori.
#'
#' @examples
#' # In the space with p = 2, there is only 2 permutations:
#' perm1 <- permutations::as.cycle(permutations::as.word(c(1, 2))) # (1)(2)
#' perm2 <- permutations::as.cycle(permutations::as.word(c(2, 1))) # (1,2)
#' S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
#' g1 <- gips(S1, 100, perm = perm1)
#' g2 <- gips(S1, 100, perm = perm2)
#' log_posteriori_of_gips(g1) # -136.6, this is the MAP Estimator
#' log_posteriori_of_gips(g2) # -140.4
#'
#' exp(log_posteriori_of_gips(g1) - log_posteriori_of_gips(g2)) # 41.3
#' # g1 is over 40 times more likely than g2.
#' # This is the expected outcome because S[1,1] significantly differs from S[2,2].
#'
#' # ========================================================================
#'
#' S2 <- matrix(c(1, 0.5, 0.5, 1.1), nrow = 2, byrow = TRUE)
#' g1 <- gips(S2, 100, perm = perm1)
#' g2 <- gips(S2, 100, perm = perm2)
#' log_posteriori_of_gips(g1) # -99.5
#' log_posteriori_of_gips(g2) # -96.9, this is the MAP Estimator
#'
#' exp(log_posteriori_of_gips(g2) - log_posteriori_of_gips(g1)) # 12.7
#' # g2 is over 12 times more likely than g1.
#' # This is the expected outcome because S[1,1] is very close to S[2,2].
log_posteriori_of_gips <- function(g) {
  validate_gips(g)

  number_of_observations <- attr(g, "number_of_observations")
  was_mean_estimated <- attr(g, "was_mean_estimated")

  if (was_mean_estimated) {
    edited_number_of_observations <- number_of_observations - 1
  } else {
    edited_number_of_observations <- number_of_observations
  }

  log_posteriori_of_perm(
    perm_proposal = g[[1]], S = attr(g, "S"),
    number_of_observations = edited_number_of_observations,
    delta = attr(g, "delta"), D_matrix = attr(g, "D_matrix")
  )
}

log_posteriori_of_perm <- function(perm_proposal, S, number_of_observations,
                                   delta, D_matrix) {
  U <- S * number_of_observations # in the paper there is U everywhere instead of S, so it is easier to use U matrix in the code
  perm_size <- dim(S)[1]

  if (!inherits(perm_proposal, "gips_perm")) {
    perm_proposal <- gips_perm(perm_proposal, perm_size)
  }

  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size) # identity matrix
  }

  structure_constants <- get_structure_constants(perm_proposal)

  # Ac_part
  Ac <- sum(structure_constants[["r"]] * structure_constants[["k"]] * log(structure_constants[["k"]])) # (20)
  Ac_part <- (-number_of_observations / 2 * Ac)

  # G_part and phi_part
  G_part <- G_function(structure_constants, delta + number_of_observations) -
    G_function(structure_constants, delta)

  # phi_part
  phi_part <- calculate_phi_part(
    perm_proposal, number_of_observations, U,
    delta, D_matrix, structure_constants
  )

  out <- Ac_part + G_part + phi_part

  if (is.infinite(out)) {
    rlang::warn("The infinite value of a posteriori was produced.")
  }
  if (is.nan(out)) {
    rlang::warn("The NaN value of a posteriori was produced.")
  }

  out
}

#' Uniformly random transposition of perm_size elements
#'
#' @param perm_size A size from which take transpositions.
runif_transposition <- function(perm_size) {
  sample(perm_size, 2, replace = FALSE)
}

#' Calculate log phi_part of log_posteriori_of_gips
#'
#' @param structure_constants An output of `get_structure_constants(perm_proposal, perm_size)`.
#' Rest of params as in `log_posteriori_of_gips()`.
#'
#' @noRd
calculate_phi_part <- function(perm_proposal, number_of_observations, U,
                               delta, D_matrix, structure_constants) {

  # projection of matrices on perm_proposal
  equal_indices <- get_equal_indices_by_perm(perm_proposal)
  Dc <- project_matrix(D_matrix, perm_proposal,
    precomputed_equal_indices = equal_indices
  )
  Uc <- project_matrix(U, perm_proposal,
    precomputed_equal_indices = equal_indices
  )

  # divide by 2 - refer to newest version of the paper
  Dc <- Dc / 2
  Uc <- Uc / 2

  # diagonalization
  # TODO add basis argument? ISSUE#6
  diagonalising_matrix <- prepare_orthogonal_matrix(perm_proposal)
  Dc_diagonalised <- t(diagonalising_matrix) %*% Dc %*% diagonalising_matrix
  DcUc_diagonalised <- t(diagonalising_matrix) %*% (Uc + Dc) %*% diagonalising_matrix

  # block part
  block_ends <- cumsum(structure_constants[["r"]] * structure_constants[["d"]])
  Dc_block_dets <- calculate_determinants_of_block_matrices(
    Dc_diagonalised,
    block_ends
  )
  DcUc_block_dets <- calculate_determinants_of_block_matrices(
    DcUc_diagonalised,
    block_ends
  )
  Dc_exponent <- (delta - 2) / 2 + structure_constants[["dim_omega"]] /
    (structure_constants[["r"]] * structure_constants[["k"]])
  DcUc_exponent <- -(number_of_observations + delta - 2) / 2 - structure_constants[["dim_omega"]] /
    (structure_constants[["r"]] * structure_constants[["k"]])

  out <- sum(log(Dc_block_dets) * Dc_exponent + log(DcUc_block_dets) * DcUc_exponent)

  out
}

#' Calculate determinants of matrices from block decomposition
#'
#' Block decomposition 1 from paper
#'
#' @param diagonalized_matrix A middle matrix from decomposition 1.
#' @param block_ends The indices of last columns of block matrices.
#' Last element equals size of matrix.
#'
#' @returns A numeric vector.
#' @noRd
calculate_determinants_of_block_matrices <- function(diagonalised_matrix,
                                                     block_ends) {
  block_starts <- c(0, block_ends[-length(block_ends)] + 1)
  sapply(1:length(block_starts), function(i) {
    slice <- block_starts[i]:block_ends[i]
    block_matrix <- diagonalised_matrix[slice, slice, drop = FALSE]
    det(block_matrix)
  })
}
