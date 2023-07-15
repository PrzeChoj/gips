#' A log of a posteriori that the covariance matrix is invariant under permutation
#'
#' More precisely, it is the logarithm of an unnormalized
#' posterior probability.
#' It is the goal function for optimization algorithms in [find_MAP()] function.
#' The `perm_proposal` that maximizes this function is
#' the Maximum A Posteriori (MAP) Estimator.
#'
#' It is calculated using
#' [formulas (33) and (27) from references](https://arxiv.org/abs/2004.03503).
#'
#' If `Inf` or `NaN` is reached, it produces a warning.
#'
#' @export
#'
#' @param g An object of a `gips` class.
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [arXiv link](https://arxiv.org/abs/2004.03503);
#' \doi{10.1214/22-AOS2174}
#'
#' @seealso
#' * [calculate_gamma_function()] - The function that calculates
#'     the value needed for `log_posteriori_of_gips()`.
#' * [get_structure_constants()] - The function that calculates
#'     the structure constants needed for `log_posteriori_of_gips()`.
#' * [find_MAP()] - The functions that optimizes
#'     the `log_posteriori_of_gips` function.
#' * [compare_posteriories_of_perms()] - Uses `log_posteriori_of_gips()`
#'     to compare a posteriori of two permutations.
#' * `vignette("Theory", package = "gips")` or its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html) - 
#'     A place to learn more about the math behind the `gips` package.
#'
#' @returns Returns a value of
#'     the logarithm of an unnormalized A Posteriori.
#'
#' @examples
#' # In the space with p = 2, there is only 2 permutations:
#' perm1 <- permutations::as.cycle("(1)(2)")
#' perm2 <- permutations::as.cycle("(1,2)")
#' S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
#' g1 <- gips(S1, 100, perm = perm1)
#' g2 <- gips(S1, 100, perm = perm2)
#' log_posteriori_of_gips(g1) # -134.1615, this is the MAP Estimator
#' log_posteriori_of_gips(g2) # -138.1695
#'
#' exp(log_posteriori_of_gips(g1) - log_posteriori_of_gips(g2)) # 55.0
#' # g1 is 55 times more likely than g2.
#' # This is the expected outcome because S[1,1] significantly differs from S[2,2].
#' 
#' compare_posteriories_of_perms(g1, g2)
#' # The same result, but presented in a more pleasant way
#'
#' # ========================================================================
#'
#' S2 <- matrix(c(1, 0.5, 0.5, 1.1), nrow = 2, byrow = TRUE)
#' g1 <- gips(S2, 100, perm = perm1)
#' g2 <- gips(S2, 100, perm = perm2)
#' log_posteriori_of_gips(g1) # -98.40984
#' log_posteriori_of_gips(g2) # -95.92039, this is the MAP Estimator
#'
#' exp(log_posteriori_of_gips(g2) - log_posteriori_of_gips(g1)) # 12.05
#' # g2 is 12 times more likely than g1.
#' # This is the expected outcome because S[1,1] is very close to S[2,2].
#' 
#' compare_posteriories_of_perms(g2, g1)
#' # The same result, but presented in a more pleasant way
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
#'
#' @noRd
runif_transposition <- function(perm_size) {
  sample(perm_size, 2, replace = FALSE)
}

#' Calculate log phi_part of log_posteriori_of_gips
#'
#' @param structure_constants An output of
#' `get_structure_constants(perm_proposal, perm_size)`.
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
  diagonalising_matrix <- prepare_orthogonal_matrix(perm_proposal)
  Dc_diagonalised <- t(diagonalising_matrix) %*% Dc %*% diagonalising_matrix
  DcUc_diagonalised <- t(diagonalising_matrix) %*% (Uc + Dc) %*% diagonalising_matrix

  # block part
  block_ends <- get_block_ends(structure_constants)
  Dc_block_log_dets <- calculate_log_determinants_of_block_matrices(
    Dc_diagonalised,
    block_ends
  )
  DcUc_block_log_dets <- calculate_log_determinants_of_block_matrices(
    DcUc_diagonalised,
    block_ends
  )
  Dc_exponent <- (delta - 2) / 2 + structure_constants[["dim_omega"]] /
    (structure_constants[["r"]] * structure_constants[["k"]])
  DcUc_exponent <- -(number_of_observations + delta - 2) / 2 - structure_constants[["dim_omega"]] /
    (structure_constants[["r"]] * structure_constants[["k"]])

  out <- sum(Dc_block_log_dets * Dc_exponent + DcUc_block_log_dets * DcUc_exponent)

  out
}

#' Calculate log of determinants of matrices from block decomposition
#'
#' Block decomposition 1 from paper
#'
#' @param diagonalized_matrix A middle matrix from decomposition 1.
#' @param block_ends The indices of last columns of block matrices.
#' Last element equals size of matrix.
#'
#' @returns A numeric vector.
#' @noRd
calculate_log_determinants_of_block_matrices <- function(diagonalised_matrix,
                                                         block_ends) {
  block_starts <- c(0, block_ends[-length(block_ends)] + 1)
  sapply(1:length(block_starts), function(i) {
    slice <- block_starts[i]:block_ends[i]
    block_matrix <- diagonalised_matrix[slice, slice, drop = FALSE]
    determinant(block_matrix, logarithm = TRUE)[["modulus"]]
  })
}

#' Compare the posteriori probabilities of 2 permutations
#'
#' Check which permutation is more likely and how much more likely.
#'
#' @param perm1,perm2 Permutations to compare.
#'     How many times `perm1` is more likely than `perm2`?
#'     Those can be provided as the `gips` objects,
#'     the `gips_perm` objects or anything that can be used as
#'     the `x` parameter in the [gips_perm()] function.
#'     They do not have to be of the same class.
#' @param S,number_of_observations,delta,D_matrix,was_mean_estimated
#'     The same parameters as in the [gips()] function.
#'     If at least one of `perm1` or `perm2` is a `gips` object,
#'     they are overwritten with those from `gips` object.
#' @param print_output A boolean.
#'     When `TRUE` (default), the computed value will be printed with
#'     additional text and returned invisibly. When `FALSE`,
#'     the computed value will be returned visibly.
#' @param digits Integer. Only used when `print_output = TRUE`.
#'     Number of digits after comma to print.
#'     Can be negative, can be `+Inf`. It is passed to `base::round()`.
#'
#' @returns The function `compare_posteriories_of_perms()` returns
#'     the value of how many times the `perm1` is more likely than `perm2`.
#'
#' @seealso
#' * [print.gips()] - The function that prints the posterior of
#'     the optimized `gips` object compared to the starting permutation.
#' * [summary.gips()] - The function that calculates the posterior of
#'     the optimized `gips` object compared to the starting permutation.
#' * [find_MAP()] - The function that finds the permutation that
#'     maximizes `log_posteriori_of_gips()`.
#' * [log_posteriori_of_gips()] - The function this 
#'     `compare_posteriories_of_perms()` calls underneath.
#'
#' @export
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 6
#' mu <- runif(6, -10, 10) # Assume we don't know the mean
#' sigma_matrix <- matrix(
#'   data = c(
#'     1.05, 0.8, 0.6, 0.4, 0.6, 0.8,
#'     0.8, 1.05, 0.8, 0.6, 0.4, 0.6,
#'     0.6, 0.8, 1.05, 0.8, 0.6, 0.4,
#'     0.4, 0.6, 0.8, 1.05, 0.8, 0.6,
#'     0.6, 0.4, 0.6, 0.8, 1.05, 0.8,
#'     0.8, 0.6, 0.4, 0.6, 0.8, 1.05
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' ) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- cov(Z) # Assume we have to estimate the mean
#'
#' g <- gips(S, number_of_observations)
#' g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
#'
#' compare_posteriories_of_perms(g_map, g, print_output = TRUE)
#' exp(compare_log_posteriories_of_perms(g_map, g, print_output = FALSE))
compare_posteriories_of_perms <- function(perm1, perm2 = "()", S = NULL,
                                          number_of_observations = NULL,
                                          delta = 3, D_matrix = NULL,
                                          was_mean_estimated = TRUE,
                                          print_output = TRUE,
                                          digits = 3) {
  compare_log <- compare_log_posteriories_of_perms(
    perm1, perm2, # if they are objects of `gips` class, they will be validated here
    S = S,
    number_of_observations = number_of_observations,
    delta = delta, D_matrix = D_matrix,
    was_mean_estimated = was_mean_estimated,
    print_output = FALSE
  )

  out <- exp(compare_log)
  
  if (!print_output) {
    return(out)
  }
  
  perm1_is_gips <- inherits(perm1, "gips")
  perm2_is_gips <- inherits(perm2, "gips")

  if (perm1_is_gips) {
    S <- attr(perm1, "S")
    perm1 <- perm1[[1]]
  }
  if (perm2_is_gips) {
    S <- attr(perm2, "S")
    perm2 <- perm2[[1]]
  }
  
  perm_size <- ncol(S)
  if (!inherits(perm1, "gips_perm")) {
    perm1 <- gips_perm(perm1, perm_size)
  }
  if (!inherits(perm2, "gips_perm")) {
    perm2 <- gips_perm(perm2, perm_size)
  }

  my_print_text <- paste0(
    "The permutation ", as.character.gips_perm(perm1), " is ", convert_log_diff_to_str(compare_log, digits),
    " times more likely than the ", as.character.gips_perm(perm2),
    " permutation."
  )
  if (out < 1) {
    my_print_text <- paste0(my_print_text, "\nThat means, the second permutation is more likely.")
  }
  cat(my_print_text)
  cat("\n")

  invisible(out)
}

#' @describeIn compare_posteriories_of_perms More stable,
#'     logarithmic version of `compare_posteriories_of_perms()`.
#'     The natural logarithm is used.
#'
#' @returns The function `compare_log_posteriories_of_perms()` returns
#'     the logarithm of how many times the `perm1` is more likely than `perm2`.
#'
#' @export
compare_log_posteriories_of_perms <- function(perm1, perm2 = "()", S = NULL,
                                              number_of_observations = NULL,
                                              delta = 3, D_matrix = NULL,
                                              was_mean_estimated = TRUE,
                                              print_output = TRUE,
                                              digits = 3) {
  perm1_is_gips <- inherits(perm1, "gips")
  perm2_is_gips <- inherits(perm2, "gips")
  
  if (perm1_is_gips) {
    validate_gips(perm1)
  }
  if (perm2_is_gips) {
    validate_gips(perm2)
  }
  
  if (perm1_is_gips && perm2_is_gips) {
    # Check the same parameters
    if (dim(attr(perm1, "S"))[1] != dim(attr(perm2, "S"))[1]) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `S` matrix! They cannot be compared!"),
                   class = "different_parameters")
    }
    if (any(abs((attr(perm1, "S") - attr(perm2, "S"))) > 0.00000001)) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `S` matrix! They cannot be compared!"),
                   class = "different_parameters")
    }
    if (attr(perm1, "number_of_observations") != attr(perm2, "number_of_observations")) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `number_of_observations`! They cannot be compared!"),
                   class = "different_parameters")
    }
    if (attr(perm1, "delta") != attr(perm2, "delta")) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `delta`! They cannot be compared!"),
                   class = "different_parameters")
    }
    if (dim(attr(perm1, "D_matrix"))[1] != dim(attr(perm2, "D_matrix"))[1]) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `D_matrix` matrix! They cannot be compared!"),
                   class = "different_parameters")
    }
    if (any(abs((attr(perm1, "D_matrix") - attr(perm2, "D_matrix"))) > 0.00000001)) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `D_matrix` matrix! They cannot be compared!"),
                   class = "different_parameters")
    }
    if (attr(perm1, "was_mean_estimated") != attr(perm2, "was_mean_estimated")) {
      rlang::abort(c("x" = "Give perm1 and perm2 are `gips` objects, but have different `was_mean_estimated`! They cannot be compared!"),
                   class = "different_parameters")
    }
  }
  
  if (perm1_is_gips) {
    S <- attr(perm1, "S")
    number_of_observations <- attr(perm1, "number_of_observations")
    delta <- attr(perm1, "delta")
    D_matrix <- attr(perm1, "D_matrix")
    was_mean_estimated <- attr(perm1, "was_mean_estimated")

    perm1 <- perm1[[1]]
  } else if (perm2_is_gips) {
    S <- attr(perm2, "S")
    number_of_observations <- attr(perm2, "number_of_observations")
    delta <- attr(perm2, "delta")
    D_matrix <- attr(perm2, "D_matrix")
    was_mean_estimated <- attr(perm2, "was_mean_estimated")

    perm2 <- perm2[[1]]
  }

  if (was_mean_estimated) {
    edited_number_of_observations <- number_of_observations - 1
  } else {
    edited_number_of_observations <- number_of_observations
  }

  perm_size <- ncol(S)
  if (!inherits(perm1, "gips_perm")) {
    perm1 <- gips_perm(perm1, perm_size)
  }
  if (!inherits(perm2, "gips_perm")) {
    perm2 <- gips_perm(perm2, perm_size)
  }

  check_correctness_of_arguments(S,
    edited_number_of_observations,
    max_iter = 5,
    start_perm = perm1, delta = delta, D_matrix = D_matrix,
    was_mean_estimated = was_mean_estimated, save_all_perms = TRUE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  )

  validate_gips_perm(perm1)
  validate_gips_perm(perm2)

  log_post1 <- log_posteriori_of_perm(
    perm1, S, edited_number_of_observations,
    delta, D_matrix
  )
  log_post2 <- log_posteriori_of_perm(
    perm2, S, edited_number_of_observations,
    delta, D_matrix
  )

  out <- log_post1 - log_post2
  
  if (!print_output) {
    return(out)
  }

  my_print_text <- paste0(
    "The permutation ", as.character.gips_perm(perm1), " is exp(", round(out, digits),
    ") times more likely than the ", as.character.gips_perm(perm2),
    " permutation."
  )
  if (out < 0) {
    my_print_text <- paste0(my_print_text, "\nThat means, the second permutation is more likely.")
  }
  cat(my_print_text)
  cat("\n")

  invisible(out)
}
