test_that("log_posteriori_of_gips returns proper values", {
  # The value of log_posteriori_of_gips on matrix should the same as the projection of the matrix
  # and S2 == pi_c(S1)
  c_perm <- permutations::as.cycle(permutations::as.word(c(2, 1)))
  id_perm <- permutations::id

  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE) / 100
  S2 <- matrix(c(1.5, 0.5, 0.5, 1.5), nrow = 2, byrow = TRUE) / 100

  D_matrix <- diag(nrow = 2)

  expect_equal(
    log_posteriori_of_gips(gips(S1, 100,
      perm = c_perm, D_matrix = D_matrix,
      was_mean_estimated = FALSE
    )),
    log_posteriori_of_gips(gips(S2, 100,
      perm = c_perm, D_matrix = D_matrix,
      was_mean_estimated = FALSE
    ))
  )

  # Those values were calculated by hand:

  expect_equal(
    exp(log_posteriori_of_gips(gips(S1 * 2, 100,
      perm = c_perm,
      D_matrix = D_matrix * 2, was_mean_estimated = FALSE
    ))),
    6^(-103 / 2) * gamma(103 / 2) * gamma(103 / 2) / (pi / 4)
  )
  expect_equal(
    exp(log_posteriori_of_gips(gips(S1 * 2, 100,
      perm = id_perm,
      D_matrix = D_matrix * 2, was_mean_estimated = FALSE
    ))),
    (23 / 4)^(-52) * gamma(52) * gamma(51.5) * sqrt(2 * pi) / (pi / sqrt(2))
  )
  expect_equal(
    exp(log_posteriori_of_gips(gips(S2 * 2, 100,
      perm = id_perm,
      D_matrix = D_matrix * 2, was_mean_estimated = FALSE
    ))),
    6^(-52) * gamma(52) * gamma(51.5) * sqrt(2 * pi) / (pi / sqrt(2))
  )
})

test_that("log_posteriori_of_perm returns proper values", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE) / 100

  expect_silent(log_posteriori_of_perm(
    perm_proposal = "(1,2)", S = S1,
    number_of_observations = 100,
    3, D_matrix = NULL
  ))
})

test_that("log_posteriori_of_gips has the desired property", {
  # Example from the paper chapter 5
  # This test is randomized.
  # It is mathematically possible the Z variables will be drawn such that
  # the test fails. However, it is highly unlikely.
  # See ISSUE#9 for discussion.

  id_perm <- permutations::id

  p <- 10
  n <- 20

  mu <- numeric(p)
  sigma_matrix <- matrix(numeric(p * p), nrow = p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma_matrix[i, j] <- 1 - min(abs(i - j), p - abs(i - j)) / p
    }
    sigma_matrix[i, i] <- 1 + 1 / p
  }

  Z <- MASS::mvrnorm(n, mu = mu, Sigma = sigma_matrix)
  S <- t(Z) %*% Z / n

  actual_permutation <- permutations::as.cycle(permutations::as.word(c(2:p, 1)))

  actual_permutation_function_value <- log_posteriori_of_gips(gips(S, n, perm = actual_permutation, was_mean_estimated = FALSE))
  another_permutation_function_value <- log_posteriori_of_gips(gips(S, n, perm = id_perm, was_mean_estimated = FALSE))

  # We want the posteriori to have a bigger value for the real permutation than for the another
  expect_gt(
    actual_permutation_function_value,
    another_permutation_function_value
  )
})

test_that("calculate phi_part returns proper values", {
  gips_example_perm <- gips_perm(example_perm, 6)
  D_matrix <- matrix(c(
    10, 1, 1, 2, 2, 3,
    1, 10, 1, 2, 2, 3,
    1, 1, 10, 2, 2, 3,
    2, 2, 2, 12, 4, 5,
    2, 2, 2, 4, 12, 5,
    3, 3, 3, 5, 5, 14
  ), byrow = TRUE, ncol = 6) * 2
  delta <- 3
  n <- 1
  U_matrix <- diag(6) * 2
  structure_constants <- get_structure_constants(gips_example_perm)
  expected_phi_part <- log(2206^(-3) * 100^(-3 / 2) * 9^(-2)) - log(1680^(-5 / 2) * 81^(-1) * 8^(-3 / 2))

  expect_equal(
    calculate_phi_part(
      gips_example_perm, n, U_matrix, delta, D_matrix,
      structure_constants
    ),
    expected_phi_part
  )
})

test_that("calculate_block_determinants returns proper values", {
  block_diagonal_matrix <- matrix(c(
    2, 1, 0, 0, 0,
    1, 2, 0, 0, 0,
    0, 0, 2, 0, 0,
    0, 0, 0, 5, 0,
    0, 0, 0, 0, 0.2
  ), byrow = TRUE, ncol = 5)

  expect_equal(
    calculate_determinants_of_block_matrices(
      block_diagonal_matrix,
      c(2, 3, 5)
    ),
    c(3, 2, 1)
  )
})

test_that("compare_posteriories_of_perms properly calculates", {
  gips_example_perm <- gips_perm(example_perm, 6)
  gips_example_perm2 <- gips_perm(example_perm2, 6)
  gips_id <- gips_perm("()", 6)

  g <- gips(matrix_invariant_by_example_perm, 14, perm = gips_example_perm)
  g2 <- gips(matrix_invariant_by_example_perm, 14, perm = gips_example_perm2)

  expect_equal(
    compare_log_posteriories_of_perms(g, example_perm, print_output = FALSE),
    0
  )
  expect_equal(
    compare_posteriories_of_perms(g, example_perm, print_output = FALSE),
    1
  )
  expect_equal(
    compare_posteriories_of_perms(g, gips_example_perm, print_output = FALSE),
    1
  )

  expect_equal(
    compare_posteriories_of_perms(g2, g, print_output = FALSE),
    94914.44395167663
  )
  expect_equal(
    compare_posteriories_of_perms(gips_example_perm2, g, print_output = FALSE),
    94914.44395167663
  )

  expect_equal(
    expect_output(
      compare_posteriories_of_perms(g, g2, print_output = TRUE),
      "times more likely than the \\(1,2,3,4,5\\) permutation\\.\\nThat means, the second permutation is more likely\\."
    ),
    1 / 94914.44395167663
  )
  expect_equal(
    expect_output(
      compare_log_posteriories_of_perms(g, g2, print_output = TRUE),
      "times more likely than the \\(1,2,3,4,5\\) permutation\\.\\nThat means, the second permutation is more likely\\."
    ),
    -log(94914.44395167663)
  )

  expect_equal(
    compare_posteriories_of_perms(g, print_output = FALSE),
    compare_posteriories_of_perms(g, gips_id, print_output = FALSE)
  )
})
