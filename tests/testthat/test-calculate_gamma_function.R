test_that("calculate_gamma_omega returns proper value", {
  delta <- 3
  gips_c <- gips_perm(permutations::as.cycle(permutations::as.word(c(2, 1))), 2)
  gips_cprim <- gips_perm(permutations::id, 2)

  # perm = (1)(2)
  structure_constants <- get_structure_constants(gips_cprim)

  dim_omega_i <- structure_constants[["dim_omega"]][1]
  r_i <- structure_constants[["r"]][1]
  d_i <- structure_constants[["d"]][1]
  k_i <- structure_constants[["k"]][1]
  expect_equal(
    c(dim_omega_i, r_i, d_i, k_i),
    c(3, 2, 1, 1)
  )

  lambda <- 1 / 2 * k_i * (delta - 2) + dim_omega_i / r_i
  expect_equal(lambda, 2)

  expect_equal(
    exp(calculate_gamma_omega(
      lambda = lambda,
      dim_omega_i = dim_omega_i,
      r_i = r_i,
      d_i = d_i
    )),
    pi / sqrt(2)
  )

  # perm = (1,2)
  structure_constants <- get_structure_constants(gips_c)

  dim_omega_i <- structure_constants[["dim_omega"]][1]
  r_i <- structure_constants[["r"]][1]
  d_i <- structure_constants[["d"]][1]
  k_i <- structure_constants[["k"]][1]
  expect_equal(
    c(dim_omega_i, r_i, d_i, k_i),
    c(1, 1, 1, 1)
  )

  lambda <- 1 / 2 * k_i * (delta - 2) + dim_omega_i / r_i
  expect_equal(lambda, 1.5)

  expect_equal(
    exp(calculate_gamma_omega(
      lambda = lambda,
      dim_omega_i = dim_omega_i,
      r_i = r_i,
      d_i = d_i
    )),
    sqrt(pi / 4)
  ) # gamma(1.5)
})

test_that("calculate_gamma_omega gives warning and infinity on divergent integral", {
  expect_warning(out <- calculate_gamma_omega(0.4, 3, 2, 1))
  expect_true(is.infinite(out))
})

test_that("calculate_gamma_function returns Inf and warning when integral diverges", {
  id_perm <- gips_perm(permutations::id, 2)
  expect_warning(
    out <- calculate_gamma_function(id_perm, 0.5),
    "Gamma integral is divergent for the given permutation and lambda value."
  )
  expect_true(is.infinite(out))
})

test_that("calculate_G_part returns the expected value when L is 1", {
  delta <- 3
  number_of_observations <- 10
  gips_cprim <- gips_perm(permutations::id, 2)

  structure_constants <- get_structure_constants(gips_cprim)
  lambda_prior <- 1 / 2 * structure_constants[["k"]][1] * (delta - 2) +
    structure_constants[["dim_omega"]][1] / structure_constants[["r"]][1]
  lambda_posterior <- lambda_prior +
    structure_constants[["k"]][1] * number_of_observations / 2

  expect_equal(
    calculate_G_part(structure_constants, delta, number_of_observations),
    calculate_gamma_omega(
      lambda_posterior,
      structure_constants[["dim_omega"]][1],
      structure_constants[["r"]][1],
      structure_constants[["d"]][1]
    ) -
      calculate_gamma_omega(
        lambda_prior,
        structure_constants[["dim_omega"]][1],
        structure_constants[["r"]][1],
        structure_constants[["d"]][1]
      )
  )
})

test_that("calculate_G_part matches precomputed reference values", {
  perms <- list(
    gips_perm("()", 5),
    gips_perm("(1,2,3)(4,5)", 6),
    gips_perm("(1,2,3,4,5,6,7)", 7)
  )
  delta <- 3
  number_of_observations <- 20
  
  res <- c(91.7179250705395, 114.157676814219, 152.553199409633)
  
  for (i in seq_along(perms)) {
    perm <- perms[[i]]
    structure_constants <- get_structure_constants(perm)
    
    expect_equal(
      calculate_G_part(structure_constants, delta, number_of_observations),
      res[i]
    )
  }
})

test_that("G part example from documentation", {
  perm_size <- 6
  perm <- permutations::as.cycle(permutations::as.word(c(2, 3, 1, 5, 4, 6)))
  gips_perm <- gips_perm(perm, 6)
  structure_constants <- get_structure_constants(gips_perm)
  delta <- 3
  number_of_observations <- 10

  G_prior <- sum(vapply(seq_len(structure_constants[["L"]]), function(i) {
    lambda <- structure_constants[["k"]][i] * (delta - 2) / 2 +
      structure_constants[["dim_omega"]][i] / structure_constants[["r"]][i]

    calculate_gamma_omega(
      lambda,
      structure_constants[["dim_omega"]][i],
      structure_constants[["r"]][i],
      structure_constants[["d"]][i]
    )
  }, numeric(1)))
  G_posterior <- sum(vapply(seq_len(structure_constants[["L"]]), function(i) {
    lambda <- structure_constants[["k"]][i] *
      (delta + number_of_observations - 2) / 2 +
      structure_constants[["dim_omega"]][i] / structure_constants[["r"]][i]

    calculate_gamma_omega(
      lambda,
      structure_constants[["dim_omega"]][i],
      structure_constants[["r"]][i],
      structure_constants[["d"]][i]
    )
  }, numeric(1)))

  expect_equal(
    calculate_G_part(structure_constants, delta, number_of_observations),
    G_posterior - G_prior
  )
})

test_that("calculate_G_part stops when the internal convergence check fails", {
  # These structure constants and delta are not allowed for valid `gips` inputs.
  # This only tests the defensive branch for corrupted/internal misuse.
  divergent_structure_constants <- list(
    L = 1,
    k = 1,
    r = 2,
    d = 1,
    dim_omega = 3
  )

  expect_error(
    calculate_G_part(divergent_structure_constants, delta = 0, number_of_observations = 1),
    "Gamma integral divergence detected"
  )
})

test_that("calculate_gamma_function has desired properties", {
  gips_cprim <- gips_perm(permutations::id, 2)

  expect_gt(
    calculate_gamma_function(gips_cprim, 0.50000001),
    calculate_gamma_function(gips_cprim, 0.5001)
  )

  expect_gt(
    calculate_gamma_function(gips_cprim, 0.500000000001),
    calculate_gamma_function(gips_cprim, 0.50000001)
  )

  expect_warning(calculate_gamma_function(gips_cprim, 0.5))
})

test_that("calculate_gamma_function can get a gips object", {
  my_perm <- gips_perm("(1,2,3)", 3)
  expect_silent(gamma1 <- calculate_gamma_function(my_perm, 0.5001))
  expect_silent(gamma2 <- calculate_gamma_function(gips(diag(3), 13, perm = my_perm), 0.5001))
  expect_equal(gamma1, gamma2)
})
