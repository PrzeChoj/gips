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

test_that("when L is 1, G_function returns the same value as calculate_gamma_omega", {
  delta <- 3
  gips_cprim <- gips_perm(permutations::id, 2)

  structure_constants <- get_structure_constants(gips_cprim)
  lambda <- 1 / 2 * structure_constants[["k"]][1] * (delta - 2) + structure_constants[["dim_omega"]][1] / structure_constants[["r"]][1]

  expect_equal(
    G_function(structure_constants, delta),
    calculate_gamma_omega(
      lambda,
      structure_constants[["dim_omega"]][1],
      structure_constants[["r"]][1],
      structure_constants[["d"]][1]
    )
  )
})

test_that("G_function example from documentation", {
  perm_size <- 6
  perm <- permutations::as.cycle(permutations::as.word(c(2, 3, 1, 5, 4, 6)))
  gips_perm <- gips_perm(perm, 6)
  structure_constants <- get_structure_constants(gips_perm)

  expect_true(abs(exp(G_function(structure_constants, 3)) - 16.443561) < 0.0001)
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
