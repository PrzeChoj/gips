test_that("get_structure_constants for examples from paper", {
  gips_example_perm <- gips_perm(example_perm, 6)
  gips_example_perm2 <- gips_perm(example_perm2, 5)
  example_structure_constants <- list(
    "r" = c(3, 1, 1),
    "d" = c(1, 2, 1),
    "k" = c(1, 2, 1),
    "L" = 3,
    "dim_omega" = c(6, 1, 1)
  )
  example_structure_constants2 <- list(
    "r" = c(1, 1, 1),
    "d" = c(1, 2, 2),
    "k" = c(1, 2, 2),
    "L" = 3,
    "dim_omega" = c(1, 1, 1)
  )

  expect_equal(
    get_structure_constants(gips_example_perm),
    example_structure_constants
  )
  expect_equal(
    get_structure_constants(gips_example_perm2),
    example_structure_constants2
  )
})

test_that("get_structure_constants handles a permutation of large order", {
  cycle_lengths <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
  p <- sum(cycle_lengths) # 328
  cycles <- unname(split(
    seq_len(p),
    rep(seq_along(cycle_lengths), times = cycle_lengths)
  ))
  perm <- gips_perm(list(cycles), size = p)

  # All 15 cycles share frequency 0, giving r = 15 and dim_omega = 15 * 16 / 2.
  # The distinct odd primes give 156 distinct frequencies strictly between 0 and 1/2.
  # Each has r = 1 and d = 2; only the 2-cycle contributes frequency 1/2, with d = 1.
  expected_constants <- list(
    r = c(15, rep(1, 157)),
    d = c(1, rep(2, 156), 1),
    k = c(1, rep(2, 156), 1),
    L = 158,
    dim_omega = c(120, rep(1, 157))
  )

  expect_silent(constants <- get_structure_constants(perm))
  expect_equal(sum(constants$r * constants$d), p, tolerance = 0)
  expect_equal(constants, expected_constants, tolerance = 0)
})

# Example 6 from the paper
test_that("calculate_r works for example from paper (without 0 values)", {
  expect_equal(
    calculate_r(c(3, 2, 1), 6),
    c(3, 1, 1)
  )
})

test_that("calculate_r works for identity", {
  expect_equal(
    calculate_r(c(1, 1, 1, 1), 1),
    4
  )
})

test_that("calculate_d works for even perm_order", {
  expect_equal(calculate_d(4, 6), c(1, 2, 2, 1))
  expect_equal(calculate_d(3, 4), c(1, 2, 1))
  expect_equal(calculate_d(2, 2), c(1, 1))
})

test_that("calculate_d works for odd perm_order", {
  expect_equal(calculate_d(3, 5), c(1, 2, 2))
  expect_equal(calculate_d(2, 3), c(1, 2))
})

test_that("calculate_d works for identity", {
  expect_equal(calculate_d(1, 1), 1)
})

test_that("get_structure_constants checks for proper argument", {
  expect_silent(get_structure_constants(gips_perm("(1,2,3)", 3)))

  expect_error(
    get_structure_constants(3),
    "You provided `perm` with `class\\(perm\\) == \\(numeric\\)`"
  )

  expect_error(
    get_structure_constants(permutations::as.cycle("(1,2,3)")),
    "You provided `perm` with `class\\(perm\\) == \\(permutation, cycle\\)`"
  )
})

test_that("get_structure_constants rejects permutations with missing fixed points", {
  incomplete <- new_gips_perm(list(1:2), 3)
  expect_error(
    get_structure_constants(incomplete), "complete permutation",
    fixed = TRUE, class = "rlang_error"
  )

  complete <- gips_perm("(1,2)", 3)
  expect_equal(sum(get_structure_constants(complete)$dim_omega), 4)
})

test_that("get_structure_constants can get a gips as perm", {
  perm <- gips_perm(permutations::as.word(c(1, 2, 3, 5, 4)), 5)
  expect_silent(constants1 <- get_structure_constants(perm))
  expect_silent(constants2 <- get_structure_constants(gips(diag(5), 5, perm = perm)))
  expect_equal(constants1, constants2)
})
