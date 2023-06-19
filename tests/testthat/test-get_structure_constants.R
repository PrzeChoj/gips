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

# Example 6 from the paper
test_that("calculate_r works for example from paper", {
  expect_equal(
    calculate_r(c(3, 2, 1), 6),
    c(3, 0, 1, 1)
  )
})

test_that("calculate_r works for identity", {
  expect_equal(
    calculate_r(c(1, 1, 1, 1), 1),
    4
  )
})

test_that("calculate_d works for even perm_order", {
  expect_equal(calculate_d(6), c(1, 2, 2, 1))
  expect_equal(calculate_d(4), c(1, 2, 1))
  expect_equal(calculate_d(2), c(1, 1))
})

test_that("calculate_d works for odd perm_order", {
  expect_equal(calculate_d(5), c(1, 2, 2))
  expect_equal(calculate_d(3), c(1, 2))
})

test_that("calculate_d works for identity", {
  expect_equal(calculate_d(1), 1)
})

test_that("get_structure_constants checks for proper argument", {
  expect_silent(get_structure_constants(gips_perm("(1,2,3)", 3)))
  
  expect_error(get_structure_constants(3),
               "You provided `perm` with `class\\(perm\\) == \\(numeric\\)`")
  
  expect_error(get_structure_constants(permutations::as.cycle("(1,2,3)")),
               "You provided `perm` with `class\\(perm\\) == \\(permutation, cycle\\)`")
})

test_that("get_structure_constants can get a gips as perm", {
  perm <- gips_perm(permutations::as.word(c(1, 2, 3, 5, 4)), 5)
  expect_silent(constants1 <- get_structure_constants(perm))
  expect_silent(constants2 <- get_structure_constants(gips(diag(5), 5, perm = perm)))
  expect_equal(constants1, constants2)
})
