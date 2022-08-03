test_that("orthogonal matrix fulfills condition from example 6", {
  example_basis <- diag(6)
  gips_example_perm <- gips_perm(example_perm, 6)

  U_Gamma <- prepare_orthogonal_matrix(gips_example_perm, example_basis)
  tested_matrix <- t(U_Gamma) %*% matrix_invariant_by_example_perm %*% U_Gamma
  tested_elements <- tested_matrix
  tested_elements[1:3, 1:3] <- 0
  tested_elements[4, 4] <- 0
  tested_elements[5, 5] <- 0
  tested_elements[6, 6] <- 0
  expect_equal(
    tested_elements,
    matrix(rep(0, 6 * 6), ncol = 6)
  )
  expect_true(isSymmetric(tested_matrix[1:3, 1:3]))
  expect_equal(tested_matrix[4, 4], tested_matrix[5, 5])
})

test_that("prepare_orthogonal_matrix works for example", {
  example_orth_matrix <- matrix(c(
    rep(1 / sqrt(3), 3), rep(0, 3),
    rep(0, 3), rep(1 / sqrt(2), 2), 0,
    rep(0, 5), 1,
    sqrt(2 / 3), rep(-sqrt(1 / 6), 2), rep(0, 3),
    0, 1 / sqrt(2), -1 / sqrt(2), rep(0, 3),
    rep(0, 3), 1 / sqrt(2), -1 / sqrt(2), 0
  ), nrow = 6)
  example_basis <- diag(6)

  gips_example_perm <- gips_perm(example_perm, 6)
  expect_equal(
    prepare_orthogonal_matrix(gips_example_perm, example_basis),
    example_orth_matrix
  )
})

test_that("get_v_matrix_for_subcycle works for 3-length cycle", {
  example_v_object <- list(
    matrix(c(
      rep(1 / sqrt(3), 3), rep(0, 3),
      sqrt(2 / 3), rep(-sqrt(1 / 6), 2), rep(0, 3),
      0, 1 / sqrt(2), -1 / sqrt(2), rep(0, 3)
    ), nrow = 6),
    matrix(c(
      rep(0, 3), 1 / sqrt(2), 1 / sqrt(2), 0,
      rep(0, 3), 1 / sqrt(2), -1 / sqrt(2), 0
    ), nrow = 6),
    matrix(c(rep(0, 5), 1), nrow = 6)
  )
  example_basis <- diag(6)

  expect_equal(
    get_v_matrix_for_subcycle(
      1:3,
      example_basis
    ),
    example_v_object[[1]]
  )
})

test_that("get_v_matrix_for_subcycle works for 2-length cycle", {
  example_v_object <- list(
    matrix(c(
      rep(1 / sqrt(3), 3), rep(0, 3),
      sqrt(2 / 3), rep(-sqrt(1 / 6), 2), rep(0, 3),
      0, 1 / sqrt(2), -1 / sqrt(2), rep(0, 3)
    ), nrow = 6),
    matrix(c(
      rep(0, 3), 1 / sqrt(2), 1 / sqrt(2), 0,
      rep(0, 3), 1 / sqrt(2), -1 / sqrt(2), 0
    ), nrow = 6),
    matrix(c(rep(0, 5), 1), nrow = 6)
  )
  example_basis <- diag(6)

  expect_equal(
    get_v_matrix_for_subcycle(
      4:5,
      example_basis
    ),
    example_v_object[[2]]
  )
})

test_that("get_v_matrix_for_subcycle works for identity", {
  example_v_object <- list(
    matrix(c(
      rep(1 / sqrt(3), 3), rep(0, 3),
      sqrt(2 / 3), rep(-sqrt(1 / 6), 2), rep(0, 3),
      0, 1 / sqrt(2), -1 / sqrt(2), rep(0, 3)
    ), nrow = 6),
    matrix(c(
      rep(0, 3), 1 / sqrt(2), 1 / sqrt(2), 0,
      rep(0, 3), 1 / sqrt(2), -1 / sqrt(2), 0
    ), nrow = 6),
    matrix(c(rep(0, 5), 1), nrow = 6)
  )
  example_basis <- diag(6)

  expect_equal(
    get_v_matrix_for_subcycle(
      6,
      example_basis
    ),
    example_v_object[[3]]
  )
})

test_that("arrange_v_object works for example", {
  example_v_object <- list(
    matrix(c(
      rep(1 / sqrt(3), 3), rep(0, 3),
      sqrt(2 / 3), rep(-sqrt(1 / 6), 2), rep(0, 3),
      0, 1 / sqrt(2), -1 / sqrt(2), rep(0, 3)
    ), nrow = 6),
    matrix(c(
      rep(0, 3), 1 / sqrt(2), 1 / sqrt(2), 0,
      rep(0, 3), 1 / sqrt(2), -1 / sqrt(2), 0
    ), nrow = 6),
    matrix(c(rep(0, 5), 1), nrow = 6)
  )
  example_orth_matrix <- matrix(c(
    rep(1 / sqrt(3), 3), rep(0, 3),
    rep(0, 3), rep(1 / sqrt(2), 2), 0,
    rep(0, 5), 1,
    sqrt(2 / 3), rep(-sqrt(1 / 6), 2), rep(0, 3),
    0, 1 / sqrt(2), -1 / sqrt(2), rep(0, 3),
    rep(0, 3), 1 / sqrt(2), -1 / sqrt(2), 0
  ), nrow = 6)

  expect_equal(
    arrange_v_object(example_v_object),
    example_orth_matrix
  )
})
