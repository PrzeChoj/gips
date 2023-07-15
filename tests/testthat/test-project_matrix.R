test_that("projected matrix is invariant by example_perm", {
  S <- matrix(rnorm(6 * 6), ncol = 6)
  S <- S %*% t(S)
  gips_example_perm <- gips_perm(example_perm, 6)

  projected <- project_matrix(S, gips_example_perm)
  expect_true(isSymmetric(projected))
  # now we need to check only lower diagonal
  expect_equal(
    rep(projected[1, 1], 2),
    c(projected[2, 2], projected[3, 3])
  )
  expect_equal(
    rep(projected[2, 1], 2),
    projected[3, 1:2]
  )
  expect_equal(
    rep(projected[4, 1], 6),
    as.vector(projected[4:5, 1:3])
  )
  expect_equal(
    rep(projected[6, 1], 2),
    projected[6, 2:3]
  )
  expect_equal(projected[4, 4], projected[5, 5])
  expect_equal(projected[6, 4], projected[6, 5])

  # other ways of computing the same outcome
  projected2 <- project_matrix(S, "(1,2,3)(4,5)")
  expect_identical(projected, projected2)

  precomputed_equal_indices <- get_equal_indices_by_perm(gips_perm(example_perm, 6))
  projected3 <- project_matrix(S, example_perm, precomputed_equal_indices)
  expect_identical(projected, projected3)
})

test_that("project_matrix gives errors", {
  gips_example_perm <- gips_perm(example_perm, 6)

  expect_error(project_matrix(7, gips_example_perm))
  expect_error(project_matrix(matrix(rnorm(5 * 6), ncol = 6), gips_example_perm))
  S <- matrix(rnorm(5 * 5), ncol = 5)
  S <- S %*% t(S)
  expect_error(project_matrix(S, gips_example_perm))

  expect_warning(project_matrix(matrix(rnorm(36), nrow = 6), gips_example_perm),
    class = "not_positive_semi_definite_matrix"
  )
})

test_that("project_matrix can get gips as per", {
  p <- 6
  my_perm <- "(14)(23)"
  number_of_observations <- 10
  X <- matrix(rnorm(p * number_of_observations), number_of_observations, p)
  S <- cov(X)
  projected_S <- project_matrix(S, perm = my_perm)

  g <- gips(S, number_of_observations, perm = my_perm)
  g_MAP <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")

  S_MAP1 <- project_matrix(attr(g, "S"), perm = g_MAP[[1]])               # gips_perm class
  S_MAP2 <- project_matrix(attr(g, "S"), perm = g_MAP)                    # gips class
  S_MAP3 <- project_matrix(attr(g, "S"), perm = as.character(g_MAP[[1]])) # character

  expect_equal(S_MAP1, S_MAP2)
  expect_equal(S_MAP1, S_MAP3)
})

test_that("project_matrix does not forget colnames or rownames", {
  p <- 9
  S <- matrix(rnorm(p * p), nrow = p)
  S <- S %*% t(S)
  rownames(S) <- LETTERS[1:p]
  colnames(S) <- LETTERS[1:p]

  S_proj <- project_matrix(S, "(123)")

  expect_equal(rownames(S_proj), rownames(S))
  expect_equal(colnames(S_proj), colnames(S))
})

test_that("get_equal_indices_by_perm works for example_perm", {
  values <- unique(as.integer(matrix_invariant_by_example_perm))
  expected_equal_indices_by_example_perm <- lapply(values, function(v) {
    which(as.integer(matrix_invariant_by_example_perm) == v)
  })
  gips_example_perm <- gips_perm(example_perm, 6)

  actual_l <- lapply(
    get_equal_indices_by_perm(gips_example_perm),
    sort
  )
  expected_l <- lapply(
    expected_equal_indices_by_example_perm,
    sort
  )
  expect_setequal(actual_l, expected_l)
})

test_that("get_equal_indices_by_perm works for identity", {
  expect_setequal(
    get_equal_indices_by_perm(gips_perm(to_perm(1:3), 3)),
    list(1, c(2, 4), c(3, 7), 5, c(6, 8), 9)
  )
})

test_that("get_single_from_double_indices works", {
  full_double_indices <- matrix(
    c(
      rep(1:4, times = 4),
      rep(1:4, each = 4)
    ),
    ncol = 2
  )

  expect_equal(
    get_single_from_double_indices(full_double_indices, 4),
    1:16
  )
})

test_that("get_single_from_double_indices works for 0 input", {
  expect_equal(
    get_single_from_double_indices(matrix(numeric(0), ncol = 2), 4),
    numeric(0)
  )
})

test_that("get_double_from_single_indices works", {
  full_double_indices <- matrix(
    c(
      rep(1:4, times = 4),
      rep(1:4, each = 4)
    ),
    ncol = 2
  )

  expect_equal(
    get_double_from_single_indices(1:16, 4),
    full_double_indices
  )
  expect_equal(
    get_double_from_single_indices(as.vector(matrix(1:16, ncol = 4)), 4),
    full_double_indices
  )
})

test_that("get_double_from_single_indices works for 0 input", {
  expect_equal(
    get_double_from_single_indices(numeric(0), 4),
    matrix(numeric(0), ncol = 2)
  )
})
