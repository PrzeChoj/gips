test_that("projected matrix is invariant by example_perm", {
  S <- matrix(rnorm(6 * 6), ncol = 6)
  S <- tcrossprod(S)
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
})

test_that("project_matrix gives errors", {
  gips_example_perm <- gips_perm(example_perm, 6)

  expect_error(project_matrix(7, gips_example_perm))
  expect_error(project_matrix(matrix(rnorm(5 * 6), ncol = 6), gips_example_perm))
  S <- matrix(rnorm(5 * 5), ncol = 5)
  S <- tcrossprod(S)
  expect_error(project_matrix(S, gips_example_perm))

  expect_warning(project_matrix(matrix(rnorm(36), nrow = 6), gips_example_perm),
    class = "not_positive_semi_definite_matrix"
  )
})

test_that("project_matrix validates permutations before checking their size", {
  incomplete <- new_gips_perm(list(1:2), 3)
  duplicate <- new_gips_perm(list(1:2, 2:3), 3)
  invalid_size <- gips_perm("(1,2)", 3)
  attr(invalid_size, "size") <- NA_real_

  for (S in list(diag(3), list(diag(3), 2 * diag(3)))) {
    expect_error(
      project_matrix(S, incomplete), "complete permutation",
      fixed = TRUE, class = "rlang_error"
    )
    expect_error(
      project_matrix(S, duplicate), "must not repeat",
      fixed = TRUE, class = "rlang_error"
    )
    expect_error(
      project_matrix(S, invalid_size), "single finite non-negative integer",
      fixed = TRUE, class = "rlang_error"
    )
    expect_error(
      project_matrix(S, gips_perm("()", 4)), "Size of `perm` must be equal",
      fixed = TRUE, class = "rlang_error"
    )
  }
})

test_that("project_matrix preserves valid single- and multi-sample projections", {
  S <- diag(c(1, 3, 5))
  expected <- diag(c(2, 2, 5))
  for (perm in list("(1,2)", gips_perm("(1,2)", 3))) {
    expect_equal(project_matrix(S, perm), expected)
    expect_equal(
      project_matrix(list(S, 2 * S), perm),
      list(expected, 2 * expected)
    )
  }
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
  S <- tcrossprod(S)
  rownames(S) <- LETTERS[1:p]
  colnames(S) <- LETTERS[1:p]

  S_proj <- project_matrix(S, "(123)")

  expect_equal(rownames(S_proj), rownames(S))
  expect_equal(colnames(S_proj), colnames(S))
})

test_that("project_matrices_cpp_ projects a list of matrices", {
  p <- 6
  S1 <- matrix(rnorm(p * p), nrow = p)
  S1 <- tcrossprod(S1)
  S2 <- matrix(rnorm(p * p), nrow = p)
  S2 <- tcrossprod(S2)
  rownames(S1) <- rownames(S2) <- LETTERS[1:p]
  colnames(S1) <- colnames(S2) <- letters[1:p]
  perm <- gips_perm("(1,2,3)(4,5)", p)

  projected_list <- project_matrices_cpp_(list(first = S1, second = S2), perm)

  expect_type(projected_list, "list")
  expect_named(projected_list, c("first", "second"))
  expect_equal(projected_list[[1]], project_matrix(S1, perm))
  expect_equal(projected_list[[2]], project_matrix(S2, perm))
  expect_equal(rownames(projected_list[[1]]), rownames(S1))
  expect_equal(colnames(projected_list[[2]]), colnames(S2))
})

test_that("project_matrix handles identity, transpositions, long cycles, and fixed points", {
  S <- matrix(seq_len(36), nrow = 6)
  S <- S + t(S) + diag(6) * 100

  identity_projected <- project_matrix(S, "()")
  expect_equal(identity_projected, S)

  transposition_projected <- project_matrix(S, "(1,2)")
  expect_equal(transposition_projected[1, 1], transposition_projected[2, 2])
  expect_equal(transposition_projected[1, 3], transposition_projected[2, 3])
  expect_equal(transposition_projected[3, 1], transposition_projected[3, 2])

  long_cycle_projected <- project_matrix(S, "(1,2,3,4,5,6)")
  for (i in 1:6) {
    for (j in 1:6) {
      expect_equal(
        long_cycle_projected[i, j],
        long_cycle_projected[ifelse(i == 6, 1, i + 1), ifelse(j == 6, 1, j + 1)]
      )
    }
  }

  fixed_point_projected <- project_matrix(S, "(1,2,3)")
  expect_equal(fixed_point_projected[4, 4], S[4, 4])
  expect_equal(fixed_point_projected[5, 5], S[5, 5])
  expect_equal(fixed_point_projected[6, 6], S[6, 6])
})
