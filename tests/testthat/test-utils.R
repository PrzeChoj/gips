test_that("shift_vector works", {
  expect_equal(
    shift_vector(c(1, 2, 3), 1),
    c(2, 3, 1)
  )

  v <- sample(60)
  expect_equal(
    shift_vector(v, 0),
    v
  )
  expect_equal(
    shift_vector(v, 60),
    v
  )
})

test_that("rearrange_vector works", {
  expect_equal(
    rearrange_vector(c(4, 2, 3, 1)),
    c(1, 4, 2, 3)
  )
})

test_that("is.wholenumber works", {
  expect_true(is.wholenumber(3L))
  expect_true(is.wholenumber(as.double(3L)))
  expect_false(is.wholenumber(3L / 2))
})

test_that("is.positive.semi.definite.matrix works", {
  p <- 7
  random_matrix <- matrix(rnorm(p * p), nrow = p)
  random_matrix <- t(random_matrix) %*% random_matrix

  random_matrix <- random_matrix - diag(p) * (eigen(random_matrix, symmetric = TRUE, only.values = TRUE)[["values"]][p] - 10)
  # the smallest eigen value of `random_matrix` is very close to 10

  expect_true(is.positive.semi.definite.matrix(random_matrix))
  expect_true(is.positive.semi.definite.matrix(random_matrix - diag(p) * 10, tolerance = 0.001)) # the smallest eigen value is very close to 0
  expect_false(is.positive.semi.definite.matrix(random_matrix - diag(p) * (10.01))) # the smallest eigen value is very close to -0.01
})

test_that("change_log_probabilities_unnorm_to_probabilities forks as intendent", {
  desired_probs <- c(0.1, 0.2, 0.3, 0.4)
  changed_probs <- log(desired_probs) + runif(1, -1000, 1000)
  est_desired_probs <- change_log_probabilities_unnorm_to_probabilities(changed_probs)

  expect_equal(desired_probs, est_desired_probs) # function can recreate the original distribution
})

test_that("pretty_plot work", {
  gg_plot_object <- pretty_plot_matrix(matrix_invariant_by_example_perm)
  expect_true(inherits(gg_plot_object, "ggplot"))

  gg_plot_object2 <- pretty_plot_block_matrix(
    matrix_invariant_by_example_perm,
    example_perm
  )
  expect_true(inherits(gg_plot_object2, "ggplot"))
})
