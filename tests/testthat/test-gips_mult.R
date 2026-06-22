test_that("gips() accepts a list of matrices (multi-sample construction)", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L))

  expect_s3_class(g, "gips")
  expect_true(is.list(attr(g, "S")))
  expect_equal(length(attr(g, "S")), 2)
  expect_equal(attr(g, "number_of_observations"), c(10L, 12L))
  expect_true(is.list(attr(g, "D_matrix")))
  expect_equal(length(attr(g, "D_matrix")), 2)
})


test_that("gips() multi-sample errors on mismatched matrix dimensions", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- diag(3)

  expect_error(gips(list(S1, S2), c(10L, 12L)))
})


test_that("gips() multi-sample errors on wrong-length number_of_observations", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  expect_error(gips(list(S1, S2), c(10L, 12L, 15L)))
  expect_error(gips(list(S1, S2), 10L))
})


test_that("gips() multi-sample errors when S list contains non-matrices", {
  expect_error(gips(list(matrix(c(1, 0.5, 0.5, 2), nrow = 2), "not a matrix"), c(10L, 12L)))
})


test_that("gips() multi-sample accepts per-group delta vector", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L), delta = c(3, 5))
  expect_equal(attr(g, "delta"), c(3, 5))
})


test_that("gips() multi-sample broadcasts scalar delta to all groups", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L), delta = 4)
  expect_equal(attr(g, "delta"), c(4, 4))
})


test_that("gips() multi-sample errors on wrong-length delta vector", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  expect_error(gips(list(S1, S2), c(10L, 12L), delta = c(3, 4, 5)))
})


test_that("log_posteriori_of_gips() respects per-group delta", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)
  n1 <- 15L; n2 <- 18L
  delta1 <- 3; delta2 <- 5

  g_mult <- gips(list(S1, S2), c(n1, n2), delta = c(delta1, delta2),
                 was_mean_estimated = FALSE)

  # Test there is also the same default for D_matrix:
  g1 <- gips(S1, n1, delta = delta1, was_mean_estimated = FALSE)
  g2 <- gips(S2, n2, delta = delta2, was_mean_estimated = FALSE)

  expect_equal(
    log_posteriori_of_gips(g_mult),
    log_posteriori_of_gips(g1) + log_posteriori_of_gips(g2),
    tolerance = 1e-10
  )
})


test_that("gips() multi-sample accepts a custom D_matrix list", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L), D_matrix = list(diag(2), 2*diag(2)))
  expect_s3_class(g, "gips")
})


test_that("log_posteriori_of_gips() on multi-sample equals sum of single-sample log-posteriors", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)
  n1 <- 15L
  n2 <- 18L

  g_mult <- gips(list(S1, S2), c(n1, n2), was_mean_estimated = FALSE)

  D1 <- attr(g_mult, "D_matrix")[[1]]
  D2 <- attr(g_mult, "D_matrix")[[2]]
  deltas <- attr(g_mult, "delta") # c(3, 3) after broadcasting

  g1 <- gips(S1, n1, delta = deltas[1], D_matrix = D1, was_mean_estimated = FALSE)
  g2 <- gips(S2, n2, delta = deltas[2], D_matrix = D2, was_mean_estimated = FALSE)

  expect_equal(
    log_posteriori_of_gips(g_mult),
    log_posteriori_of_gips(g1) + log_posteriori_of_gips(g2),
    tolerance = 1e-10
  )
})


test_that("find_MAP() with BF optimizer works on multi-sample gips (p=3, G=2)", {
  p <- 3
  Sigma <- diag(p)
  S1 <- rWishart(1, df = 10, Sigma = Sigma)[, , 1] / 10
  S2 <- rWishart(1, df = 12, Sigma = Sigma)[, , 1] / 12

  g <- gips(list(S1, S2), c(10L, 12L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)

  expect_s3_class(g_map, "gips")
  expect_true(is.list(attr(g_map, "S")))
  expect_false(is.null(attr(g_map, "optimization_info")))
})


test_that("find_MAP() with MH optimizer and continue works on larger multi-sample gips (p=30, G=2)", {
  p <- 30
  Sigma <- diag(p)
  S1 <- rWishart(1, df = 100, Sigma = Sigma)[, , 1] / 100
  S2 <- rWishart(1, df = 120, Sigma = Sigma)[, , 1] / 120
  
  g <- gips(list(S1, S2), c(10L, 12L), perm = "(1,2,3,4,5,6,7,8,9,10,11,12,13,14)")
  
  g_map <- find_MAP(g, optimizer = "MH", max_iter = 20, show_progress_bar = FALSE)
  g_map_2 <- find_MAP(g_map, optimizer = "continue", max_iter = 20, show_progress_bar = FALSE)
  
  expect_s3_class(g_map, "gips")
  expect_s3_class(g_map_2, "gips")
  expect_true(is.list(attr(g_map, "S")))
  expect_true(is.list(attr(g_map_2, "S")))
  expect_false(is.null(attr(g_map, "optimization_info")))
  expect_false(is.null(attr(g_map_2, "optimization_info")))
  # The continued optimization should have more log_posteriori_values
  expect_gt(
    length(attr(g_map_2, "optimization_info")[["log_posteriori_values"]]),
    length(attr(g_map, "optimization_info")[["log_posteriori_values"]])
  )
})


test_that("project_matrix() returns a list when S is a list", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L))
  projected <- project_matrix(list(S1, S2), g[[1]])

  expect_true(is.list(projected))
  expect_equal(length(projected), 2L)
  expect_true(is.matrix(projected[[1]]))
  expect_true(is.matrix(projected[[2]]))
  # Each element must equal the single-group projection
  expect_equal(projected[[1]], project_matrix(S1, g[[1]]))
  expect_equal(projected[[2]], project_matrix(S2, g[[1]]))
})


test_that("logLik() works on multi-sample gips", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(15L, 18L))
  ll <- logLik(g)

  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  # df = G * single-group n_parameters
  expect_equal(attr(ll, "df"), 2L * sum(get_structure_constants(g[[1]])[["dim_omega"]]))
  # nobs = sum of all observations
  expect_equal(attr(ll, "nobs"), 15L + 18L)
})


test_that("logLik() returns NULL when min(n_g) < n0", {
  # For perm = "()" with p=3, n0 = 3; use n=2 per group
  S1 <- diag(3)
  S2 <- diag(3)

  g <- gips(list(S1, S2), c(2L, 2L), was_mean_estimated = FALSE)
  expect_warning(ll <- logLik(g))
  expect_null(ll)
})


test_that("AIC() and BIC() work on multi-sample gips", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(15L, 18L))

  expect_true(is.finite(AIC(g)))
  expect_true(is.finite(BIC(g)))
})


test_that("print() does not error on multi-sample gips", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L))
  expect_no_error(capture.output(print(g)))
})


test_that("summary() does not error on multi-sample gips and returns correct fields", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L))
  s <- summary(g)

  expect_s3_class(s, "summary.gips")
  expect_true(isTRUE(s[["is_multi_sample"]]))
  expect_equal(s[["G"]], 2L)
  expect_equal(s[["number_of_observations"]], c(10L, 12L))
  expect_no_error(capture.output(print(s)))
})


test_that("summary() displays correct degrees of freedom for multi-sample with different group sizes", {
  # Regression test for bug where degrees of freedom only considered first group
  S1 <- diag(3)
  S2 <- diag(3)
  n1 <- 15L
  n2 <- 18L

  # Test with mean estimated (default)
  g_est <- gips(list(S1, S2), c(n1, n2), was_mean_estimated = TRUE)
  s_est <- summary(g_est)
  output_est <- capture.output(print(s_est))
  
  # Should show degrees of freedom for each group (n - 1)
  expect_match(paste(output_est, collapse = "\n"), "Degrees of freedom left \\(per group\\): 14, 17")

  # Test with mean not estimated
  g_not_est <- gips(list(S1, S2), c(n1, n2), was_mean_estimated = FALSE)
  s_not_est <- summary(g_not_est)
  output_not_est <- capture.output(print(s_not_est))
  
  # Should show actual sample sizes (15, 18)
  expect_match(paste(output_not_est, collapse = "\n"), "all degrees of freedom were preserved \\(15, 18\\)")
})


test_that("plot() errors with a not-yet-implemented message on multi-sample gips", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L))
  expect_error(plot(g, type = "heatmap"), regexp = "not yet implemented")
})


test_that("compare_log_posteriories_of_perms() works on multi-sample gips", {
  p <- 5
  Sigma <- diag(p)
  S1 <- rWishart(1, df = 100, Sigma = Sigma)[, , 1] / 100
  S2 <- rWishart(1, df = 120, Sigma = Sigma)[, , 1] / 120

  g_1 <- gips(list(S1, S2), c(10L, 12L), perm = "(1,2,3)")
  g_2 <- gips(list(S1, S2), c(10L, 12L), perm = "(1,2,5)")

  # The function should work without error
  expect_output(expect_no_error(compare_log_posteriories_of_perms(g_1, g_2)))

  # Manually check the ratio matches the log-posterior difference
  log_post_1 <- log_posteriori_of_gips(g_1)
  log_post_2 <- log_posteriori_of_gips(g_2)
  expected_ratio <- log_post_1 - log_post_2

  expect_output(result <- compare_log_posteriories_of_perms(g_1, g_2))
  expect_equal(result, expected_ratio, tolerance = 1e-10)
})


test_that("single-sample gips() is unaffected by multi-sample changes", {
  S <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  g <- gips(S, 10L)

  expect_s3_class(g, "gips")
  expect_true(is.matrix(attr(g, "S")))
  expect_equal(attr(g, "number_of_observations"), 10L)

  ll <- logLik(g)
  expect_s3_class(ll, "logLik")
})
