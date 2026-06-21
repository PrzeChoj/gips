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


test_that("gips() multi-sample accepts a custom D_matrix list", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L), D_matrix = list(diag(2), diag(2)))
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
  delta <- attr(g_mult, "delta")

  g1 <- gips(S1, n1, delta = delta, D_matrix = D1, was_mean_estimated = FALSE)
  g2 <- gips(S2, n2, delta = delta, D_matrix = D2, was_mean_estimated = FALSE)

  expect_equal(
    log_posteriori_of_gips(g_mult),
    log_posteriori_of_gips(g1) + log_posteriori_of_gips(g2),
    tolerance = 1e-10
  )
})


test_that("find_MAP() with BF optimizer works on multi-sample gips (p=3, G=2)", {
  set.seed(7531)
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


test_that("plot() errors with a not-yet-implemented message on multi-sample gips", {
  S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
  S2 <- matrix(c(2, 0.3, 0.3, 1.5), nrow = 2, byrow = TRUE)

  g <- gips(list(S1, S2), c(10L, 12L))
  expect_error(plot(g, type = "heatmap"), regexp = "not yet implemented")
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
