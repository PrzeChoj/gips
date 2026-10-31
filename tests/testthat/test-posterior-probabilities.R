# These tests specify the 1.3.1 contract: posterior weights normalized over
# distinct visited cyclic groups, rather than estimates from visit counts.
posterior_test_bf <- function(g) {
  fit <- find_MAP(g,
    optimizer = "BF", save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )
  get_probabilities_from_gips(fit)
}

posterior_test_visited_groups <- function(g) {
  perms <- attr(g, "optimization_info")$visited_perms
  unique(vapply(perms, function(perm) {
    as.character(get_group_representative(perm))
  }, character(1)))
}

test_that("MH probabilities equal BF after visiting every cyclic group", {
  g <- gips(diag(c(1, 4, 12)), 100, was_mean_estimated = FALSE)
  reference <- posterior_test_bf(g)
  fit <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 50, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))

  # This short run visits all five groups; no sampling-convergence tolerance
  # is appropriate when probabilities are calculated from their log weights.
  groups <- posterior_test_visited_groups(fit)
  expect_length(reference, 5)
  expect_setequal(groups, names(reference))
  expect_gt(unname(reference["()"]), 0.999)

  actual <- get_probabilities_from_gips(fit)
  expect_type(actual, "double")
  expect_length(actual, length(groups))
  expect_setequal(names(actual), groups)
  expect_true(all(is.finite(actual)))
  expect_true(all(actual >= 0))
  expect_equal(sum(actual), 1, tolerance = 1e-12)
  expect_true(all(diff(unname(actual)) <= 0))
  expect_equal(actual[names(reference)], reference, tolerance = 1e-12)
})

test_that("calculating probabilities preserves the MH path and its plateaus", {
  g <- gips(diag(c(1, 4, 12)), 100, was_mean_estimated = FALSE)
  run_mh <- function(return_probabilities) {
    withr::with_seed(1, suppressMessages(find_MAP(g,
      optimizer = "MH", max_iter = 50, save_all_perms = TRUE,
      return_probabilities = return_probabilities, show_progress_bar = FALSE
    )))
  }
  without_probabilities <- run_mh(FALSE)
  with_probabilities <- run_mh(TRUE)
  without_info <- attr(without_probabilities, "optimization_info")
  with_info <- attr(with_probabilities, "optimization_info")
  without_states <- vapply(without_info$visited_perms, as.character, character(1))
  with_states <- vapply(with_info$visited_perms, as.character, character(1))

  # The fixture includes both moves and repeated states, making preservation
  # of the numeric plot data meaningful. Execution times vary between runs.
  expect_gt(length(unique(without_states)), 1L)
  expect_true(any(rle(without_states)$lengths > 1L))
  expect_identical(with_states, without_states)
  expect_equal(with_info$log_posteriori_values, without_info$log_posteriori_values,
    tolerance = 1e-12
  )
  expect_identical(with_info$all_n0, without_info$all_n0)
  expect_identical(as.character(with_info$last_perm), as.character(without_info$last_perm))
  expect_identical(as.character(with_probabilities), as.character(without_probabilities))
})

test_that("MH probabilities condition BF on the visited subset of groups", {
  g <- gips(diag(3), 10, was_mean_estimated = FALSE)
  reference <- posterior_test_bf(g)
  fit <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 3, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))

  groups <- posterior_test_visited_groups(fit)
  expect_gt(length(groups), 1)
  expect_lt(length(groups), length(reference))
  expect_lt(sum(reference[groups]), 1 - 1e-6)
  expected <- reference[groups] / sum(reference[groups])

  actual <- get_probabilities_from_gips(fit)
  expect_length(actual, length(groups))
  expect_setequal(names(actual), groups)
  expect_equal(actual[groups], expected, tolerance = 1e-12)
  expect_equal(sum(actual), 1, tolerance = 1e-12)
})

test_that("multi-sample MH probabilities use the joint posterior weights", {
  matrices <- list(
    diag(c(1, 2, 4)),
    matrix(c(2, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 3), nrow = 3)
  )
  n <- c(12, 17)
  delta <- c(3, 4)
  D <- list(diag(3), diag(c(2, 3, 4)))
  g <- gips(matrices, n,
    delta = delta, D_matrix = D, was_mean_estimated = TRUE
  )

  # Independent samples share a group. With the uniform group prior, their
  # individual BF probabilities multiply to give the joint unnormalized weights.
  individual <- lapply(seq_along(matrices), function(i) {
    posterior_test_bf(gips(matrices[[i]], n[i],
      delta = delta[i], D_matrix = D[[i]], was_mean_estimated = TRUE
    ))
  })
  group_names <- names(individual[[1]])
  joint_weights <- individual[[1]] * individual[[2]][group_names]

  fit <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 50, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))
  groups <- posterior_test_visited_groups(fit)
  expect_gt(length(groups), 1)
  expected <- joint_weights[groups] / sum(joint_weights[groups])

  actual <- get_probabilities_from_gips(fit)
  expect_length(actual, length(groups))
  expect_setequal(names(actual), groups)
  expect_equal(actual[groups], expected, tolerance = 1e-12)
})

test_that("the getter does not calculate unsaved MH probabilities in 1.3.1", {
  g <- gips(diag(c(1, 1.4, 2)), 8, was_mean_estimated = FALSE)
  fit <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 3, save_all_perms = TRUE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  )))

  # TODO (1.4.0): calculate probabilities on demand from saved models and
  # log weights, without running the optimizer again; update this expectation.
  actual <- suppressMessages(get_probabilities_from_gips(fit))
  expect_null(actual)
})

test_that("forget_perms preserves probabilities already calculated after MH", {
  g <- gips(diag(c(1, 1.4, 2)), 8, was_mean_estimated = FALSE)
  fit <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 10, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))
  before <- get_probabilities_from_gips(fit)
  expect_type(before, "double")
  expect_gt(length(before), 1)

  forgotten <- forget_perms(fit)
  expect_equal(get_probabilities_from_gips(forgotten), before)
})
