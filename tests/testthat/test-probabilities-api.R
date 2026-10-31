probabilities_api_fit <- function(optimizer = "MH", return_probabilities = TRUE) {
  g <- gips(diag(c(1, 1.4, 2)), 8, was_mean_estimated = FALSE)
  withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = optimizer, max_iter = 10, save_all_perms = TRUE,
    return_probabilities = return_probabilities, show_progress_bar = FALSE
  )))
}

probabilities_api_as_legacy <- function(g) {
  info <- attr(g, "optimization_info")
  info$visited_log_posteriori_values <- NULL
  info$post_probabilities_method <- NULL
  info$post_probabilities_scope <- NULL
  attr(g, "optimization_info") <- info
  g
}

test_that("new MH and BF probabilities identify their method and scope", {
  for (optimizer in c("MH", "BF")) {
    fit <- probabilities_api_fit(optimizer)
    info <- attr(fit, "optimization_info")

    expect_identical(info$post_probabilities_method, "posterior_weights")
    expect_identical(info$post_probabilities_scope,
      if (optimizer == "BF") "all_groups" else "visited_groups"
    )
    if (optimizer == "MH") {
      expect_type(info$visited_log_posteriori_values, "double")
      expect_length(info$visited_log_posteriori_values, length(info$visited_perms))
    }
    expect_silent(actual <- get_probabilities_from_gips(fit))
    expect_identical(actual, info$post_probabilities)
  }
})

test_that("legacy MH probabilities warn and are returned unchanged", {
  fit <- probabilities_api_as_legacy(probabilities_api_fit())
  old_probabilities <- c("()" = 0.8, "(1,2)" = 0.2)

  for (n_runs in c(1L, 2L)) {
    info <- attr(fit, "optimization_info")
    info$optimization_algorithm_used <- rep("Metropolis_Hastings", n_runs)
    info$post_probabilities <- old_probabilities
    attr(fit, "optimization_info") <- info

    expect_warning(actual <- get_probabilities_from_gips(fit),
      regexp = "visit-frequency estimates", class = "gips_legacy_probabilities"
    )
    expect_identical(actual, old_probabilities)
    expect_identical(attr(fit, "optimization_info"), info)
  }
})

test_that("legacy BF probabilities remain silently readable", {
  fit <- probabilities_api_as_legacy(probabilities_api_fit("BF"))
  before <- attr(fit, "optimization_info")

  expect_silent(actual <- get_probabilities_from_gips(fit))
  expect_identical(actual, before$post_probabilities)
  expect_identical(attr(fit, "optimization_info"), before)
})

test_that("missing MH probabilities are neither computed nor marked as legacy", {
  fit <- probabilities_api_fit(return_probabilities = FALSE)

  for (candidate in list(fit, probabilities_api_as_legacy(fit))) {
    before <- attr(candidate, "optimization_info")
    expect_gt(length(before$visited_perms), 1L)
    expect_warning(
      actual <- suppressMessages(get_probabilities_from_gips(candidate)), NA
    )
    expect_null(actual)
    expect_identical(attr(candidate, "optimization_info"), before)
  }
})

test_that("MH rejects probabilities without saved states before optimization", {
  inputs <- list(
    gips(diag(3), 8, was_mean_estimated = FALSE),
    gips(list(diag(3), 2 * diag(3)), c(8, 10), was_mean_estimated = FALSE)
  )

  for (g in inputs) {
    withr::with_seed(1, {
      seed_before <- .Random.seed
      expect_error(suppressMessages(find_MAP(g,
        optimizer = "MH", max_iter = 3, save_all_perms = FALSE,
        return_probabilities = TRUE, show_progress_bar = FALSE
      )), "save_all_perms == FALSE")
      expect_identical(.Random.seed, seed_before)
    })
  }
})

test_that("forget_perms drops saved states and weights but keeps probabilities and plot data", {
  for (optimizer in c("MH", "BF")) {
    fit <- probabilities_api_fit(optimizer)
    before <- attr(fit, "optimization_info")
    if (optimizer == "MH") {
      expect_gt(length(before$visited_log_posteriori_values), 1L)
    }

    expect_silent(forgotten <- forget_perms(fit))
    after <- attr(forgotten, "optimization_info")
    expect_true(all(is.na(after$visited_perms)))
    expect_null(after$visited_log_posteriori_values)
    expect_identical(after$post_probabilities, before$post_probabilities)
    expect_identical(after$post_probabilities_method, before$post_probabilities_method)
    expect_identical(after$post_probabilities_scope, before$post_probabilities_scope)
    expect_identical(after$log_posteriori_values, before$log_posteriori_values)
    expect_identical(after$all_n0, before$all_n0)
    expect_silent(actual <- get_probabilities_from_gips(forgotten))
    expect_identical(actual, before$post_probabilities)
  }
})
