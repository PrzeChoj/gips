legacy_probability_object <- function(g) {
  info <- attr(g, "optimization_info")
  for (field in c("visited_log_posteriori_values",
    "post_probabilities_method", "post_probabilities_scope")) {
    info[[field]] <- NULL
  }
  attr(g, "optimization_info") <- info
  g
}

# Enumerate small cyclic groups independently of the probability implementation.
legacy_probability_group_key <- function(perm, p) {
  if (!inherits(perm, "gips_perm")) {
    perm <- gips_perm(perm, p)
  }
  word <- seq_len(p)
  for (cycle in unclass(perm)) {
    if (length(cycle) > 1L) {
      word[cycle] <- c(cycle[-1L], cycle[1L])
    }
  }
  power <- seq_len(p)
  powers <- character()
  repeat {
    powers <- c(powers, paste(power, collapse = ","))
    power <- word[power]
    if (all(power == seq_len(p))) {
      break
    }
  }
  paste(sort(powers), collapse = ";")
}

legacy_probabilities_by_group <- function(probabilities, p) {
  if (is.null(probabilities)) {
    return(NULL)
  }
  keys <- vapply(names(probabilities), legacy_probability_group_key,
    character(1), p = p
  )
  values <- setNames(as.numeric(probabilities), keys)
  values[order(names(values))]
}

# Count only compatibility recovery, leaving optimizer evaluations unchanged.
# Cloned function environments avoid modifying the package namespace.
legacy_probability_recovery_spy <- function() {
  calls <- character()
  original_posterior <- log_posteriori_of_perm
  restore <- restore_visited_log_weights
  environment(restore) <- list2env(list(
    log_posteriori_of_perm = function(perm, ...) {
      calls <<- c(calls, as.character(perm))
      original_posterior(perm, ...)
    }
  ), parent = environment(restore_visited_log_weights))
  run <- find_MAP
  environment(run) <- list2env(list(
    restore_visited_log_weights = restore
  ), parent = environment(find_MAP))
  list(run = run, calls = function() calls)
}

test_that("continuing an old MH history reuses its aligned evaluation log", {
  g <- gips(diag(c(1, 2, 4, 7)), 10, was_mean_estimated = FALSE)
  bf <- find_MAP(g, optimizer = "BF", save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )
  full_probabilities <- legacy_probabilities_by_group(
    get_probabilities_from_gips(bf), 4
  )
  first <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 8, save_all_perms = TRUE,
    show_progress_bar = FALSE
  )))
  legacy <- legacy_probability_object(first)
  spy <- legacy_probability_recovery_spy()
  continued <- withr::with_seed(4, suppressMessages(spy$run(legacy,
    optimizer = "continue", max_iter = 8, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))

  visited <- attr(continued, "optimization_info")[["visited_perms"]]
  keys <- sort(unique(vapply(visited, legacy_probability_group_key,
    character(1), p = 4
  )))
  expected <- full_probabilities[keys]
  expected <- expected / sum(expected)
  actual <- legacy_probabilities_by_group(
    get_probabilities_from_gips(continued), 4
  )

  expect_length(spy$calls(), 0L)
  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("old HC histories recover weights with original multi-sample settings", {
  S1 <- matrix(c(
    1.00, 0.20, 0.12, 0.18,
    0.20, 1.10, 0.25, 0.16,
    0.12, 0.25, 0.90, 0.21,
    0.18, 0.16, 0.21, 1.05
  ), nrow = 4, byrow = TRUE)
  S2 <- matrix(c(
    1.20, 0.16, 0.23, 0.11,
    0.16, 0.95, 0.18, 0.20,
    0.23, 0.18, 1.05, 0.15,
    0.11, 0.20, 0.15, 1.10
  ), nrow = 4, byrow = TRUE)
  g <- gips(list(S1, S2), c(5L, 7L), was_mean_estimated = TRUE,
    delta = c(2.7, 4.3),
    D_matrix = list(diag(c(1, 0.9, 1.1, 1.2)), diag(c(1.4, 0.8, 0.9, 1.1)))
  )
  bf <- find_MAP(g, optimizer = "BF", save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )
  full_probabilities <- legacy_probabilities_by_group(
    get_probabilities_from_gips(bf), 4
  )
  hc <- find_MAP(g, optimizer = "HC", max_iter = Inf,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )
  hc_visited <- attr(hc, "optimization_info")[["visited_perms"]]
  expect_gte(length(hc_visited), 3L)
  legacy <- legacy_probability_object(hc)
  spy <- legacy_probability_recovery_spy()

  without_probabilities <- withr::with_seed(1, suppressMessages(spy$run(legacy,
    optimizer = "MH", max_iter = 3, save_all_perms = TRUE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  )))
  expect_length(spy$calls(), 0L)
  expect_null(suppressMessages(get_probabilities_from_gips(without_probabilities)))

  combined <- withr::with_seed(1, suppressMessages(spy$run(legacy,
    optimizer = "MH", max_iter = 3, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))
  visited <- attr(combined, "optimization_info")[["visited_perms"]]
  keys <- sort(unique(vapply(visited, legacy_probability_group_key,
    character(1), p = 4
  )))
  expected <- full_probabilities[keys]
  expected <- expected / sum(expected)
  actual <- legacy_probabilities_by_group(
    get_probabilities_from_gips(combined), 4
  )

  # The new MH segment already supplies the last HC state's weight. Only
  # older HC visits without a matching known weight need a fresh evaluation.
  known_perms <- unique(vapply(tail(visited, 3), as.character, character(1)))
  missing_perms <- setdiff(vapply(hc_visited, as.character, character(1)), known_perms)
  expect_gt(length(missing_perms), 0L)
  expect_setequal(spy$calls(), missing_perms)
  expect_length(spy$calls(), length(missing_perms))
  expect_equal(actual, expected, tolerance = 1e-10)
})
