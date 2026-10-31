# Identify the entire cyclic subgroup independently of the production
# representative-selection helper. These fixtures have only four variables.
history_probability_group_key <- function(perm, p) {
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

history_probability_keys <- function(perms, p) {
  unique(vapply(perms, history_probability_group_key, character(1), p = p))
}

history_probabilities_by_group <- function(probabilities, p) {
  if (is.null(probabilities)) {
    return(NULL)
  }
  keys <- vapply(names(probabilities), history_probability_group_key,
    character(1), p = p
  )
  values <- setNames(as.numeric(probabilities), keys)
  values[order(names(values))]
}

test_that("MH continuation uses posterior weights over both saved histories", {
  S <- diag(c(1, 2, 4, 7))
  g <- gips(S, 10, was_mean_estimated = FALSE)
  bf <- find_MAP(g, optimizer = "BF", save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )
  full_probabilities <- history_probabilities_by_group(
    get_probabilities_from_gips(bf), 4
  )

  first <- withr::with_seed(1, suppressMessages(find_MAP(g,
    optimizer = "MH", max_iter = 8, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))
  first_info <- attr(first, "optimization_info")

  # Independently obtain the second segment, so the reference union does not
  # depend on whether continuation retains or removes its repeated boundary.
  restart <- gips(S, 10, perm = first_info[["last_perm"]],
    was_mean_estimated = FALSE
  )
  second <- withr::with_seed(4, suppressMessages(find_MAP(restart,
    optimizer = "MH", max_iter = 8, save_all_perms = TRUE,
    show_progress_bar = FALSE
  )))
  continued <- withr::with_seed(4, suppressMessages(find_MAP(first,
    optimizer = "continue", max_iter = 8, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))

  first_keys <- history_probability_keys(first_info[["visited_perms"]], 4)
  second_keys <- history_probability_keys(
    attr(second, "optimization_info")[["visited_perms"]], 4
  )
  visited_keys <- sort(union(first_keys, second_keys))
  expected <- full_probabilities[visited_keys]
  expected <- expected / sum(expected)

  # Both segments contribute distinct groups, and the reference is a proper
  # subset of the BF posterior, not a full-posterior comparison in disguise.
  expect_gt(length(setdiff(second_keys, first_keys)), 0L)
  expect_gt(sum(expected[setdiff(first_keys, second_keys)]), 0.001)
  expect_lt(sum(full_probabilities[visited_keys]), 0.99)

  actual <- history_probabilities_by_group(
    get_probabilities_from_gips(continued), 4
  )
  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("MH after several HC improvements includes HC-only groups", {
  S <- matrix(c(
    1.00, 0.20, 0.12, 0.18,
    0.20, 1.10, 0.25, 0.16,
    0.12, 0.25, 0.90, 0.21,
    0.18, 0.16, 0.21, 1.05
  ), nrow = 4, byrow = TRUE)
  g <- gips(S, 4, was_mean_estimated = FALSE)
  bf <- find_MAP(g, optimizer = "BF", save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )
  full_probabilities <- history_probabilities_by_group(
    get_probabilities_from_gips(bf), 4
  )

  hc <- find_MAP(g, optimizer = "HC", max_iter = Inf,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )
  hc_visited <- attr(hc, "optimization_info")[["visited_perms"]]
  hc_log_weights <- vapply(hc_visited, function(perm) {
    log_posteriori_of_gips(gips(S, 4, perm = perm,
      was_mean_estimated = FALSE
    ))
  }, numeric(1))
  expect_gte(length(hc_visited), 3L)
  expect_true(all(diff(hc_log_weights) > 1e-6))

  restart <- gips(S, 4, perm = hc, was_mean_estimated = FALSE)
  mh_only <- withr::with_seed(1, suppressMessages(find_MAP(restart,
    optimizer = "MH", max_iter = 3, save_all_perms = TRUE,
    show_progress_bar = FALSE
  )))
  combined <- withr::with_seed(1, suppressMessages(find_MAP(hc,
    optimizer = "MH", max_iter = 3, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))

  hc_keys <- history_probability_keys(hc_visited, 4)
  mh_keys <- history_probability_keys(
    attr(mh_only, "optimization_info")[["visited_perms"]], 4
  )
  visited_keys <- sort(union(hc_keys, mh_keys))

  # From identity, HC evaluates all six transpositions in its first step.
  # Evaluating a candidate must not add its group to the visited support.
  expect_identical(history_probability_group_key(hc_visited[[1]], 4),
    history_probability_group_key("()", 4)
  )
  first_step_neighbors <- utils::combn(seq_len(4), 2, function(pair) {
    paste0("(", paste(pair, collapse = ","), ")")
  }, simplify = FALSE)
  evaluated_only_keys <- setdiff(
    history_probability_keys(first_step_neighbors, 4), visited_keys
  )
  expect_gt(length(evaluated_only_keys), 0L)

  expected <- full_probabilities[visited_keys]
  expected <- expected / sum(expected)

  # A calculation from the MH segment alone must give a materially different
  # answer: earlier HC groups retain substantial posterior mass in this union.
  expect_gt(sum(expected[setdiff(hc_keys, mh_keys)]), 0.1)
  expect_lt(sum(full_probabilities[visited_keys]), 0.99)

  actual <- history_probabilities_by_group(
    suppressMessages(get_probabilities_from_gips(combined)), 4
  )
  expect_length(intersect(names(actual), evaluated_only_keys), 0L)
  expect_equal(actual, expected, tolerance = 1e-10)
})
