test_that("Handle improper parameters", {
  expect_error(hill_climbing_optimizer(
    S = matrix_invariant_by_example_perm, number_of_observations = 13,
    max_iter = Inf, start_perm = permutations::id, show_progress_bar = TRUE
  ), "rogress bar is not yet supported for infinite max_iter")

  expect_error(Metropolis_Hastings_optimizer(
    S = matrix_invariant_by_example_perm, number_of_observations = 13, max_iter = Inf,
    start_perm = permutations::permutation("(1,2,3)(4,5)"), show_progress_bar = TRUE
  ), "`max_iter` in `Metropolis_Hastings_optimizer` must")

  g1 <- gips(matrix_invariant_by_example_perm, 13, was_mean_estimated = FALSE)

  expect_error(
    find_MAP(g1, 10, return_probabilities = TRUE, optimizer = "HC"),
    "ilities can only be returned with the `optimizer == 'Metropoli"
  )
  expect_error(
    find_MAP(g1, 10, optimizer = "BD"),
    "imizer\\` must be one of: c\\('MH', 'Me"
  )
  expect_error(
    find_MAP(g1, 10,
      optimizer = c("MH", "BG")
    ),
    "`optimizer` must be the character vector of length 1."
  )
  expect_error(find_MAP(g1, 10, optimizer = "BD"), "You provided `optimizer == 'BD'`")
  expect_error(
    find_MAP(g1, max_iter = NA, optimizer = "MH"),
    "You provided `optimizer == MH` and `max_iter = NA`."
  )
  expect_error(
    find_MAP(g1, Inf, optimizer = "MH"),
    "You provided `max_iter == Inf`."
  )
  expect_error(
    find_MAP(g1, 10, optimizer = "continue"),
    "You provided `optimizer == 'continue'`, but the gips object `g` is not optimized."
  )

  g_small <- gips(matrix(c(1, 0.5, 0.5, 5), ncol = 2), 13, was_mean_estimated = FALSE)
  g_BF <- find_MAP(g_small, optimizer = "brute_force", show_progress_bar = FALSE)
  expect_error(
    find_MAP(g_BF, max_iter = 10, optimizer = "continue"),
    "optimized.*brute_force"
  )

  matrix_big_space <- matrix(runif(35 * 35), nrow = 35)
  matrix_big_space <- crossprod(matrix_big_space)
  g_big <- gips(matrix_big_space, number_of_observations = 13, was_mean_estimated = FALSE)

  expect_error(
    find_MAP(g_big, optimizer = "brute_force", show_progress_bar = FALSE),
    "Optimizer 'brute_force' cannot browse such a big permutational space."
  )
})

bf_final_inputs <- function() {
  S <- diag(c(1, 2))
  list(
    gips(S, 8, was_mean_estimated = FALSE),
    gips(list(S, 2 * S), c(8, 10), was_mean_estimated = FALSE)
  )
}

test_that("find_MAP() rejects every optimizer when the last optimization was BF", {
  optimizers <- c(
    "MH", "Metropolis_Hastings", "HC", "hill_climbing",
    "BF", "brute_force", "full", "continue", NA_character_
  )

  for (g in bf_final_inputs()) {
    g_bf <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
    expect_identical(
      tail(attr(g_bf, "optimization_info")[["optimization_algorithm_used"]], 1L),
      "brute_force"
    )

    withr::with_seed(1, {
      for (optimizer in optimizers) {
        expect_error(suppressMessages(find_MAP(g_bf,
          optimizer = optimizer, max_iter = 3, show_progress_bar = FALSE
        )), "optimized.*brute_force")
      }
    })
  }
})

test_that("completed BF is reported before a missing max_iter", {
  for (g in bf_final_inputs()) {
    g_bf <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)

    # The max_iter error also mentions brute_force; match the reason for rejection.
    expect_error(suppressMessages(find_MAP(g_bf)), "optimized.*brute_force")
    for (optimizer in c("MH", "HC", "continue")) {
      expect_error(find_MAP(g_bf,
        optimizer = optimizer, show_progress_bar = FALSE
      ), "optimized.*brute_force")
    }
  }
})

test_that("forget_perms() does not allow further optimization after BF", {
  for (g in bf_final_inputs()) {
    g_bf <- find_MAP(g, optimizer = "BF", save_all_perms = TRUE,
      show_progress_bar = FALSE
    )
    forgotten <- forget_perms(g_bf)
    expect_true(all(is.na(attr(forgotten, "optimization_info")[["visited_perms"]])))
    expect_identical(
      tail(attr(forgotten, "optimization_info")[["optimization_algorithm_used"]], 1L),
      "brute_force"
    )

    withr::with_seed(1, {
      expect_error(suppressMessages(find_MAP(forgotten,
        optimizer = "MH", max_iter = 3, show_progress_bar = FALSE
      )), "optimized.*brute_force")
    })
  }
})

test_that("missing probabilities do not allow another BF on its result", {
  for (g in bf_final_inputs()) {
    g_bf <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
    expect_null(attr(g_bf, "optimization_info")[["post_probabilities"]])

    expect_error(find_MAP(g_bf, optimizer = "BF", save_all_perms = TRUE,
      return_probabilities = TRUE, show_progress_bar = FALSE
    ), "optimized.*brute_force")

    # The original, unoptimized input can still be used to request probabilities.
    expect_no_error(g_bf_with_probabilities <- find_MAP(g,
      optimizer = "BF", save_all_perms = TRUE,
      return_probabilities = TRUE, show_progress_bar = FALSE
    ))
    probabilities <- get_probabilities_from_gips(g_bf_with_probabilities)
    expect_type(probabilities, "double")
    expect_equal(sum(probabilities), 1)
  }
})

test_that("MH and HC results can be optimized with BF as their final stage", {
  for (g in bf_final_inputs()) {
    for (optimizer in c("MH", "HC")) {
      g_iterative <- withr::with_seed(1, suppressMessages(find_MAP(g,
        optimizer = optimizer, max_iter = 3, show_progress_bar = FALSE
      )))
      expect_no_error(g_bf <- find_MAP(g_iterative,
        optimizer = "BF", show_progress_bar = FALSE
      ))
      expect_identical(
        tail(attr(g_bf, "optimization_info")[["optimization_algorithm_used"]], 1L),
        "brute_force"
      )

      withr::with_seed(1, {
        expect_error(suppressMessages(find_MAP(g_bf,
          optimizer = "MH", max_iter = 3, show_progress_bar = FALSE
        )), "optimized.*brute_force")
      })
    }
  }
})

mh_boundary_fixture <- function(n = 10000) {
  S <- diag(c(1, 100))
  start <- gips(S, n, perm = "(1,2)", was_mean_estimated = FALSE)
  identity <- gips(S, n, perm = "()", was_mean_estimated = FALSE)
  log_weights <- c(log_posteriori_of_gips(start), log_posteriori_of_gips(identity))
  expect_true(all(is.finite(log_weights)))
  expect_gt(log_weights[2], log_weights[1])

  list(start = start, identity = identity, log_weights = log_weights)
}

run_mh_boundary_segment <- function(g, optimizer = "MH",
                                    save_all_perms = TRUE, seed = 1) {
  withr::with_seed(seed, suppressMessages(find_MAP(g,
    optimizer = optimizer, max_iter = 2, save_all_perms = save_all_perms,
    show_progress_bar = FALSE
  )))
}

test_that("MH continuation removes one shared boundary from all saved histories", {
  fixture <- mh_boundary_fixture()
  # At p = 2 there is only one proposal, and the return move is always rejected.
  expect_identical(exp(fixture$log_weights[1] - fixture$log_weights[2]), 0)
  first <- run_mh_boundary_segment(fixture$start)
  continued <- run_mh_boundary_segment(first, optimizer = "continue", seed = 2)
  info <- attr(continued, "optimization_info")

  expect_equal(vapply(info$visited_perms, as.character, character(1)),
    c("(1,2)", "()", "()")
  )
  expect_equal(info$log_posteriori_values, fixture$log_weights[c(1, 2, 2)])
  expect_equal(info$visited_log_posteriori_values, fixture$log_weights[c(1, 2, 2)])
  expect_equal(info$all_n0, c(1, 2, 2))
  expect_equal(sum(info$iterations_performed), 2)
})

test_that("repeated MH continuation keeps genuine repetitions after rejection", {
  fixture <- mh_boundary_fixture()
  first <- run_mh_boundary_segment(fixture$start)
  second <- run_mh_boundary_segment(first, optimizer = "continue", seed = 2)
  # Omitting optimizer must also continue from the last state.
  third <- withr::with_seed(3, suppressMessages(find_MAP(second,
    max_iter = 2, save_all_perms = TRUE, show_progress_bar = FALSE
  )))
  info <- attr(third, "optimization_info")

  expect_equal(vapply(info$visited_perms, as.character, character(1)),
    c("(1,2)", "()", "()", "()")
  )
  expect_equal(info$log_posteriori_values, fixture$log_weights[c(1, 2, 2, 2)])
  expect_equal(info$visited_log_posteriori_values, fixture$log_weights[c(1, 2, 2, 2)])
  expect_equal(info$all_n0, c(1, 2, 2, 2))
  expect_equal(sum(info$iterations_performed), 3)
})

test_that("MH continuation removes boundaries without saving permutations", {
  fixture <- mh_boundary_fixture()
  saved <- run_mh_boundary_segment(fixture$start)
  unsaved <- run_mh_boundary_segment(fixture$start, save_all_perms = FALSE)

  for (stage in 2:3) {
    saved <- run_mh_boundary_segment(saved, optimizer = "continue", seed = stage)
    unsaved <- run_mh_boundary_segment(unsaved, optimizer = "continue",
      save_all_perms = FALSE, seed = stage
    )
    saved_info <- attr(saved, "optimization_info")
    info <- attr(unsaved, "optimization_info")

    expect_equal(info$log_posteriori_values,
      c(fixture$log_weights[1], rep(fixture$log_weights[2], stage))
    )
    expect_equal(info$all_n0, c(1, rep(2, stage)))
    expect_equal(info$log_posteriori_values, saved_info$log_posteriori_values)
    expect_equal(info$all_n0, saved_info$all_n0)
    expect_true(all(is.na(info$visited_perms)))
    expect_null(info$visited_log_posteriori_values)
  }
})

test_that("explicit MH restart retains its start when the best and last states differ", {
  fixture <- mh_boundary_fixture(n = 20)
  return_probability <- exp(fixture$log_weights[1] - fixture$log_weights[2])
  expect_gt(return_probability, 0)
  expect_lt(return_probability, 1)
  first <- run_mh_boundary_segment(fixture$start)
  info <- attr(first, "optimization_info")
  expect_equal(vapply(info$visited_perms, as.character, character(1)),
    c("(1,2)", "()")
  )

  # Construct a possible accepted move from the better identity to the swap.
  # The best model remains identity; this history is not claimed to be sampled.
  fields <- c("visited_perms", "visited_log_posteriori_values",
    "log_posteriori_values", "all_n0")
  for (field in fields) {
    info[[field]] <- rev(info[[field]])
  }
  info$original_perm <- first[[1]]
  info$start_perm <- first[[1]]
  info$last_perm <- fixture$start[[1]]
  info$last_perm_log_posteriori <- fixture$log_weights[1]
  attr(first, "optimization_info") <- info
  expect_false(identical(as.character(first[[1]]), as.character(info$last_perm)))

  continued <- run_mh_boundary_segment(first, optimizer = "continue", seed = 2)
  continued_info <- attr(continued, "optimization_info")
  expect_equal(vapply(continued_info$visited_perms, as.character, character(1)),
    c("()", "(1,2)", "()")
  )
  expect_equal(continued_info$log_posteriori_values, fixture$log_weights[c(2, 1, 2)])
  expect_equal(continued_info$visited_log_posteriori_values,
    fixture$log_weights[c(2, 1, 2)]
  )
  expect_equal(continued_info$all_n0, c(2, 1, 2))

  # Reproduce only the new segment, with the same start and seed as the restart.
  independent <- run_mh_boundary_segment(fixture$identity, seed = 2)
  restarted <- run_mh_boundary_segment(first, optimizer = "MH", seed = 2)
  independent_info <- attr(independent, "optimization_info")
  restarted_info <- attr(restarted, "optimization_info")
  for (field in fields) {
    expect_equal(restarted_info[[field]], c(info[[field]], independent_info[[field]]))
    expect_length(restarted_info[[field]], 4L)
  }
  expect_equal(as.character(restarted_info$visited_perms[[3]]), "()")
})

test_that("MH acceptance rate counts proposals without the initial state", {
  fixture <- mh_boundary_fixture()
  expect_identical(exp(fixture$log_weights[1] - fixture$log_weights[2]), 0)

  for (save in c(FALSE, TRUE)) {
    accepted <- run_mh_boundary_segment(fixture$start, save_all_perms = save)
    rejected <- run_mh_boundary_segment(fixture$identity, save_all_perms = save)
    expect_equal(attr(accepted, "optimization_info")$acceptance_rate, 1)
    expect_equal(attr(rejected, "optimization_info")$acceptance_rate, 0)
  }
})

test_that("forced moves and blocked returns are excluded from MH acceptance rate", {
  run_controlled_mh <- function(n_blocked, save) {
    proposals <- c(
      list(c(1L, 2L), c(1L, 3L), c(2L, 3L)),
      rep(list(c(1L, 2L)), n_blocked),
      list(c(1L, 3L))
    )
    proposal_index <- 0L
    evaluated <- character()
    optimizer <- Metropolis_Hastings_optimizer
    # Override only this copy's dependencies, without changing the namespace
    # or requiring mocking helpers newer than testthat 3.1.5.
    environment(optimizer) <- list2env(list(
      runif_transposition = function(perm_size) {
        proposal_index <<- proposal_index + 1L
        proposals[[proposal_index]]
      },
      log_posteriori_of_perm = function(perm, ...) {
        key <- as.character(perm)
        evaluated <<- c(evaluated, key)
        switch(key,
          "()" = 0,
          "(1,2)" = -1000,
          "(1,3)" = -2000,
          "(2,3)" = -2000,
          -500 # Either 3-cycle.
        )
      }
    ), parent = environment(optimizer))

    # All three initial proposals are rejected; forced-move weights are
    # exactly c(1, 0, 0). The final, unblocked proposal is always accepted.
    result <- withr::with_seed(1, optimizer(
      S = diag(3), number_of_observations = 20,
      start_perm = gips_perm("()", 3), max_iter = length(proposals) + 1L,
      save_all_perms = save, show_progress_bar = FALSE
    ))
    # Forced moves reuse cached weights; blocked returns evaluate nothing.
    expect_equal(evaluated, c("()", "(1,2)", "(1,3)", "(2,3)", last_state))
    expect_equal(attr(result, "optimization_info")$number_of_log_posteriori_calls, length(evaluated))
    expect_equal(attr(result, "optimization_info")$log_posteriori_calls_at_best, 1)
    result
  }

  last_state <- as.character(compose_with_transposition(gips_perm("(1,2)", 3), c(1L, 3L)))
  for (n_blocked in c(0L, 1L, 3L)) {
    for (save in c(FALSE, TRUE)) {
      info <- attr(run_controlled_mh(n_blocked, save), "optimization_info")
      # Returning to identity would improve the posterior, but is forbidden.
      # Check the path so the ratio cannot pass without exercising the moves.
      expect_equal(info$log_posteriori_values,
        c(rep(0, 3), rep(-1000, n_blocked + 1L), -500)
      )
      if (save) {
        expect_equal(vapply(info$visited_perms, as.character, character(1)),
          c(rep("()", 3), rep("(1,2)", n_blocked + 1L), last_state)
        )
      }
      expect_equal(info$acceptance_rate, 1 / 4)
    }
  }
})

test_that("consecutive forced MH moves exclude returns and preserve posterior weights", {
  proposals <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L),
    c(1L, 3L), c(2L, 3L), c(1L, 2L))
  proposal_index <- 0L
  evaluated <- character()
  optimizer <- Metropolis_Hastings_optimizer
  environment(optimizer) <- list2env(list(
    runif_transposition = function(perm_size) {
      proposal_index <<- proposal_index + 1L
      proposals[[proposal_index]]
    },
    log_posteriori_of_perm = function(perm, ...) {
      key <- as.character(perm)
      evaluated <<- c(evaluated, key)
      switch(key, "()" = 0, "(1,2)" = -1000, -2000)
    }
  ), parent = environment(optimizer))
  run <- find_MAP
  environment(run) <- list2env(list(Metropolis_Hastings_optimizer = optimizer),
    parent = environment(run))

  # The first full cache forces (1,2). Two rejected 3-cycles do not yet fill
  # the next cache; only drawing the blocked return triggers the second move.
  # Either cycle may then be chosen, but not identity.
  fit <- withr::with_seed(1, suppressMessages(run(
    gips(diag(3), 20, was_mean_estimated = FALSE),
    optimizer = "MH", max_iter = 7, save_all_perms = TRUE,
    return_probabilities = TRUE, show_progress_bar = FALSE
  )))
  info <- attr(fit, "optimization_info")
  cycles <- vapply(list(c(1L, 3L), c(2L, 3L)), function(e) {
    as.character(compose_with_transposition(gips_perm("(1,2)", 3), e))
  }, character(1))
  expect_equal(evaluated, c("()", "(1,2)", "(1,3)", "(2,3)", cycles))
  expect_equal(info$number_of_log_posteriori_calls, length(evaluated))
  expect_equal(info$log_posteriori_calls_at_best, 1)
  states <- vapply(info$visited_perms, as.character, character(1))
  expect_equal(states[1:6], c(rep("()", 3), rep("(1,2)", 3)))
  expect_true(states[7] %in% cycles)
  expected_weights <- c(rep(0, 3), rep(-1000, 3), -2000)
  expect_equal(info$log_posteriori_values, expected_weights)
  expect_equal(info$visited_log_posteriori_values, expected_weights)

  # Normalize true weights of the three visited groups, not the shifted
  # cache or the technical -Inf used to prevent the return to identity.
  expected_probabilities <- c("()" = 1, "(1,2)" = 0, "(1,2,3)" = 0)
  actual <- get_probabilities_from_gips(fit)
  expect_setequal(names(actual), names(expected_probabilities))
  expect_equal(actual[names(expected_probabilities)], expected_probabilities)
})

# Isolate objective-call accounting from numerical posterior calculations.
posterior_counter_fixture <- function(transposition = c(1L, 2L),
                                      identity_log_weight = 0, other_log_weight = -1000) {
  evaluated <- numeric()
  goal <- function(perm, ...) {
    value <- if (as.character(perm) == "()") identity_log_weight else other_log_weight
    evaluated <<- c(evaluated, value)
    value
  }
  optimizers <- list(
    Metropolis_Hastings_optimizer = Metropolis_Hastings_optimizer,
    hill_climbing_optimizer = hill_climbing_optimizer,
    brute_force_optimizer = brute_force_optimizer
  )
  optimizers <- lapply(optimizers, function(optimizer) {
    environment(optimizer) <- list2env(list(
      log_posteriori_of_perm = goal,
      runif_transposition = function(perm_size) transposition
    ), parent = environment(optimizer))
    optimizer
  })
  find <- find_MAP
  environment(find) <- list2env(optimizers, parent = environment(find))
  list(
    run = function(g, optimizer, max_iter = 6, save = FALSE) {
      withr::with_seed(1, suppressMessages(find(g,
        optimizer = optimizer, max_iter = max_iter,
        save_all_perms = save, show_progress_bar = FALSE
      )))
    },
    values = function() evaluated
  )
}

test_that("each optimizer records actual objective calls and the first best evaluation", {
  for (save in c(FALSE, TRUE)) {
    for (optimizer in c("MH", "HC", "BF")) {
      fixture <- posterior_counter_fixture()
      fit <- fixture$run(gips(diag(3), 20, was_mean_estimated = FALSE),
        optimizer, save = save
      )
      info <- attr(fit, "optimization_info")
      # MH repeatedly proposes one rejected neighbor, never filling the cache.
      # HC checks all three neighbors; BF checks all five cyclic groups.
      expect_length(fixture$values(), c(MH = 2L, HC = 4L, BF = 5L)[[optimizer]])
      expect_equal(info$number_of_log_posteriori_calls, length(fixture$values()))
      expect_equal(info$log_posteriori_calls_at_best, which.max(fixture$values()))
      if (optimizer == "MH") expect_length(info$log_posteriori_values, 6L)
    }
  }
})

test_that("MH continuation and restart count reevaluated boundaries and retain earlier best", {
  for (save in c(FALSE, TRUE)) {
    for (optimizer in c("continue", "MH")) {
      fixture <- posterior_counter_fixture()
      first <- fixture$run(gips(diag(3), 20, was_mean_estimated = FALSE),
        "MH", save = save
      )
      combined <- fixture$run(first, optimizer, max_iter = 4, save = save)
      info <- attr(combined, "optimization_info")
      expect_equal(fixture$values(), c(0, -1000, 0, -1000))
      expect_equal(info$number_of_log_posteriori_calls, 4)
      expect_equal(info$log_posteriori_calls_at_best, 1)
      expect_length(info$log_posteriori_values, if (optimizer == "continue") 9L else 10L)
    }
  }
})

test_that("the best evaluation counter distinguishes small strict improvements", {
  fixture <- posterior_counter_fixture(identity_log_weight = 5e-8, other_log_weight = 0)
  fit <- fixture$run(gips(diag(3), 20, perm = "(1,2)", was_mean_estimated = FALSE),
    "MH", max_iter = 2
  )
  expect_equal(fixture$values(), c(0, 5e-8))
  expect_equal(attr(fit, "optimization_info")$number_of_log_posteriori_calls, 2)
  expect_equal(attr(fit, "optimization_info")$log_posteriori_calls_at_best, 2)
})

test_that("mixed histories offset a later best by actual calls and preserve it on ties", {
  fixture <- posterior_counter_fixture(c(1L, 3L))
  # MH starts at a swap and accepts an equally weighted 3-cycle.
  # HC starts from the saved best swap and discovers identity as its first neighbor.
  mh <- fixture$run(gips(diag(3), 20, perm = "(1,2)", was_mean_estimated = FALSE),
    "MH", max_iter = 2, save = TRUE
  )
  expect_equal(attr(mh, "optimization_info")$log_posteriori_calls_at_best, 1)
  mixed <- fixture$run(mh, "HC", save = TRUE)
  info <- attr(mixed, "optimization_info")
  expect_length(fixture$values(), 9L) # Two MH calls, seven HC calls.
  expect_equal(info$number_of_log_posteriori_calls, 9)
  expect_equal(info$log_posteriori_calls_at_best, 4)
  expect_equal(info$log_posteriori_calls_at_best, which.max(fixture$values()))

  restarted <- fixture$run(mixed, "MH", save = TRUE)
  info <- attr(restarted, "optimization_info")
  expect_length(fixture$values(), 11L)
  expect_equal(info$number_of_log_posteriori_calls, 11)
  expect_equal(info$log_posteriori_calls_at_best, 4)
})

test_that("MH acceptance rates combine by proposal counts across unequal stages", {
  fixture <- mh_boundary_fixture()
  for (save in c(FALSE, TRUE)) {
    first <- run_mh_boundary_segment(fixture$start, save_all_perms = save)
    for (optimizer in c("continue", "MH")) {
      combined <- withr::with_seed(2, suppressMessages(find_MAP(first,
        optimizer = optimizer, max_iter = 4,
        save_all_perms = save, show_progress_bar = FALSE
      )))
      # One accepted proposal in the first stage, three rejections in the next.
      expect_equal(attr(combined, "optimization_info")$acceptance_rate, 1 / 4)
    }
  }
})

test_that("HC and BF have no acceptance rate and HC does not affect MH rates", {
  fixture <- mh_boundary_fixture()
  hc <- find_MAP(fixture$start, optimizer = "HC", max_iter = Inf,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )
  bf <- find_MAP(fixture$start, optimizer = "BF", show_progress_bar = FALSE)
  hc_continued <- find_MAP(hc, optimizer = "continue", max_iter = 2,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )
  for (g in list(hc, hc_continued, bf)) {
    expect_null(attr(g, "optimization_info")$acceptance_rate)
    result <- summary(g)
    expect_null(result$acceptance_rate)
    expect_false(any(grepl("Acceptance rate", capture.output(print(result)), fixed = TRUE)))
  }

  hc_mh <- run_mh_boundary_segment(hc)
  expect_equal(attr(hc_mh, "optimization_info")$acceptance_rate, 0)
  expect_equal(summary(hc_mh)$acceptance_rate, 0)
  expect_output(print(summary(hc_mh)), "Acceptance rate (MH proposals):", fixed = TRUE)

  mh <- run_mh_boundary_segment(fixture$start)
  mh_hc <- find_MAP(mh, optimizer = "HC", max_iter = Inf,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )
  expect_equal(attr(mh_hc, "optimization_info")$acceptance_rate, 1)
  mh_hc_mh <- withr::with_seed(2, suppressMessages(find_MAP(mh_hc,
    optimizer = "MH", max_iter = 4,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )))
  expect_equal(attr(mh_hc_mh, "optimization_info")$acceptance_rate, 1 / 4)
})

test_that("removing MH continuation boundaries preserves acceptance rate aggregation", {
  fixture <- mh_boundary_fixture()
  first <- run_mh_boundary_segment(fixture$start)
  independent_second <- run_mh_boundary_segment(fixture$identity, seed = 2)
  independent_third <- run_mh_boundary_segment(fixture$identity, seed = 3)
  independent_rates <- vapply(list(first, independent_second, independent_third),
    function(g) attr(g, "optimization_info")[["acceptance_rate"]], numeric(1)
  )
  expect_gt(independent_rates[1], independent_rates[2])
  expect_equal(independent_rates[2], independent_rates[3])

  second <- run_mh_boundary_segment(first, optimizer = "continue", seed = 2)
  third <- run_mh_boundary_segment(second, optimizer = "continue", seed = 3)
  # Equal-length segments must keep their equal weights when boundaries are removed.
  # Derive rates from separate runs instead of fixing the current per-run convention.
  expect_equal(attr(second, "optimization_info")[["acceptance_rate"]],
    mean(independent_rates[1:2])
  )
  expect_equal(attr(third, "optimization_info")[["acceptance_rate"]],
    mean(independent_rates)
  )
})

test_that("MH continuation after HC removes only the shared boundary", {
  for (save in c(FALSE, TRUE)) {
    g <- find_MAP(mh_boundary_fixture()$start,
      optimizer = "HC", max_iter = Inf,
      save_all_perms = save, show_progress_bar = FALSE
    )
    g <- run_mh_boundary_segment(g, save_all_perms = save)

    fields <- c("all_n0", "log_posteriori_values")
    if (save) {
      fields <- c(fields, "visited_perms", "visited_log_posteriori_values")
    }

    for (optimizer in c("continue", NA_character_)) {
      before <- attr(g, "optimization_info")
      g <- run_mh_boundary_segment(g,
        optimizer = optimizer, save_all_perms = save
      )
      after <- attr(g, "optimization_info")

      # One MH proposal retains one new record after rejection.
      for (field in fields) {
        expect_equal(after[[field]],
          c(before[[field]], tail(before[[field]], 1L))
        )
      }
    }
  }
})

test_that("explicit MH restart after HC retains its initial state", {
  g <- find_MAP(mh_boundary_fixture()$start,
    optimizer = "HC", max_iter = Inf,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )
  g <- run_mh_boundary_segment(g)
  before <- attr(g, "optimization_info")

  restarted <- run_mh_boundary_segment(g, optimizer = "MH")
  after <- attr(restarted, "optimization_info")

  # A restart records its initial state and the state after rejection.
  for (field in c("all_n0", "log_posteriori_values",
                 "visited_perms", "visited_log_posteriori_values")) {
    last <- tail(before[[field]], 1L)
    expect_equal(after[[field]], c(before[[field]], last, last))
  }
})

test_that("HC records n0 for its initial and accepted states without requiring saved permutations", {
  for (was_mean_estimated in c(FALSE, TRUE)) {
    for (start_perm in c("(1,2)", "()")) {
      g <- gips(diag(c(1, 100)), 10000,
        perm = start_perm, was_mean_estimated = was_mean_estimated
      )
      saved <- find_MAP(g, optimizer = "HC", max_iter = Inf,
        save_all_perms = TRUE, show_progress_bar = FALSE
      )
      unsaved <- find_MAP(g, optimizer = "HC", max_iter = Inf,
        save_all_perms = FALSE, show_progress_bar = FALSE
      )
      info <- attr(saved, "optimization_info")
      expected <- if (start_perm == "(1,2)") c(1, 2) else 2

      # The final unsuccessful HC step must not add a repeated state.
      expect_equal(info$all_n0, expected + was_mean_estimated)
      expect_equal(info$all_n0,
        vapply(info$visited_perms, length, integer(1)) + was_mean_estimated
      )
      expect_length(info$all_n0, length(info$visited_log_posteriori_values))
      expect_equal(attr(unsaved, "optimization_info")$all_n0, info$all_n0)
      expect_equal(attr(forget_perms(saved), "optimization_info")$all_n0, info$all_n0)
    }
  }
})

test_that("BF records n0 in evaluation order with or without saved permutations", {
  # p = 2 enumerates every permutation; p = 3 uses prepared group representatives.
  for (p in c(2L, 3L)) {
    for (was_mean_estimated in c(FALSE, TRUE)) {
      g <- gips(diag(seq_len(p)), 20, was_mean_estimated = was_mean_estimated)
      saved <- find_MAP(g, optimizer = "BF",
        save_all_perms = TRUE, show_progress_bar = FALSE
      )
      unsaved <- find_MAP(g, optimizer = "BF",
        save_all_perms = FALSE, show_progress_bar = FALSE
      )
      info <- attr(saved, "optimization_info")
      # BF stores permutations::cycle objects, which omit fixed points.
      cycle_counts <- vapply(seq_along(info$visited_perms), function(i) {
        length(gips_perm(info$visited_perms[i], p))
      }, integer(1))

      expect_equal(info$all_n0, cycle_counts + was_mean_estimated)
      expect_equal(range(info$all_n0), c(1, p) + was_mean_estimated)
      expect_length(info$all_n0, length(info$log_posteriori_values))
      expect_equal(attr(unsaved, "optimization_info")$all_n0, info$all_n0)
      expect_equal(attr(forget_perms(saved), "optimization_info")$all_n0, info$all_n0)
    }
  }
})

test_that("mixed HC and MH histories retain n0 for every saved state", {
  for (was_mean_estimated in c(FALSE, TRUE)) {
    g <- gips(diag(c(1, 100)), 10000,
      perm = "(1,2)", was_mean_estimated = was_mean_estimated
    )
    hc <- find_MAP(g, optimizer = "HC", max_iter = Inf,
      save_all_perms = TRUE, show_progress_bar = FALSE
    )
    hc_mh <- run_mh_boundary_segment(hc)
    mh <- run_mh_boundary_segment(g)
    mh_hc <- find_MAP(mh, optimizer = "HC", max_iter = Inf,
      save_all_perms = TRUE, show_progress_bar = FALSE
    )

    expect_equal(attr(hc_mh, "optimization_info")$all_n0,
      c(1, 2, 2, 2) + was_mean_estimated
    )
    expect_equal(attr(mh_hc, "optimization_info")$all_n0,
      c(1, 2, 2) + was_mean_estimated
    )
    for (result in list(hc_mh, mh_hc)) {
      info <- attr(result, "optimization_info")
      expect_false(anyNA(info$all_n0))
      expect_equal(info$all_n0,
        vapply(info$visited_perms, length, integer(1)) + was_mean_estimated
      )
      expect_length(info$all_n0, length(info$visited_log_posteriori_values))
    }
  }
})

test_that("Handle proper parameters", {
  g1 <- gips(matrix_invariant_by_example_perm, 13, was_mean_estimated = FALSE)
  
  expect_silent(
    g_map <- find_MAP(g1, 10, optimizer = "MH", show_progress_bar = FALSE)
  )

  expect_silent(out <- Metropolis_Hastings_optimizer(
    S = matrix_invariant_by_example_perm,
    number_of_observations = 13, max_iter = 2,
    start_perm = permutations::id,
    show_progress_bar = FALSE
  ))
  expect_output(
    out <- Metropolis_Hastings_optimizer(
      S = matrix_invariant_by_example_perm,
      number_of_observations = 13, max_iter = 3,
      start_perm = permutations::permutation("(1,2,3)(4,5)"),
      show_progress_bar = TRUE
    ),
    "===="
  )
  expect_output(
    out <- hill_climbing_optimizer(
      S = matrix_invariant_by_example_perm,
      number_of_observations = 13, max_iter = 8,
      start_perm = permutations::id,
      show_progress_bar = TRUE
    ),
    "===="
  )
  expect_warning(expect_output(
    out <- hill_climbing_optimizer(
      S = matrix_invariant_by_example_perm,
      number_of_observations = 13, max_iter = 2,
      start_perm = permutations::id,
      show_progress_bar = TRUE
    ),
    "===="
  )) # Can we stack the `expect`s like that? Looks neat!

  g_small <- gips(matrix(c(1, 0.5, 0.5, 5), ncol = 2), 13, was_mean_estimated = FALSE)
  expect_message(
    find_MAP(g_small, max_iter = 10, optimizer = "MH", show_progress_bar = FALSE),
    "Consider using `optimizer = 'brute_force'`, because it will check only 2 permutations and will definitely find the maximum a posteriori estimator."
  )

  expect_silent(g_MAP_small <- find_MAP(g_small,
    max_iter = 10,
    optimizer = "HC",
    show_progress_bar = FALSE,
    save_all_perms = TRUE
  ))
  expect_message(find_MAP(g_MAP_small,
    max_iter = 10,
    optimizer = "MH",
    show_progress_bar = FALSE,
    save_all_perms = TRUE
  ))

  expect_warning(
    g_limited <- find_MAP(g1, max_iter = 2, optimizer = "HC",
      save_all_perms = TRUE, show_progress_bar = FALSE
    ),
    "Hill Climbing algorithm did not converge in 2 iterations!"
  )
  limited_info <- attr(g_limited, "optimization_info")
  expect_length(limited_info$all_n0, 3L)
  expect_equal(limited_info$all_n0,
    vapply(limited_info$visited_perms, length, integer(1))
  )
  expect_length(limited_info$all_n0, length(limited_info$visited_log_posteriori_values))
  
  
  
  
  g2 <- gips(matrix_invariant_by_example_perm[1:4, 1:4], 13, was_mean_estimated = FALSE)
  expect_message(
    g_map <- find_MAP(g2, show_progress_bar = FALSE),
    "The 'optimizer = NA' was automatically changed to 'optimizer = \"BF\"."
  )
})

test_that("Warns when found group has n0 > n", {
  number_of_observations <- 2
  Z <- MASS::mvrnorm(number_of_observations,
    mu = numeric(6),
    Sigma = matrix_invariant_by_example_perm
  )
  S <- crossprod(Z) / number_of_observations

  g <- gips(S, number_of_observations, was_mean_estimated = FALSE)

  expect_warning(
    find_MAP(g,
      max_iter = 2, show_progress_bar = FALSE,
      optimizer = "MH"
    ), # after 2 iterations the n0 can either be 6, or 5, or 4. It cannot be 1 or 2.
    "n0"
  )
})

test_that("find_MAP() can guess the correct optimizer and message the user", {
  g <- gips(matrix_invariant_by_example_perm, 13, was_mean_estimated = FALSE)

  # can guess:
  expect_message(
    find_MAP(g,
      max_iter = 2, show_progress_bar = FALSE,
      optimizer = "Me"
    ),
    "is will be changed to `optimizer == 'Metropolis_Hastings'`"
  )

  # cannot guess, can be "MH" or "Metropolis_Hastings":
  expect_error(
    find_MAP(g,
      max_iter = 2, show_progress_bar = FALSE,
      optimizer = "M"
    ),
    "You provided `optimizer == 'M'`"
  )
  # cannot guess
  expect_error(
    find_MAP(g,
      max_iter = 2, show_progress_bar = FALSE,
      optimizer = "etropolis_Hastings"
    ),
    "You provided `optimizer == 'etropolis_Hastings'`"
  )
})

test_that("find_MAP() will remember the right number of observations and was_mean_estimated", {
  number_of_observations <- 13

  # em - estimated mean
  g_em_FALSE <- gips(
    matrix_invariant_by_example_perm,
    number_of_observations,
    was_mean_estimated = FALSE
  )
  g_em_TRUE <- gips(
    matrix_invariant_by_example_perm,
    number_of_observations,
    was_mean_estimated = TRUE
  )

  g_map_em_FALSE <- find_MAP(
    g_em_FALSE,
    max_iter = 10,
    show_progress_bar = FALSE, optimizer = "MH"
  )
  g_map_em_TRUE <- find_MAP(
    g_em_TRUE,
    max_iter = 10,
    show_progress_bar = FALSE, optimizer = "MH"
  )

  expect_equal(
    attr(g_map_em_FALSE, "number_of_observations"),
    number_of_observations
  )
  expect_equal(
    attr(g_map_em_TRUE, "number_of_observations"),
    number_of_observations
  )
})

test_that("find_MAP() with calculate exact probabilities will return probability", {
  g <- gips(
    S = matrix_invariant_by_example_perm[1:4, 1:4],
    number_of_observations = 13,
    D_matrix = diag(1, 4)
  )
  g_map <- find_MAP(g,
    max_iter = 10, show_progress_bar = FALSE,
    optimizer = "brute_force", return_probabilities = TRUE, save_all_perms = TRUE
  )

  # note: Tied probabilities have no guaranteed relative order.
  my_post_prob <- c(
    `(1,2,3,4)` = 0.326304997567465, `(1,2,4,3)` = 0.326304997567465,
    `(1,3,2,4)` = 0.326304997567465, `(1,2,3)` = 0.00709444735420907,
    `(2,3,4)` = 0.00466132503808372, `(1,2,4)` = 0.00466132503808372,
    `(1,3,4)` = 0.00466132503808372, `(1,2)(3,4)` = 2.13328433453314e-06,
    `(1,3)(2,4)` = 2.13328433453314e-06, `(1,4)(2,3)` = 2.13328433453314e-06,
    `(2,3)` = 3.5989590475781e-08, `(1,2)` = 3.5989590475781e-08,
    `(1,3)` = 3.5989590475781e-08, `(3,4)` = 2.56690835859916e-08,
    `(2,4)` = 2.56690835859916e-08, `(1,4)` = 2.56690835859916e-08,
    `()` = 1.18881017985689e-13
  )
  
  post_probabilities <- attr(g_map, "optimization_info")[["post_probabilities"]]
  
  expect_length(post_probabilities, length(my_post_prob))
  expect_setequal(names(post_probabilities), names(my_post_prob))
  expect_equal(
    post_probabilities[names(my_post_prob)],
    my_post_prob
  )
  expect_true(all(diff(unname(post_probabilities)) <= 0))
  expect_equal(sum(post_probabilities), 1)
})

test_that("BF keeps deduplication when its permutation count differs from OEIS", {
  optimizer <- brute_force_optimizer
  # Exercise the full-permutation path in a small space, without a p = 11 run.
  environment(optimizer) <- list2env(list(
    perm_group_generators_list = list(rep(TRUE, factorial(3)))
  ), parent = environment(optimizer))
  arguments <- list(
    S = diag(c(1, 2, 3)), number_of_observations = 20,
    return_probabilities = TRUE, save_all_perms = TRUE,
    show_progress_bar = TRUE
  )

  expected_output <- capture.output(
    expected_info <- attr(do.call(brute_force_optimizer, arguments), "optimization_info")
  )
  actual_output <- capture.output(
    actual_info <- attr(do.call(optimizer, arguments), "optimization_info")
  )
  expect_length(expected_output, 2)
  expect_length(actual_output, 3)
  expect_length(expected_info$visited_perms, 5)
  expect_length(actual_info$visited_perms, 6)
  expected <- expected_info$post_probabilities
  actual <- actual_info$post_probabilities

  expect_length(actual, 5)
  expect_setequal(names(actual), names(expected))
  expect_equal(actual[names(expected)], expected)
})

test_that("there is proper number of generators", {
  num_of_generators <- sapply(perm_group_generators_list, sum)
  expect_equal(length(perm_group_generators_list), 8)
  for (i in seq_along(num_of_generators) + 2) {
    expect_equal(num_of_generators[i - 2], OEIS_A051625[i])
  }
})

test_that("permutations package preserves the allperms order", {
  all_perms <- permutations::as.cycle(permutations::allperms(4))
  expect_identical(as.character(all_perms[2]), "(3,4)")
  expect_identical(as.character(all_perms[14]), "(1,3,4,2)")
  expect_identical(as.character(all_perms[17]), "(1,3)(2,4)")
  
  all_perms <- permutations::as.cycle(permutations::allperms(5))
  expect_identical(as.character(all_perms[15]), "(2,4)")
  expect_identical(as.character(all_perms[58]), "(1,3,4,5)")
  expect_identical(as.character(all_perms[101]), "(1,5,3,4,2)")
  expect_identical(as.character(all_perms[102]), "(1,5,2)(3,4)")
})
