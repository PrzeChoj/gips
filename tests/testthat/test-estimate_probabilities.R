probability_test_perms <- function(x) {
  lapply(x, gips_perm, 3)
}

probability_test_calculate <- function(perms, log_weights) {
  calculate_visited_probabilities(perms, log_weights)
}

test_that("supplied posterior weights are normalized over distinct groups", {
  perms <- probability_test_perms(c("()", "(1,2)", "(1,2,3)"))

  probabilities <- probability_test_calculate(
    perms, log_weights = log(c(2, 3, 5))
  )

  expect_equal(probabilities, c("(1,2,3)" = 0.5, "(1,2)" = 0.3, "()" = 0.2))
})

test_that("repeated visits do not change weight-based probabilities", {
  perms <- probability_test_perms(c("()", "(1,2)", "(1,2,3)"))
  log_weights <- log(c(2, 3, 5))
  visits <- c(rep(1L, 9), 2L, 3L, 3L)

  probabilities <- probability_test_calculate(
    perms[visits], log_weights = log_weights[visits]
  )

  expect_equal(probabilities, c("(1,2,3)" = 0.5, "(1,2)" = 0.3, "()" = 0.2))
})

test_that("equivalent generators contribute one group weight without a phi correction", {
  perms <- probability_test_perms(c("(1,2,3)", "(1,3,2)", "(1,2)"))

  probabilities <- probability_test_calculate(
    perms, log_weights = log(c(6, 6, 3))
  )

  expect_equal(probabilities, c("(1,2,3)" = 2 / 3, "(1,2)" = 1 / 3))
})

test_that("different groups with equal weights remain distinct", {
  perms <- probability_test_perms(c("(1,2)", "(1,3)", "(2,3)"))
  expected <- c("(1,2)" = 1 / 3, "(1,3)" = 1 / 3, "(2,3)" = 1 / 3)

  probabilities <- probability_test_calculate(
    perms, log_weights = rep(log(7), 3)
  )

  expect_setequal(names(probabilities), names(expected))
  expect_equal(probabilities[names(expected)], expected)
})

test_that("a single supplied group has conditional probability one", {
  probabilities <- probability_test_calculate(
    probability_test_perms("(1,2,3)"), log_weights = -123.4
  )

  expect_equal(probabilities, c("(1,2,3)" = 1))
})

test_that("probabilities are stable under large common shifts of log weights", {
  perms <- probability_test_perms(c("()", "(1,2)", "(1,2,3)"))
  expected <- c("()" = 1, "(1,2)" = exp(-1), "(1,2,3)" = exp(-2))
  expected <- expected / sum(expected)

  # Direct exponentiation overflows at +1000 and underflows at -1000.
  for (shift in c(-1000, 0, 1000)) {
    probabilities <- probability_test_calculate(
      perms, log_weights = c(0, -1, -2) + shift
    )

    expect_equal(probabilities, expected, tolerance = 1e-12,
                 info = paste("shift =", shift))
  }
})

test_that("widely separated finite log weights produce finite probabilities", {
  perms <- probability_test_perms(c("()", "(1,2)", "(1,2,3)"))
  expected <- c("()" = 1, "(1,2)" = 0, "(1,2,3)" = 0)

  probabilities <- probability_test_calculate(
    perms, log_weights = c(1000, -1000, -500)
  )

  expect_setequal(names(probabilities), names(expected))
  expect_equal(probabilities[names(expected)], expected)
})

test_that("visited probability inputs must be valid and aligned", {
  perms <- probability_test_perms(c("()", "(1,2)"))

  expect_error(calculate_visited_probabilities(list(), numeric()), "nonempty list")
  expect_error(calculate_visited_probabilities("()", 0), "nonempty list")
  expect_error(calculate_visited_probabilities(list("()"), 0), "gips_perm")
  expect_error(calculate_visited_probabilities(
    list(structure("()", class = "gips_perm", size = 3)), 0
  ), "gips_perm.*list")
  for (invalid_size in list(NULL, c(3, 3), NA_real_, Inf, "3", 3 + 1i, 0, 3.5)) {
    invalid_perm <- perms[[1]]
    attr(invalid_perm, "size") <- invalid_size
    expect_error(calculate_visited_probabilities(list(invalid_perm), 0),
                 "finite positive integer.*size")
  }
  expect_error(calculate_visited_probabilities(
    list(gips_perm("()", 3), gips_perm("()", 4)), c(0, 0)
  ), "same size")
  expect_error(calculate_visited_probabilities(
    list(new_gips_perm(list(c(1L, 2L)), 3)), 0
  ), "complete permutation")

  expect_error(calculate_visited_probabilities(perms, 0), "log_weights")
  expect_error(calculate_visited_probabilities(perms, c("0", "1")), "log_weights")
  expect_error(calculate_visited_probabilities(perms, c(0, 1i)), "log_weights")
  expect_error(calculate_visited_probabilities(perms, matrix(c(0, 1))), "log_weights")
  for (nonfinite in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(calculate_visited_probabilities(perms, c(0, nonfinite)), "log_weights")
  }
})

test_that("exact permutation keys and group deduplication ignore display labels", {
  perms <- probability_test_perms(c("(1,2)", "(1,3)", "(2,3)"))
  withr::with_options(list(perm_set = rep("x", 3), print_in_length_order = TRUE), {
    expect_identical(vapply(perms, get_permutation_key, character(1)),
                     c("2,1,3", "3,2,1", "1,3,2"))
    probabilities <- calculate_visited_probabilities(perms, log(c(2, 3, 5)))
    expect_equal(unname(probabilities), c(0.5, 0.3, 0.2))
  })
  expect_identical(get_permutation_key(gips_perm("()", 3)), "1,2,3")
  expect_identical(get_permutation_key(gips_perm("()", 4)), "1,2,3,4")
})

test_that("group representatives are calculated once per distinct permutation", {
  perms <- probability_test_perms(c("(1,2,3)", "(1,3,2)", "(1,2)"))
  visits <- c(rep(1L, 20), rep(2L, 10), rep(3L, 5))
  log_weights <- log(c(6, 6, 3))
  representative_calls <- 0L
  original_representative <- get_group_representative
  calculate_with_counter <- calculate_visited_probabilities
  environment(calculate_with_counter) <- list2env(list(
    get_group_representative = function(perm) {
      representative_calls <<- representative_calls + 1L
      original_representative(perm)
    }
  ), parent = environment(calculate_visited_probabilities))

  probabilities <- calculate_with_counter(perms[visits], log_weights[visits])

  expect_identical(representative_calls, 3L)
  expect_equal(probabilities, c("(1,2,3)" = 2 / 3, "(1,2)" = 1 / 3))
})

test_that("get_representative examples", {
  p1 <- gips_perm(permutations::as.word(c(2, 1, 3, 6, 4, 5)), 6)
  expected_p <- gips_perm("(1,2)(4,5,6)", 6)

  expect_equal(
    get_group_representative(p1),
    expected_p
  )

  gips_nullperm <- gips_perm(permutations::nullword, 0)

  expect_equal(
    get_group_representative(gips_nullperm),
    gips_nullperm
  )
})
