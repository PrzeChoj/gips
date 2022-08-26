test_examples("../..") # example for Metropolis_Hastings and hill_climbing are here

test_that("Handle improper parameters", {
  expect_error(hill_climbing(
    S = matrix_invariant_by_example_perm, number_of_observations = 13,
    max_iter = Inf, show_progress_bar = TRUE
  ))

  expect_error(Metropolis_Hastings(
    S = matrix_invariant_by_example_perm, number_of_observations = 13, max_iter = Inf,
    start_perm = permutations::permutation("(1,2,3)(4,5)"), show_progress_bar = TRUE
  ))

  g1 <- gips(matrix_invariant_by_example_perm, 13, estimated_mean = FALSE)

  expect_error(find_MAP(g1, 10, return_probabilities = TRUE, optimizer = "HC"))
  expect_error(find_MAP(g1, 10, optimizer = "BD"))
  expect_error(
    find_MAP(g1, 10,
      optimizer = c("MH", "BG")
    ),
    "`optimizer` must be the character vector of length 1."
  )
  expect_error(find_MAP(g1, 10, optimizer = "BD"))
  expect_error(
    find_MAP(g1, max_iter = NA, optimizer = "MH"),
    "You provided `optimizer == MH` and `max_iter = NA`."
  )
  expect_error(
    find_MAP(g1, Inf, optimizer = "MH"),
    "`max_iter` in `Metropolis_Hastings\\(\\)` must be finite."
  )
  expect_error(
    find_MAP(g1, 10, optimizer = "continue"),
    "You provided `optimizer == 'continue'`, but the gips object `g` is not optimized."
  )

  g_small <- gips(matrix(c(1, 0.5, 0.5, 5), ncol = 2), 13, estimated_mean = FALSE)
  g_BF <- find_MAP(g_small, optimizer = "full", show_progress_bar = FALSE)
  expect_error(
    find_MAP(g_BF, max_iter = 10, optimizer = "continue"),
    "`optimizer == 'continue'` cannot be provided after optimizating with `optimizer == 'brute_force'`"
  )

  matrix_big_space <- matrix(runif(35 * 35), nrow = 35)
  matrix_big_space <- t(matrix_big_space) %*% matrix_big_space
  g_big <- gips(matrix_big_space, number_of_observations = 13, estimated_mean = FALSE)

  expect_error(
    find_MAP(g_big, optimizer = "full", show_progress_bar = FALSE),
    "Optimizer 'brute_force' cannot browse such a big permutional space."
  )
})

test_that("Handle proper parameters", {
  g1 <- gips(matrix_invariant_by_example_perm, 13, estimated_mean = FALSE)

  expect_message(
    g_map <- find_MAP(g1, 10, show_progress_bar = FALSE),
    "The 'optimizer = NA' was automatically changed to 'optimizer = \"MH\"."
  )

  expect_silent(out <- Metropolis_Hastings(
    S = matrix_invariant_by_example_perm,
    number_of_observations = 13, max_iter = 2,
    show_progress_bar = FALSE
  ))
  expect_output(
    out <- Metropolis_Hastings(
      S = matrix_invariant_by_example_perm,
      number_of_observations = 13, max_iter = 3,
      start_perm = permutations::permutation("(1,2,3)(4,5)"),
      show_progress_bar = TRUE
    ),
    "===="
  )
  expect_output(
    out <- hill_climbing(
      S = matrix_invariant_by_example_perm,
      number_of_observations = 13, max_iter = 8,
      show_progress_bar = TRUE
    ),
    "===="
  )
  expect_warning(expect_output(
    out <- hill_climbing(
      S = matrix_invariant_by_example_perm,
      number_of_observations = 13, max_iter = 2,
      show_progress_bar = TRUE
    ),
    "===="
  )) # Can we stack the `expect`s like that? Looks neat!

  g_small <- gips(matrix(c(1, 0.5, 0.5, 5), ncol = 2), 13, estimated_mean = FALSE)
  expect_message(
    find_MAP(g_small, max_iter = 10, optimizer = "MH", show_progress_bar = FALSE),
    "Consider using `optimizer = 'brute_force'`, because it will use 2! \\(factorial\\) = 2 iterations and will browse all permutations, therefore it will definitely find the maximum posteriori estimator."
  )

  expect_warning(
    find_MAP(g1, max_iter = 2, optimizer = "HC", show_progress_bar = FALSE),
    "Hill Climbing algorithm did not converge in 2 iterations!"
  )
})

test_that("Warns when found group has n0 > n", {
  number_of_observations <- 2
  Z <- MASS::mvrnorm(number_of_observations,
    mu = numeric(6),
    Sigma = matrix_invariant_by_example_perm
  )
  S <- (t(Z) %*% Z) / number_of_observations

  g <- gips(S, number_of_observations, estimated_mean = FALSE)

  expect_warning(
    find_MAP(g,
      max_iter = 2, show_progress_bar = FALSE,
      optimizer = "MH"
    ), # after 2 iterations the n0 can either be 6, or 5, or 4. It cannot be 1 or 2.
    "n0"
  )
})

test_that("find_MAP() can gues the correct optimizer and message the user", {
  g <- gips(matrix_invariant_by_example_perm, 13, estimated_mean = FALSE)

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

test_that("find_MAP will remember the right number of observations and estimated_mean", {
  number_of_observations <- 13
  
  g_em_FALSE <- gips(matrix_invariant_by_example_perm, number_of_observations,
                     estimated_mean = FALSE)
  g_em_TRUE <- gips(matrix_invariant_by_example_perm, number_of_observations,
                    estimated_mean = TRUE)

  g_map_em_FALSE <- find_MAP(g_em_FALSE, max_iter = 10,
                             show_progress_bar = FALSE, optimizer = "MH")
  g_map_em_TRUE <- find_MAP(g_em_TRUE, max_iter = 10,
                            show_progress_bar = FALSE, optimizer = "MH")
  
  expect_equal(attr(g_map_em_FALSE, "number_of_observations"),
               number_of_observations)
  expect_equal(attr(g_map_em_TRUE, "number_of_observations"),
               number_of_observations)
})
