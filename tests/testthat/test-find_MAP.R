test_examples("../..") # example for Metropolis_Hastings and hill_climbing are here

test_that("Handle improper parameters", {
  expect_error(hill_climbing(
    S = matrix_invariant_by_example_perm, number_of_observations = 13,
    max_iter = Inf, show_progress_bar = TRUE
  ))
  expect_error(find_MAP(gips(matrix_invariant_by_example_perm, number_of_observations), 10, return_probabilities = TRUE, optimizer = "HC"))
  expect_error(find_MAP(gips(matrix_invariant_by_example_perm, number_of_observations), 10, optimizer = "BD"))
  expect_error(Metropolis_Hastings(
    S = matrix_invariant_by_example_perm, number_of_observations = 13, max_iter = Inf,
    start_perm = permutations::permutation("(1,2,3)(4,5)"), show_progress_bar = TRUE
  ))
})

test_that("Handle proper parameters", {
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
})

test_that("Warns when found group has n0 > n", {
  number_of_observations <- 2
  Z <- MASS::mvrnorm(number_of_observations,
    mu = numeric(6),
    Sigma = matrix_invariant_by_example_perm
  )
  S <- (t(Z) %*% Z) / number_of_observations

  g <- gips(S, number_of_observations)

  expect_warning(
    find_MAP(g,
      max_iter = 2, show_progress_bar = FALSE,
      optimizer = "MH"
    ), # after 2 iterations the n0 can either be 6, or 5, or 4. It cannot be 1 or 2.
    "n0"
  )
})

test_that("find_MAP() can gues the correct optimizer and message the user", {
  g <- gips(matrix_invariant_by_example_perm, 13)

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
