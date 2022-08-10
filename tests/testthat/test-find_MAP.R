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
