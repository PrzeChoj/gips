test_examples("../..") # example for Metropolis_Hastings and best_growth are here

test_that("Process proper parameters", {
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                             start_perm=example_perm, delta=3, D_matrix=NULL,
                                             return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                             start_perm=gips_perm(example_perm, 6), delta=3, D_matrix=NULL,
                                             return_probabilities=FALSE, show_progress_bar=FALSE))
})

test_that("Handle improper parameters", {
  expect_error(best_growth(S=matrix_invariant_by_example_perm, number_of_observations=number_of_observations, max_iter=Inf,
                           start_perm=example_perm, show_progress_bar=TRUE))
  expect_error(find_gips(gips(matrix_invariant_by_example_perm/5, number_of_observations), 10, return_probabilities = TRUE, optimizer = "BG"))  # TODO(Do not divide this by 5)
  
  # A single problem at the same time:
  expect_error(check_correctness_of_arguments(S=NULL, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(S=matrix(c(1:30), nrow = 6),
                                            number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations=NULL, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations=0, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations=15.5, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10.5,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=0,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=1,  # TODO(Make it work for max_iter = 1)
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=NULL, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=NULL, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=2, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=7,
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3,
                                            D_matrix=matrix(c(1:30), nrow = 6),
                                            return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=7, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=NULL, show_progress_bar=FALSE))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=7))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                            start_perm=example_perm, delta=3, D_matrix=NULL,
                                            return_probabilities=FALSE, show_progress_bar=NULL))
  
  
  # A number of problems at the same time. Not all are printed:
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations=-1, max_iter=1,
                                            start_perm=example_perm, delta=2, D_matrix=7,
                                            return_probabilities=NULL, show_progress_bar=NULL))
})
