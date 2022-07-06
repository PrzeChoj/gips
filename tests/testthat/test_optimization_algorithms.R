test_examples("../..") # example for Metropolis_Hastings and best_growth are here

test_that("Infinite max_iter and progress bar is not supported on best_growth", {
  expect_error(best_growth(U=U, number_of_observations=number_of_observations, max_iter=Inf,
                           start_perm=start_perm, show_progress_bar=TRUE),
               "Progress bar is not yet supported for infinite max_iter. Rerun the algorithm with show_progress_bar=FALSE or finite max_iter. For more information see ISSUE#8.")
})

# TODO()
