# examples already tested

perm_size <- 6
mu <- numeric(perm_size)
# sigma is a matrix invariant under permutation (1,2,3,4,5,6)
sigma_matrix <- matrix(data = c(1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
                                0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
                                0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
                                0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
                                0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
                                0.8, 0.6, 0.4, 0.6, 0.8, 1.0),
                       nrow=perm_size, byrow = TRUE)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- (t(Z) %*% Z) / number_of_observations  # the theoretical mean is 0

custom_perm1 <- gips_perm('(1,2)(3,4,5,6)', 6)
g1 <- gips(S, number_of_observations, perm = custom_perm1)

custom_perm2 <- permutations::permutation('(1,2)(3,4,5)(6)')  # notice, the `custom_perm2` does not remember the 6
g2 <- gips(S, number_of_observations, perm = custom_perm2)

test_that('Setting custom permutation in gips constructor works',{
  expect_identical(custom_perm1, g1[[1]])
  
  expect_identical(gips_perm(custom_perm2, ncol(S)), g2[[1]])
})


test_that('Properly validate the gips class',{
  custom_perm1 <- gips_perm('(1,2)(3,4,5,6)', 6)
  g1 <- gips(S, number_of_observations, perm = custom_perm1)
  
  g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  g3 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="BG", return_probabilities = FALSE)
  
  
  # tests:
  g_err <- g2
  g_err[[2]] <- g_err[[1]]
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  class(g_err[[1]]) <- "test"
  expect_error(validate_gips(g_err))
  
  
  # test of "optimization_info" validation:
  g_err <- g2
  attr(g_err, "optimization_info") <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["log_likelihood_values"]] <- as.character(attr(g_err, "optimization_info")[["log_likelihood_values"]])
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["visited_perms"]] <- NA
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  class(attr(g_err, "optimization_info")[["visited_perms"]][[1]]) <- "test"
  expect_error(validate_gips(g_err))
  
  # TODO(Test for compering the last_perm with the last element of visited_perms)
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- 7
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- "7"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- c(2,2)
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 4  # there were 2 iterations + init
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "MH"  # Even if MH was used, it would produce the text "Metropolis_Hastings"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "best_growth"  # best_growth is legal, but the post_probabilities are not with this optimization_algorithm_used
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + 0.1
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  len_post_prob <- length(attr(g_err, "optimization_info")[["post_probabilities"]])
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1/(len_post_prob-1), len_post_prob-1))
  expect_error(validate_gips(g_err))  # TODO(investigate --->>>  once this `validate` did not throw the expected error!)
  
  g_err <- g2
  attr(g_err, "optimization_info")[["did_converge"]] <- TRUE
  expect_error(validate_gips(g_err))
  
  g_err <- g3
  attr(g_err, "optimization_info")[["did_converge"]] <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g3
  attr(g_err, "optimization_info")[["did_converge"]] <- NA
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["best_perm_log_likelihood"]] <- 7
  expect_error(validate_gips(g_err))
  
  # more than 5 problems at the time:
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  attr(g_err, "optimization_info")[["log_likelihood_values"]] <- as.character(attr(g_err, "optimization_info")[["log_likelihood_values"]])
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- 7
  attr(g_err, "optimization_info")[["iterations_performed"]] <- c(2,2)
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1/(len_post_prob-1), len_post_prob-1))
  attr(g_err, "optimization_info")[["best_perm_log_likelihood"]] <- 7
  expect_error(validate_gips(g_err))
})


test_that('Plot does not work for non_optimized gips', {
  expect_error(plot.gips(custom_perm1))
  
  expect_error(plot(g1))
  g1_found <- find_gips(g1, 3, show_progress_bar = FALSE)
  expect_silent(plot(g1_found))
})



