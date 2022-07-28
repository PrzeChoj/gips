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

test_that('new_gips works or throws an erron on wrong arguments', {
  expect_silent(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                         S, number_of_observations, 3, diag(nrow = ncol(S)), NULL))
  
  expect_error(new_gips(gips_perm('(1,2)(3,4,5,6)', 6),
                        S, number_of_observations, 3,
                        diag(nrow = ncol(S)), NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        7, number_of_observations, 3,
                        diag(nrow = ncol(S)), NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        S, number_of_observations+0.1, 3,
                        diag(nrow = ncol(S)), NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        S, 'r', 3,
                        diag(nrow = ncol(S)), NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        S, number_of_observations, 'r',
                        diag(nrow = ncol(S)), NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        S, number_of_observations, 3,
                        7, NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        S, number_of_observations, 3,
                        'diag(nrow = ncol(S))', NULL))
  expect_error(new_gips(list(gips_perm('(1,2)(3,4,5,6)', 6)),
                        S, number_of_observations, 3,
                        diag(nrow = ncol(S)), 'NULL'))
})


test_that('Properly validate the gips class with no optimization or after a single optimization',{
  custom_perm1 <- gips_perm('(1,2)(3,4,5,6)', 6)
  g1 <- gips(S, number_of_observations, perm = custom_perm1)
  
  g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  while (attr(g2, "optimization_info")$acceptance_rate == 0) {  # Around 4% of time, in the optimization all permutations were rejected. Is such a case, try again.
    g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  }
  g3 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="BG", return_probabilities = FALSE)
  
  
  # tests:
  g_err <- g2
  class(g_err[[1]]) <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  g_err[[2]] <- g_err[[1]]
  expect_error(validate_gips(g_err))
  
  g_err <- '(1,2)(3,4,5)'
  attributes(g_err) <- attributes(g2)
  expect_error(validate_gips(g_err))
  
  
  # test of "optimization_info" validation:
  g_err <- g2
  attr(g_err, "optimization_info") <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- NULL
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- NULL
  expect_error(validate_gips(g_err))  # this one shows an error that one have the list of 10 elements, which is actually expected, but the names of the fields are not expected.
  
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
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- 7
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- gips_perm('', 6)
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- 7
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- "7"
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
  expect_error(validate_gips(g_err))
  
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
  
  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- time_now
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- 1
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- (time_now - time_now)
  expect_error(validate_gips(g_err))
  
  
  
  # more than 5 problems at the time:
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  attr(g_err, "optimization_info")[["log_likelihood_values"]] <- as.character(attr(g_err, "optimization_info")[["log_likelihood_values"]])
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- 7
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 40
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1/(len_post_prob-1), len_post_prob-1))
  attr(g_err, "optimization_info")[["best_perm_log_likelihood"]] <- 7
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(validate_gips(g_err))
})

test_that('Properly validate the gips class after multiple optimizations',{
  custom_perm1 <- gips_perm('(1,2)(3,4,5,6)', 6)
  g1 <- gips(S, number_of_observations, perm = custom_perm1)
  
  g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  g2 <- find_gips(g2, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  g2 <- find_gips(g2, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  while (attr(g2, "optimization_info")$acceptance_rate == 0) {  # Around 4% of time, in the optimization all permutations were rejected. Is such a case, try again.
    g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
    g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
    g2 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="MH", return_probabilities = TRUE)
  }
  g3 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="BG", return_probabilities = FALSE)
  g3 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="BG", return_probabilities = FALSE)
  g3 <- find_gips(g1, max_iter=3, show_progress_bar=FALSE, optimizer="BG", return_probabilities = FALSE)
  
  
  
  
  # tests:
  g_err <- g2
  class(g_err[[1]]) <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  g_err[[2]] <- g_err[[1]]
  expect_error(validate_gips(g_err))
  
  g_err <- '(1,2)(3,4,5)'
  attributes(g_err) <- attributes(g2)
  expect_error(validate_gips(g_err))
  
  
  # test of "optimization_info" validation:
  g_err <- g2
  attr(g_err, "optimization_info") <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- NULL
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- NULL
  expect_error(validate_gips(g_err))  # this one shows an error that one have the list of 10 elements, which is actually expected, but the names of the fields are not expected.
  
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
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- 7
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- gips_perm('', 6)
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- 7
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- "7"
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 40
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
  expect_error(validate_gips(g_err))
  
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
  
  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- time_now
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- 1
  expect_error(validate_gips(g_err))
  
  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- (time_now - time_now)
  expect_error(validate_gips(g_err))
  
  
  
  # more than 5 problems at the time:
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  attr(g_err, "optimization_info")[["log_likelihood_values"]] <- as.character(attr(g_err, "optimization_info")[["log_likelihood_values"]])
  attr(g_err, "optimization_info")[["last_perm_log_likelihood"]] <- 7
  attr(g_err, "optimization_info")[["iterations_performed"]] <- c(2,2)
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1/(len_post_prob-1), len_post_prob-1))
  attr(g_err, "optimization_info")[["best_perm_log_likelihood"]] <- 7
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(validate_gips(g_err))
})

test_that("Process proper parameters", {
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                               start_perm=example_perm, delta=3, D_matrix=NULL,
                                               return_probabilities=FALSE, show_progress_bar=FALSE))
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations, max_iter=10,
                                               start_perm=gips_perm(example_perm, 6), delta=3, D_matrix=NULL,
                                               return_probabilities=FALSE, show_progress_bar=FALSE))
})

test_that('check_correctness_of_arguments properly validates arguments',{
  expect_silent(check_correctness_of_arguments(S, number_of_observations, 30,
                                               permutations::permutation('(1,3)(2,4)(5,6)'),
                                               3, diag(nrow = ncol(S)), FALSE, FALSE))
  
  expect_error(check_correctness_of_arguments(6, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(matrix(1:30, ncol=5),
                                              number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(matrix(c(LETTERS, LETTERS)[1:36], ncol=6),
                                              number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  S_nonsymetric <- S
  S_nonsymetric[1,2] <- S[1,2]-1
  expect_error(check_correctness_of_arguments(S_nonsymetric,
                                              number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  S_non_positive_semi_definite <- S - diag(ncol(S))*eigen(S, symmetric = TRUE, only.values = TRUE)$values[2]
  expect_error(check_correctness_of_arguments(S_non_positive_semi_definite,
                                              number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, NULL, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, 0, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations+0.1, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30.1,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 1,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              '(1,3)(2,4)(5,6)',
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              gips_perm('(1,3)(2,4)(5,6)', 7),
                                              3, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              NULL, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              1.9, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              2, diag(nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)+1), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, 7, FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, matrix(1:30, nrow = ncol(S)), FALSE, FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), 'FALSE', FALSE))
  expect_error(check_correctness_of_arguments(S, number_of_observations, 30,
                                              permutations::permutation('(1,3)(2,4)(5,6)'),
                                              3, diag(nrow = ncol(S)), FALSE, 'FALSE'))
  
  expect_error(check_correctness_of_arguments(S, number_of_observations+0.1, 1,
                                              '(1,3)(2,4)(5,6)',
                                              2, diag(nrow = ncol(S)), 'FALSE', 'FALSE'))
  
  # old tests:
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


test_that('print.gips() works', {
  g <- gips(S, number_of_observations)
  expect_output(print(g), "The permutation \\(\\) has log likelihood")
  
  g <- find_gips(g, 10, show_progress_bar = FALSE, optimizer = "MH")
  expect_output(print(g), "which was found after 10 log_likelihood calculations.")
})

test_that('plot.gips works or abords for wrong arguments', {
  expect_error(plot.gips(custom_perm1))
  
  expect_error(plot(g1, type="both"))
  expect_message(plot(g1))
  
  g1_found <- find_gips(g1, 3, show_progress_bar = FALSE, optimizer = "MH")
  expect_message(plot(g1_found))
  expect_silent(plot(g1_found, type = "both"))
  expect_silent(plot(g1_found, type = "all", logarithmic_y = FALSE))
  expect_silent(plot(g1_found, type = "best", logarithmic_x = TRUE,
                     ylim = range(attr(g1_found, "optimization_info")["log_likelihood_values"])*2))
  
  expect_error(plot.gips(g1_found, type = "non_existing"))
})



