# examples already tested

# Setup the S matrix for every test
perm_size <- 6
mu <- numeric(perm_size)
# sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
sigma_matrix <- matrix(
  data = c(
    1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.0
  ),
  nrow = perm_size, byrow = TRUE
)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix) # this place is random
S <- (t(Z) %*% Z) / number_of_observations # the theoretical mean is 0




test_that("Setting custom permutation in gips constructor works", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
             was_mean_estimated = FALSE, perm = custom_perm1)

  custom_perm2 <- permutations::permutation("(1,2)(3,4,5)(6)")
  g2 <- gips(S, number_of_observations,
             was_mean_estimated = FALSE, perm = custom_perm2)

  expect_identical(custom_perm1, g1[[1]])

  expect_identical(gips_perm(custom_perm2, ncol(S)), g2[[1]])
})

test_that("new_gips works or throws an erron on wrong arguments", {
  expect_silent(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations, 3, diag(nrow = ncol(S)), FALSE, NULL
  ))

  expect_error(new_gips(
    gips_perm("(1,2)(3,4,5,6)", 6),
    S, number_of_observations, 3,
    diag(nrow = ncol(S)), FALSE, NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    7, number_of_observations, 3,
    diag(nrow = ncol(S)), FALSE, NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations + 0.1, 3,
    diag(nrow = ncol(S)), FALSE, NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, "r", 3,
    diag(nrow = ncol(S)), FALSE, NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations, "r",
    diag(nrow = ncol(S)), NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations, 3,
    7, FALSE, NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations, 3,
    "diag(nrow = ncol(S))", FALSE, NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations, 3, diag(nrow = ncol(S)), "FALSE", NULL
  ))
  expect_error(new_gips(
    list(gips_perm("(1,2)(3,4,5,6)", 6)),
    S, number_of_observations, 3,
    diag(nrow = ncol(S)), FALSE, "NULL"
  ))
})


test_that("Properly validate the gips class with no optimization or after a single optimization", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations, was_mean_estimated = FALSE,
             perm = custom_perm1)

  g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE,
                 optimizer = "MH", return_probabilities = TRUE)
  while (attr(g2, "optimization_info")[["acceptance_rate"]] == 0) { # Around 4% of time, in the optimization all permutations were rejected. Is such a case, try again. We want g2 to have at least 1 accepted transposition.
    g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE,
                   optimizer = "MH", return_probabilities = TRUE)
  }
  g3 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE,
                 optimizer = "HC", return_probabilities = FALSE)

  g_BF <- find_MAP(gips(
    matrix(c(1, 0.5, 0.5, 5), nrow = 2),
    number_of_observations, was_mean_estimated = FALSE
  ),
  optimizer = "BF", show_progress_bar = FALSE
  )


  # tests:
  g_err <- g2
  class(g_err[[1]]) <- "test"
  expect_error(validate_gips(g_err))

  g_err <- g2
  g_err[[2]] <- g_err[[1]]
  expect_error(validate_gips(g_err))

  g_err <- "(1,2)(3,4,5)"
  attributes(g_err) <- attributes(g2)
  expect_error(validate_gips(g_err))

  g_err <- list("(1,2)(3,4,5)")
  class(g_err[[1]]) <- "gips_perm"
  attributes(g_err) <- attributes(g2)
  expect_error(
    validate_gips(g_err),
    "You provided `g\\[\\[1\\]\\]` with `class\\(g\\[\\[1\\]\\]\\) == 'gips_perm'`, but your g\\[\\[1\\]\\] does not pass `validate_gips_perm\\(g\\[\\[1\\]\\]\\)`."
  )


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
  expect_error(validate_gips(g_err)) # this one shows an error that one have the list of 10 elements, which is actually expected, but the names of the fields are not expected.

  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "brute_force" # in general, brute_force is all right, but then the acceptance_rate has to be NULL
  expect_error(
    validate_gips(g_err),
    "When brute force algorithm was used for optimization, `attr\\(g, 'optimization_info'\\)\\[\\['acceptance_rate'\\]\\]` must be a NULL"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["log_posteriori_values"]] <- as.character(attr(g_err, "optimization_info")[["log_posteriori_values"]])
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["visited_perms"]] <- NA
  expect_error(validate_gips(g_err))

  g_err <- g2
  class(attr(g_err, "optimization_info")[["visited_perms"]][[1]]) <- "test"
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["visited_perms"]] <- list()
  expect_error(validate_gips(g_err), "is a list, but of a length")

  g_err <- g_BF
  attr(g_err, "optimization_info")[["visited_perms"]] <- list()
  expect_error(validate_gips(g_err), "is a list, but of a length")

  g_err <- g_BF
  attr(g_err, "optimization_info")[["visited_perms"]] <- "()"
  expect_error(
    validate_gips(g_err),
    "You have `attr\\(g, 'optimization_info'\\)\\[\\['visited_perms'\\]\\]` of type character."
  )

  g_err <- g_BF
  attr(g_err, "optimization_info")[["last_perm"]] <- "()"
  expect_error(
    validate_gips(g_err),
    "You have `attr\\(g, 'optimization_info'\\)\\[\\['last_perm'\\]\\]` of type character."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- 7
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- gips_perm("", 6)
  expect_error(validate_gips(g_err))

  g_err <- g2
  fake_gips_perm <- "()"
  class(fake_gips_perm) <- "gips_perm"
  attr(g_err, "optimization_info")[["last_perm"]] <- fake_gips_perm
  expect_error(
    validate_gips(g_err),
    ", but your attr\\(g, 'optimization_info'\\)\\[\\['last_perm'\\]\\] does not pass `validate_gips_perm"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- 7
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- "7"
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 4 # there were 2 iterations + init (+ additional 2 iterations if repeated); so it cannot be 4
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 4.5
  expect_error(
    validate_gips(g_err),
    "attr\\(g, 'optimization_info'\\)\\[\\['iterations_performed'\\]\\]` must be a vector of whole numbers"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "MH" # Even if MH was used, it would produce the text "Metropolis_Hastings"
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "hill_climbing" # hill_climbing is legal, but the post_probabilities are not with this optimization_algorithm_used
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + 0.1
  expect_error(validate_gips(g_err))

  g_err <- g2
  len_post_prob <- length(attr(g_err, "optimization_info")[["post_probabilities"]])
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1 / (len_post_prob - 1), len_post_prob - 1))
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
  attr(g_err, "optimization_info")[["best_perm_log_posteriori"]] <- 7
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

  g_err <- g2
  attr(g_err, "optimization_info")[["full_optimization_time"]] <- NA
  expect_error(
    validate_gips(g_err),
    "You have `is.na\\(attr\\(g, 'optimization_info'\\)\\[\\['full_optimization_time'\\]\\]\\)"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["full_optimization_time"]] <- 7
  expect_error(
    validate_gips(g_err),
    "You have `attr\\(g, 'optimization_info'\\)\\[\\['full_optimization_time'\\]\\]` of a class \\(numeric\\)"
  )

  time_now1 <- Sys.time()
  g_err <- g2
  time_now2 <- Sys.time()
  attr(g_err, "optimization_info")[["full_optimization_time"]] <- time_now1 - time_now2
  expect_error(
    validate_gips(g_err),
    "`attr\\(g, 'optimization_info'\\)\\[\\['full_optimization_time'\\]\\]` has to be a time difference bigger than 0."
  )


  # more than 5 problems at the time:
  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  attr(g_err, "optimization_info")[["log_posteriori_values"]] <- as.character(attr(g_err, "optimization_info")[["log_posteriori_values"]])
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- 7
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 40
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1 / (len_post_prob - 1), len_post_prob - 1))
  attr(g_err, "optimization_info")[["best_perm_log_posteriori"]] <- 7
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(validate_gips(g_err))
})

test_that("Properly validate the gips class after multiple optimizations", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
             was_mean_estimated = FALSE, perm = custom_perm1)

  g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE)
  g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE)
  g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE)
  while (attr(g2, "optimization_info")[["acceptance_rate"]] == 0) { # Around 4% of time, in the optimization all permutations were rejected. Is such a case, try again.
    g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE)
    g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE)
    g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE)
  }
  g3 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "HC", return_probabilities = FALSE)
  g3 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "HC", return_probabilities = FALSE)
  g3 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "HC", return_probabilities = FALSE)




  # tests:
  g_err <- g2
  class(g_err[[1]]) <- "test"
  expect_error(validate_gips(g_err))

  g_err <- g2
  g_err[[2]] <- g_err[[1]]
  expect_error(validate_gips(g_err))

  g_err <- "(1,2)(3,4,5)"
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
  expect_error(validate_gips(g_err)) # this one shows an error that one have the list of 10 elements, which is actually expected, but the names of the fields are not expected.

  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["log_posteriori_values"]] <- as.character(attr(g_err, "optimization_info")[["log_posteriori_values"]])
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
  attr(g_err, "optimization_info")[["last_perm"]] <- gips_perm("", 6)
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- 7
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- "7"
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 40
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "MH" # Even if MH was used, it would produce the text "Metropolis_Hastings"
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "hill_climbing" # hill_climbing is legal, but the post_probabilities are not with this optimization_algorithm_used
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + 0.1
  expect_error(validate_gips(g_err))

  g_err <- g2
  len_post_prob <- length(attr(g_err, "optimization_info")[["post_probabilities"]])
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1 / (len_post_prob - 1), len_post_prob - 1))
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
  attr(g_err, "optimization_info")[["best_perm_log_posteriori"]] <- 7
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
  attr(g_err, "optimization_info")[["log_posteriori_values"]] <- as.character(attr(g_err, "optimization_info")[["log_posteriori_values"]])
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- 7
  attr(g_err, "optimization_info")[["iterations_performed"]] <- c(2, 2)
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1 / (len_post_prob - 1), len_post_prob - 1))
  attr(g_err, "optimization_info")[["best_perm_log_posteriori"]] <- 7
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(validate_gips(g_err))
})

test_that("Process proper parameters", {
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = number_of_observations, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = number_of_observations, max_iter = 10,
    start_perm = gips_perm(example_perm, 6), delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
})

test_that("check_correctness_of_arguments properly validates arguments", {
  expect_silent(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))

  expect_error(check_correctness_of_arguments(
    6, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    matrix(1:30, ncol = 5),
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    matrix(c(LETTERS, LETTERS)[1:36], ncol = 6),
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))

  S_nonsymetric <- S
  S_nonsymetric[1, 2] <- S[1, 2] - 1
  expect_error(check_correctness_of_arguments(
    S_nonsymetric,
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))

  S_non_positive_semi_definite <- S - diag(ncol(S)) * eigen(S, symmetric = TRUE, only.values = TRUE)[["values"]][2]
  expect_error(check_correctness_of_arguments(
    S_non_positive_semi_definite,
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, NULL, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, 0, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations + 0.1, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30.1,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 1,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    "(1,3)(2,4)(5,6)",
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    gips_perm("(1,3)(2,4)(5,6)", 7),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    NULL, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    1.9, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    2, diag(nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S) + 1), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, 7, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, matrix(1:30, nrow = ncol(S)), FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), "FALSE", FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, "FALSE", FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations + 0.1, 1,
    "(1,3)(2,4)(5,6)",
    2, diag(nrow = ncol(S)), FALSE, FALSE, "FALSE"
  ))
  
  # A lot of problems ot once
  expect_error(check_correctness_of_arguments(
    S, number_of_observations + 0.1, 1,
    "(1,3)(2,4)(5,6)",
    2, diag(nrow = ncol(S)), "FALSE", "FALSE", "FALSE"
  ), "7 problems identified with provided arguments")

  # old tests:
  # A single problem at the same time:
  expect_error(check_correctness_of_arguments(
    S = NULL, number_of_observations, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S = matrix(c(1:30), nrow = 6),
    number_of_observations, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = NULL, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = 0, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = 15.5, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10.5,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 0,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 1, # TODO(Make it work for max_iter = 1)
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = NULL, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = NULL, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 2, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = 7,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3,
    D_matrix = matrix(c(1:30), nrow = 6),
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3,
    D_matrix = NULL,
    was_mean_estimated = "FALSE",
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3,
    D_matrix = NULL,
    was_mean_estimated = NA,
    return_probabilities = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = 7, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = NULL, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = 7
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE,
    return_probabilities = FALSE, show_progress_bar = NULL
  ))


  # A number of problems at the same time. Not all are printed:
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = -1, max_iter = 1,
    start_perm = example_perm, delta = 2, D_matrix = 7,
    was_mean_estimated = NA,
    return_probabilities = NULL, show_progress_bar = NULL
  ), "\\.\\.\\. and 2 more problems")
})


test_that("print.gips() works", {
  g <- gips(S, number_of_observations, was_mean_estimated = FALSE)
  expect_output(print(g), "The permutation \\(\\) is 1 times more likely")

  g <- find_MAP(g, 10, show_progress_bar = FALSE, optimizer = "MH")
  expect_output(print(g), "was found after 10 log_posteriori calculations, is")
})

test_that("plot.gips() works or abords for wrong arguments", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
             was_mean_estimated = FALSE, perm = custom_perm1)

  expect_error(plot.gips(custom_perm1))

  expect_error(plot(g1, type = "both"))
  expect_message(
    plot(g1),
    "`type = NA` was automatically changed to `type = 'heatmap'`"
  )

  g1_found <- find_MAP(g1, 3, show_progress_bar = FALSE, optimizer = "MH")
  expect_message(
    plot(g1_found),
    "`type = NA` was automatically changed to `type = 'both'`"
  )
  expect_silent(plot(g1_found, type = "both"))
  expect_silent(plot(g1_found, type = "all", logarithmic_y = FALSE))
  expect_silent(plot(g1_found, type = "best"))
  expect_silent(plot(g1_found,
    type = "best", logarithmic_x = TRUE,
    ylim = range(attr(g1_found, "optimization_info")["log_posteriori_values"]) + c(-10, 10)
  ))

  expect_error(
    plot.gips(g1_found, type = "non_existing"),
    "Did You misspell the 'type' argument?"
  )
})

test_that("summary.gips() works", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
             was_mean_estimated = FALSE, perm = custom_perm1)
  
  start_permutation_log_posteriori <- log_posteriori_of_gips(g1)
  log_posteriori_id <- log_posteriori_of_perm(
    perm_proposal = "", S = S,
    number_of_observations = number_of_observations,
    delta = attr(g1, "delta"), D_matrix = attr(g1, "D_matrix")
  )

  my_sum <- structure(list(
    optimized = FALSE, start_permutation = structure(list(
      c(1, 2), c(3, 4, 5, 6)
    ), size = 6, class = "gips_perm"),
    start_permutation_log_posteriori = start_permutation_log_posteriori,
    times_more_likely_than_id = exp(start_permutation_log_posteriori - log_posteriori_id),
    n0 = 2, S_matrix = S, number_of_observations = 13,
    was_mean_estimated = FALSE,
    delta = 3, D_matrix = structure(c(
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1
    ), .Dim = c(6L, 6L))
  ), class = "summary.gips")

  expect_identical(summary(g1), my_sum)
  
  expect_output(
    print(summary(g1)),
    "Number of observations is bigger than n0 for this permutaion,\nso "
  )

  g2 <- find_MAP(g1, max_iter = 10, optimizer = "MH", show_progress_bar = FALSE)

  expect_output(
    print(summary(g2)),
    "Optimization time:"
  )

  g3 <- find_MAP(g2,
    max_iter = 10, optimizer = "continue",
    show_progress_bar = FALSE
  )

  expect_output(
    print(summary(g3)),
    "Optimization algorithms:"
  )
})

test_that("start_permutation_log_posteriori was calculated correctly", {
  g <- gips(S, number_of_observations, was_mean_estimated = TRUE)
  g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE,
                    optimizer = "MH")
  
  optimization_info <- attr(g_map, "optimization_info")
  expect_equal(optimization_info[["log_posteriori_values"]][1],
               log_posteriori_of_perm(
                 "",
                 attr(g_map, "S"),
                 attr(g_map, "number_of_observations") - attr(g_map, "was_mean_estimated"),
                 attr(g_map, "delta"),
                 attr(g_map, "D_matrix")
               ))
})
