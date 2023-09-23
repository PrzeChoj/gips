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
  g1 <- gips(
    S, number_of_observations,
    was_mean_estimated = FALSE, perm = custom_perm1
  )

  custom_perm2 <- permutations::permutation("(1,2)(3,4,5)(6)")
  g2 <- gips(
    S, number_of_observations,
    was_mean_estimated = FALSE, perm = custom_perm2
  )

  expect_identical(custom_perm1, g1[[1]])

  expect_identical(gips_perm(custom_perm2, ncol(S)), g2[[1]])

  # gips as the `perm` parameter
  g3 <- gips(
    S, number_of_observations,
    was_mean_estimated = FALSE, perm = g2
  )

  expect_equal(g2, g3)
})

test_that("new_gips() works or throws an erron on wrong arguments", {
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
  g1 <- gips(
    S, number_of_observations,
    was_mean_estimated = FALSE,
    perm = custom_perm1
  )

  g2 <- find_MAP(g1,
    max_iter = 3, show_progress_bar = FALSE,
    optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE
  )
  while (attr(g2, "optimization_info")[["acceptance_rate"]] == 0) { # Around 4% of time, in the optimization all permutations were rejected. Is such a case, try again. We want g2 to have at least 1 accepted transposition.
    g2 <- find_MAP(g1,
      max_iter = 3, show_progress_bar = FALSE,
      optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE
    )
  }
  g3 <- find_MAP(g1,
    max_iter = 3, show_progress_bar = FALSE,
    optimizer = "HC", return_probabilities = FALSE
  )

  g_BF <- find_MAP(
    gips(
      matrix(c(
        1, 0.5, 0.5,
        0.5, 1, 0.5,
        0.5, 0.5, 1
      ), nrow = 3),
      number_of_observations,
      was_mean_estimated = FALSE
    ),
    optimizer = "BF", show_progress_bar = FALSE
  )

  g_BF_prob <- find_MAP(
    gips(
      matrix(c(
        1, 0.5, 0.5,
        0.5, 1, 0.5,
        0.5, 0.5, 1
      ), nrow = 3),
      number_of_observations,
      was_mean_estimated = FALSE
    ),
    optimizer = "BF", show_progress_bar = FALSE,
    return_probabilities = TRUE, save_all_perms = TRUE
  )


  # tests:
  expect_silent(validate_gips(g1))
  expect_silent(validate_gips(g2))
  expect_silent(validate_gips(g3))
  expect_silent(validate_gips(g_BF))
  expect_silent(validate_gips(g_BF_prob))



  # NaNs or Infs in D_matrix
  expect_error(gips(
    S[1:4, 1:4], number_of_observations,
    was_mean_estimated = FALSE, D_matrix = diag(4) * 1e350
  ))
  expect_error(gips(
    S[1:4, 1:4], number_of_observations,
    was_mean_estimated = FALSE, D_matrix = diag(Inf, 4)
  ))

  # Other tests
  g_err <- g2
  class(g_err[[1]]) <- "test"
  expect_error(
    validate_gips(g_err),
    "must be an object of a `gips_perm` class."
  )

  g_err <- g2
  g_err[[2]] <- g_err[[1]]
  expect_error(
    validate_gips(g_err),
    "The `length\\(g\\)` must be `1`."
  )

  g_err <- "(1,2)(3,4,5)"
  attributes(g_err) <- attributes(g2)
  expect_error(
    validate_gips(g_err),
    "The `g` must be a list."
  )

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
  expect_error(
    validate_gips(g_err),
    "The `optimization_info` value must be either a `NULL`, or a list."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  expect_error(
    validate_gips(g_err),
    "You have a list of 15 elements."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- NULL
  expect_error(
    validate_gips(g_err),
    "You have a list of 13 elements."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["non_existing"]] <- "test"
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- NULL
  expect_error(
    validate_gips(g_err),
    "You have a list of 14 elements."
  )
  # this one showed an error that one have the list of 13 elements, which is actually expected, but the names of the fields are not expected.

  g_err <- g2
  attr(g_err, "optimization_info")[["acceptance_rate"]] <- -0.1
  expect_error(
    validate_gips(g_err),
    "must be a number in range \\[0, 1\\]."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "brute_force" # in general, brute_force is all right, but then the acceptance_rate has to be NULL
  expect_error(
    validate_gips(g_err),
    "When brute force algorithm was used for optimization, `attr\\(g, 'optimization_info'\\)\\[\\['acceptance_rate'\\]\\]` must be a `NULL`"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["log_posteriori_values"]] <- as.character(attr(g_err, "optimization_info")[["log_posteriori_values"]])
  expect_error(
    validate_gips(g_err),
    "must be a vector of numbers."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["visited_perms"]] <- "text"
  expect_error(
    validate_gips(g_err),
    "must be a list or `NA`."
  )

  g_err <- g2
  class(attr(g_err, "optimization_info")[["visited_perms"]][[1]]) <- "test"
  expect_error(
    validate_gips(g_err),
    "must be of a `gips_perm` class."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["visited_perms"]] <- list()
  expect_error(
    validate_gips(g_err),
    "is a list, but of a length 0."
  )

  g_err <- g_BF
  attr(g_err, "optimization_info")[["visited_perms"]] <- list()
  expect_error(
    validate_gips(g_err),
    "is a list, but of a length 0"
  )

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
  expect_error(
    validate_gips(g_err),
    "must be the last element of"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm"]] <- gips_perm("", 6)
  expect_error(
    validate_gips(g_err),
    "must be the log_posteriori of"
  )

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
  expect_error(
    validate_gips(g_err),
    "must be the log_posteriori of"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["last_perm_log_posteriori"]] <- "7"
  expect_error(
    validate_gips(g_err),
    "must be an object of a `gips_perm` class."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 4 # there were 2 iterations + init (+ additional 2 iterations if repeated); so it cannot be 4
  expect_error(
    validate_gips(g_err),
    "In every iteration at least one value of log_posteriori is calculated."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["iterations_performed"]] <- 4.5
  expect_error(
    validate_gips(g_err),
    "attr\\(g, 'optimization_info'\\)\\[\\['iterations_performed'\\]\\]` must be a vector of whole numbers"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "MH" # Even if MH was used, it would produce the text "Metropolis_Hastings"
  expect_error(
    validate_gips(g_err),
    "The available optimization algorithms are 'Metropolis_Hastings', 'hill_climbing' and 'brute_force'."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_algorithm_used"]] <- "hill_climbing" # hill_climbing is legal, but the post_probabilities are not with this optimization_algorithm_used
  expect_error(
    validate_gips(g_err),
    "`post_probabilities` can only be obtained with 'Metropolis_Hastings' or 'brute_force' optimization method."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + 0.1
  expect_error(
    validate_gips(g_err),
    "must have properties of probability. All elements in range"
  )

  g_err <- g2
  len_post_prob <- length(attr(g_err, "optimization_info")[["post_probabilities"]])
  attr(g_err, "optimization_info")[["post_probabilities"]] <- attr(g_err, "optimization_info")[["post_probabilities"]] + c(0, rep(1 / (len_post_prob - 1), len_post_prob - 1))
  expect_error(
    validate_gips(g_err),
    "must have properties of probability. All elements in range"
  )

  g_good <- g2
  attr(g_good, "optimization_info")[["post_probabilities"]] <- c(1, rep(0, length(attr(g_good, "optimization_info")[["post_probabilities"]]) - 1)) # this could underflow to 0
  names(attr(g_good, "optimization_info")[["post_probabilities"]]) <- names(attr(g2, "optimization_info")[["post_probabilities"]])
  expect_silent(
    validate_gips(g_good)
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["did_converge"]] <- TRUE
  expect_error(
    validate_gips(g_err),
    "`did_converge` can only be obtained with 'hill_climbing' or 'brute_force' optimization method."
  )

  g_err <- g3
  attr(g_err, "optimization_info")[["did_converge"]] <- "test"
  expect_error(
    validate_gips(g_err),
    "When 'hill_climbing' optimization method, the `did_converge` must be `TRUE` or `FALSE`."
  )

  g_err <- g3
  attr(g_err, "optimization_info")[["did_converge"]] <- NA
  expect_error(
    validate_gips(g_err),
    "When 'hill_climbing' optimization method, the `did_converge` must be `TRUE` or `FALSE`."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["best_perm_log_posteriori"]] <- 7
  expect_error(validate_gips(g_err))

  g_err <- g2
  attr(g_err, "optimization_info")[["optimization_time"]] <- I(NA)
  expect_error(
    validate_gips(g_err),
    "is initially set to `NA`, but that state of the gips object should not be available to the user."
  )

  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- time_now
  expect_error(
    validate_gips(g_err),
    "has to be of a class 'difftime'."
  )

  g_err <- g2
  time_now <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- 1
  expect_error(
    validate_gips(g_err),
    " has to be of a class 'difftime'."
  )

  g_err <- g2
  time_now2 <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- (time_now - time_now2)
  expect_error(
    validate_gips(g_err),
    "has to be a non negative time difference."
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["whole_optimization_time"]] <- NA
  expect_error(
    validate_gips(g_err),
    "You have `is.na\\(attr\\(g, 'optimization_info'\\)\\[\\['whole_optimization_time'\\]\\]\\)"
  )

  g_err <- g2
  attr(g_err, "optimization_info")[["whole_optimization_time"]] <- 7
  expect_error(
    validate_gips(g_err),
    "You have `attr\\(g, 'optimization_info'\\)\\[\\['whole_optimization_time'\\]\\]` of a class \\(numeric\\)"
  )

  g_err <- g2
  time_now1 <- Sys.time()
  attr(g_err, "optimization_info")[["whole_optimization_time"]] <- time_now2 - time_now1
  expect_error(
    validate_gips(g_err),
    "`attr\\(g, 'optimization_info'\\)\\[\\['whole_optimization_time'\\]\\]` has to be a non negative time difference"
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
  expect_error(
    validate_gips(g_err),
    "and 3 more problems"
  )
})

test_that("Properly validate the gips class after multiple optimizations", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
    was_mean_estimated = FALSE, perm = custom_perm1
  )

  g2 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE)
  g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE)
  g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE)
  while (attr(g2, "optimization_info")[["acceptance_rate"]] == 0) { # Around 4% of time, in the optimization all permutations were rejected. Is such a case, try again.
    g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE)
    g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE)
    g2 <- find_MAP(g2, max_iter = 3, show_progress_bar = FALSE, optimizer = "MH", return_probabilities = TRUE, save_all_perms = TRUE)
  }
  g3 <- find_MAP(g1, max_iter = 3, show_progress_bar = FALSE, optimizer = "HC", return_probabilities = FALSE)
  g3 <- find_MAP(g3, max_iter = 3, show_progress_bar = FALSE, optimizer = "HC", return_probabilities = FALSE)
  g3 <- find_MAP(g3, max_iter = 3, show_progress_bar = FALSE, optimizer = "HC", return_probabilities = FALSE)

  expect_warning(expect_message(expect_message(
    g_MH_MH <- find_MAP(
      find_MAP(
        gips(
          matrix(c(1, 0.5, 0.5, 5), nrow = 2),
          number_of_observations,
          was_mean_estimated = FALSE
        ),
        optimizer = "MH", show_progress_bar = FALSE,
        return_probabilities = FALSE, max_iter = 3
      ),
      optimizer = "MH", show_progress_bar = FALSE,
      return_probabilities = TRUE, max_iter = 3,
      save_all_perms = TRUE
    )
  ))) # 2 messages and a warning




  # tests:
  expect_silent(validate_gips(g1))
  expect_silent(validate_gips(g2))
  expect_silent(validate_gips(g3))
  expect_silent(validate_gips(g_MH_MH))

  expect_true(length(attr(g2, "opt")[["visited_perms"]]) >= 9)
  expect_false(is.null(attr(g2, "optimization_info")[["post_probabilities"]]))
  expect_true(is.null(attr(g_MH_MH, "optimization_info")[["post_probabilities"]])) # first one was optimized without saving, so for the second one it was forgotten


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
  attr(g_err, "optimization_info")[["visited_perms"]] <- "text"
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
  time_now2 <- Sys.time()
  attr(g_err, "optimization_info")[["optimization_time"]] <- (time_now - time_now2)
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
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_silent(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = number_of_observations, max_iter = 10,
    start_perm = gips_perm(example_perm, 6), delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
})

test_that("check_correctness_of_arguments() properly validates arguments", {
  expect_silent(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))

  expect_error(check_correctness_of_arguments(
    6, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    matrix(1:30, ncol = 5),
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    matrix(c(LETTERS, LETTERS)[1:36], ncol = 6),
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))

  S_nonsymetric <- S
  S_nonsymetric[1, 2] <- S[1, 2] - 1
  expect_error(check_correctness_of_arguments(
    S_nonsymetric,
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))

  S_non_positive_semi_definite <- S - diag(ncol(S)) * eigen(S, symmetric = TRUE, only.values = TRUE)[["values"]][2]
  expect_error(check_correctness_of_arguments(
    S_non_positive_semi_definite,
    number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, NULL, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, 0, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations + 0.1, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30.1,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 1,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    "(1,3)(2,4)(5,6)",
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    gips_perm("(1,3)(2,4)(5,6)", 7),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    NULL, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_silent(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    1.1, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    0.9, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    1, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S) + 1), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, 7, FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, matrix(1:30, nrow = ncol(S)), FALSE, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), "FALSE", FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, "FALSE", FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, "FALSE", FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, FALSE, FALSE, "FALSE"
  ))

  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, diag(nrow = ncol(S)), FALSE, TRUE, FALSE, FALSE
  )) # return_probabilities can be TRUE only when save_all_perms is also TRUE
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, matrix(1:30, nrow = ncol(S)), NA, FALSE, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, matrix(1:30, nrow = ncol(S)), FALSE, NA, FALSE, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, matrix(1:30, nrow = ncol(S)), FALSE, FALSE, NA, FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S, number_of_observations, 30,
    permutations::permutation("(1,3)(2,4)(5,6)"),
    3, matrix(1:30, nrow = ncol(S)), FALSE, FALSE, FALSE, NA
  ))


  # A lot of problems at once
  expect_error(check_correctness_of_arguments(
    S, number_of_observations + 0.1, 1,
    "(1,3)(2,4)(5,6)",
    1, diag(nrow = ncol(S)), "FALSE", "FALSE", "FALSE", "FALSE"
  ), "8 problems identified with the provided arguments")

  # old tests:
  # A single problem at the same time:
  expect_error(check_correctness_of_arguments(
    S = NULL, number_of_observations, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(
    S = matrix(c(1:30), nrow = 6),
    number_of_observations, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = NULL, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = 0, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = 15.5, max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10.5,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 0,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 1,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = NULL, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = NULL, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 1, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = 7,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3,
    D_matrix = matrix(c(1:30), nrow = 6),
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3,
    D_matrix = NULL,
    was_mean_estimated = "FALSE", return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3,
    D_matrix = NULL,
    was_mean_estimated = NA, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = 7,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = NULL,
    save_all_perms = FALSE, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = "FALSE", show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = NA, show_progress_bar = FALSE
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = 7
  ))
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm, number_of_observations,
    max_iter = 10,
    start_perm = example_perm, delta = 3, D_matrix = NULL,
    was_mean_estimated = FALSE, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = NULL
  ))


  # A number of problems at the same time. Not all are printed:
  expect_error(check_correctness_of_arguments(matrix_invariant_by_example_perm,
    number_of_observations = -1, max_iter = 1,
    start_perm = example_perm, delta = 1, D_matrix = 7,
    was_mean_estimated = NA, return_probabilities = "FALSE",
    save_all_perms = 7, show_progress_bar = NULL
  ), "\\.\\.\\. and 3 more problems")
})


test_that("print.gips() works", {
  expect_identical(
    convert_log_diff_to_str(1009.5, 3),
    "2.632e+438"
  )
  expect_identical(
    convert_log_diff_to_str(16.1, 3),
    "9820670.922"
  )
  expect_identical(
    convert_log_diff_to_str(16.2, 3),
    "1.085e+7"
  )
  expect_identical(
    convert_log_diff_to_str(-7.677, 3),
    "4.634e-4"
  )
  expect_identical(
    convert_log_diff_to_str(Inf, 3),
    "Inf"
  )
  expect_identical(
    convert_log_diff_to_str(-Inf, 3),
    "-Inf"
  )
  expect_identical(
    convert_log_diff_to_str(0, 3),
    "1"
  )

  g <- gips(S, number_of_observations, was_mean_estimated = FALSE)
  g_map <- find_MAP(g, 10, show_progress_bar = FALSE, optimizer = "MH")

  expect_output(
    print(g),
    "The permutation \\(\\):\n - is 1 times more likely than the \\(\\) permutation"
  )
  expect_output(
    print(g, log_value = TRUE),
    "The permutation \\(\\):\n - is 1 times more likely than the \\(\\) permutation;\n - has log posteriori"
  )
  expect_output(
    print(g_map),
    "\n - was found after 10 posteriori calculations;\n - is"
  )

  # oneline:
  expect_output(
    print(g, oneline = TRUE),
    "The permutation \\(\\): is 1 times more likely"
  )
  expect_output(
    print(g, oneline = TRUE, log_value = TRUE),
    "permutation; has log posteriori "
  )
  expect_output(
    print(g_map, oneline = TRUE),
    ": was found after 10 posteriori calculations; is"
  )
})

test_that("plot.gips() works or abords for wrong arguments", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
    was_mean_estimated = FALSE, perm = custom_perm1
  )

  expect_error(plot.gips(custom_perm1))
  expect_error(plot(g1, type = "both"))
  expect_error(plot(g1, type = c("du", "pa")))
  expect_message(
    plot(g1),
    "`type = NA` was automatically changed to `type = 'heatmap'`"
  )
  expect_silent(my_ggplot1 <- plot(g1, type = "heatmap"))
  expect_silent(my_ggplot2 <- plot(g1, type = "MLE"))
  expect_true(all.equal(my_ggplot1, my_ggplot2)) # cannot use expect_equal(), because it checks for equal enviroments, but in R all enviroments are different even if have the same elements in itself

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

test_that("plot.gips() works for books example", {
  Z <- DAAG::oddbooks[, c(1, 2, 3)]
  Z$height <- Z$height / sqrt(2)

  g <- gips(cov(Z), 7, D_matrix = 1 * diag(3))
  expect_silent(plot(g, type = "heatmap"))
})

test_that("get_diagonalized_matrix_for_heatmap() works", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5)(6)", 6)
  g1 <- gips(S, number_of_observations,
    was_mean_estimated = FALSE, perm = custom_perm1
  )
  actual <- get_diagonalized_matrix_for_heatmap(g1)
  # block_ends: 3,5,6
  actual[1:3, 1:3] <- NA
  actual[4:5, 4:5] <- NA
  actual[6, 6] <- NA
  expect_true(all(is.na(actual)))
})

test_that("summary.gips() works", {
  custom_perm1 <- gips_perm("(1,2)(3,4,5,6)", 6)
  g1 <- gips(S, number_of_observations,
    was_mean_estimated = FALSE, perm = custom_perm1,
    D_matrix = diag(1, 6)
  )

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
    log_times_more_likely_than_id = start_permutation_log_posteriori - log_posteriori_id,
    n0 = 2, S_matrix = S, number_of_observations = 13,
    was_mean_estimated = FALSE,
    delta = 3, D_matrix = structure(c(
      1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1
    ), .Dim = c(6L, 6L)),
    n_parameters = 7,
    AIC = AIC(g1),
    BIC = BIC(g1)
  ), class = "summary.gips")

  expect_identical(summary(g1), my_sum) # all in my_sum is calculated by the same functions the summary should use

  expect_output(
    print(summary(g1)),
    "The number of observations is bigger than n0 for this permutation,\nso "
  )

  expect_output(
    print(summary(g1)), "free parameters"
  )

  expect_output(
    print(summary(g1)), "AIC"
  )

  expect_output(
    print(summary(g1)), "BIC"
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

  # Optimized with BF; those 3 are computed differently for optimized with BF:
  g4 <- find_MAP(gips(S[1:4, 1:4], number_of_observations), optimizer = "BF", show_progress_bar = FALSE)
  expect_true(is.null(summary(g4)[["when_was_best"]]))
  expect_true(is.null(summary(g4)[["log_posteriori_calls_after_best"]]))
  expect_s3_class(summary(g4)[["start_permutation"]], "gips_perm")

  # proper log_posteriori_calls_after_best:
  # g_fake started in (1,2,3) and then comed to (1,3,2), which had bigger log_posteriori value, but only a little bigger
  g_fake <- structure(list(structure(list(c(1, 2, 3)), size = 3L, class = "gips_perm")), S = structure(c(
    1, 0, 0, 0, 1, 0, 0, 0, 1
  ), dim = c(3L, 3L)), number_of_observations = 10, delta = 3, D_matrix = structure(c(
    1, 0, 0, 0, 1, 0, 0, 0, 1
  ), dim = c(3L, 3L)), was_mean_estimated = TRUE, optimization_info = list(
    original_perm = gips_perm("", 3),
    acceptance_rate = 0.333333333333333, log_posteriori_values = c(
      -16.01209771488621000000, -18.9125198086359,
      -16.01209771488621421085
    ), visited_perms = list(
      structure(list(c(1, 2, 3)), size = 3L, class = "gips_perm"),
      structure(list(1, c(2, 3)), size = 3L, class = "gips_perm"),
      structure(list(c(1, 3, 2)), size = 3L, class = "gips_perm")
    ),
    start_perm = structure(list(c(1, 2, 3)), size = 3L, class = "gips_perm"),
    last_perm = structure(list(c(1, 3, 2)), size = 3L, class = "gips_perm"),
    last_perm_log_posteriori = -16.0120977148862, iterations_performed = 2L,
    optimization_algorithm_used = "Metropolis_Hastings", post_probabilities = NULL,
    did_converge = NULL, best_perm_log_posteriori = -16.0120977148862,
    optimization_time = structure(0.00564193725585938, class = "difftime", units = "secs"),
    whole_optimization_time = structure(0.00564193725585938, class = "difftime", units = "secs")
  ), class = "gips")

  expect_equal(
    summary(g_fake)[["log_posteriori_calls_after_best"]],
    2
  )
})

test_that("start_permutation_log_posteriori was calculated correctly", {
  g <- gips(S, number_of_observations, was_mean_estimated = TRUE)
  g_map <- find_MAP(g,
    max_iter = 10, show_progress_bar = FALSE,
    optimizer = "MH"
  )

  optimization_info <- attr(g_map, "optimization_info")
  expect_equal(
    optimization_info[["log_posteriori_values"]][1],
    log_posteriori_of_perm(
      "",
      attr(g_map, "S"),
      attr(g_map, "number_of_observations") - attr(g_map, "was_mean_estimated"),
      attr(g_map, "delta"),
      attr(g_map, "D_matrix")
    )
  )
})

test_that("summary.gips() returns proper n0 for estimated and unestimated mean", {
  g_no_em <- gips(S, number_of_observations, was_mean_estimated = FALSE)
  g_em <- gips(S, number_of_observations, was_mean_estimated = TRUE)

  expect_equal(summary.gips(g_no_em)[["n0"]], ncol(S)) # for known mean and perm id, one needs n >= p
  expect_equal(summary.gips(g_em)[["n0"]], ncol(S) + 1) # for estimated mean and perm id, one needs n >= (p + 1)
})

test_that("summary.gips() returns proper Times more likely than identity permutation", {
  g_no_em <- gips(S, number_of_observations, was_mean_estimated = FALSE)
  g_em <- gips(S, number_of_observations, was_mean_estimated = TRUE)

  expect_equal(
    summary.gips(g_no_em)[["times_more_likely_than_id"]],
    1 # for known mean
  )
  expect_equal(
    summary.gips(g_em)[["times_more_likely_than_id"]],
    1 # for estimated mean
  )
})

test_that("print.summary.gips() will properly print BIC and AIC when n < n0", {
  p <- 6
  number_of_observations <- 2
  Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix) # this place is random
  S <- (t(Z) %*% Z) / number_of_observations # the theoretical mean is 0
  
  g <- gips(S, number_of_observations, perm = "()")
  expect_silent(my_sum <- summary.gips(g))
  expect_null(my_sum$AIC)
  expect_null(my_sum$BIC)
  
  expect_output(print(my_sum), "BIC:\n The number of observations is smaller than n0 for this permutation,\n so the gips model based on the found permutation does not exist.")
  expect_output(print(my_sum), "AIC:\n The number of observations is smaller than n0 for this permutation,\n so the gips model based on the found permutation does not exist.")
})


test_that("get_probabilities_from_gips() works", {
  g <- gips(matrix(c(1, 0.5, 0.5, 1.3), nrow = 2), 13, was_mean_estimated = FALSE)
  g_map <- find_MAP(g,
    optimizer = "BF", show_progress_bar = FALSE,
    return_probabilities = TRUE, save_all_perms = TRUE
  )

  expect_silent(probs <- get_probabilities_from_gips(g_map))
  expect_type(probs, "double")
  expect_named(probs)

  expect_error(
    get_probabilities_from_gips(g),
    "Did You forget to optimize `g`?"
  )

  g_map_no_prob <- find_MAP(g,
    optimizer = "BF", show_progress_bar = FALSE,
    return_probabilities = FALSE, save_all_perms = TRUE
  )
  expect_message(
    out <- get_probabilities_from_gips(g_map_no_prob),
    "on the `gips` object that does not have saved probabilities."
  )
  expect_null(out)

  g_map_no_prob_no_save <- find_MAP(g,
    optimizer = "BF", show_progress_bar = FALSE,
    return_probabilities = FALSE, save_all_perms = FALSE
  )
  expect_message(
    out <- get_probabilities_from_gips(g_map_no_prob_no_save),
    "on the `gips` object that does not have saved probabilities."
  )
  expect_null(out)

  # sorted
  probs <- get_probabilities_from_gips(g_map)
  expect_equal(order(probs, decreasing = TRUE), c(1, 2))
})

test_that("forget_perms() works properly", {
  g <- gips(S, number_of_observations, was_mean_estimated = FALSE)
  g_map <- find_MAP(g, 10,
    show_progress_bar = FALSE, optimizer = "MH",
    save_all_perms = TRUE
  )

  expect_silent(validate_gips(g_map))
  expect_message(
    forget_perms(g),
    "Provided \\`g\\` is a \\`gips\\` object, but it was not optimized yet\\."
  )
  expect_false(all(is.na(attr(g_map, "optimization_info")[["visited_perms"]])))
  expect_silent(g_map_no_perms <- forget_perms(g_map))
  expect_true(all(is.na(attr(g_map_no_perms, "optimization_info")[["visited_perms"]])))
  expect_silent(validate_gips(g_map_no_perms))

  expect_message(
    g_map_no_perms_again <- forget_perms(g_map_no_perms),
    "Provided \\`g\\` is an optimized \\`gips\\` object that already has forgotten all permutations\\."
  )
})

test_that("logLik.gips() works", {
  # logLik calculated by hand:
  Z <- matrix(c(
    -1, 2, -3, 4,
    1, -1, 1, 1,
    -1, 2, 2, 3,
    2, -3, 2, -3,
    3, -2, 2, -3
  ), nrow = 5, byrow = TRUE)

  p <- ncol(Z) # 4
  n <- nrow(Z) # 5


  # ==================
  # Mean is (0,0,0,0)
  U <- t(Z) %*% Z
  perm <- gips_perm("(12)(34)", 4)
  S <- project_matrix(U, perm) / n

  # logLik from definition:
  loglikelihoods <- mvtnorm::dmvnorm(
    Z, rep(0, 4), S, log = TRUE
  )
  logLik_definition <- sum(loglikelihoods)

  expect_equal(logLik_definition, -35.2883973048347)

  attr(logLik_definition, "df") <- 6 # (choose(p, 2) + p) - 4 parameters; choose(p, 2) + p is standard CoV; 4 is how much equalities is with (12)(34)
  attr(logLik_definition, "nobs") <- n
  class(logLik_definition) <- "logLik"

  # logLik.gips:
  expect_equal(
    logLik(gips(U / n, n, perm = perm, was_mean_estimated = FALSE)),
    logLik_definition
  )


  # ==================
  # mean was estimated
  U <- cov(Z) * (n - 1)
  perm <- gips_perm("(12)(34)", 4)
  S <- project_matrix(U, perm) / (n - 1)

  logLik_expected <- structure(-28.8015774226105, df = 6, nobs = 5L, class = "logLik")
  
  # logLik.gips:
  expect_equal(
    logLik(gips(U / (n - 1), n, perm = perm, was_mean_estimated = TRUE)),
    logLik_expected
  )


  # ==================
  # NUll:
  U <- t(Z) %*% Z
  expect_warning(
    expect_equal(
      logLik(gips(U / n, 2, perm = perm)), NULL
    ),
    class = "likelihood_does_not_exists"
  )

  # ==================
  # -Inf:
  g <- gips(diag(0, 4), n)
  expect_warning(
    expect_equal(logLik(g), -Inf),
    class = "singular_matrix"
  )
})

test_that("AIC.gips() works", {
  Z <- matrix(c(
    -1, 2, -3, 4,
    1, -1, 1, 1,
    -1, 2, 2, 3,
    2, -3, 2, -3,
    3, -2, 2, -3
  ), nrow = 5, byrow = TRUE)

  p <- ncol(Z) # 4
  n <- nrow(Z) # 5

  g <- gips(t(Z) %*% Z / n, n, perm = "(12)(34)", was_mean_estimated = FALSE)

  expect_equal(AIC(g), 82.5767946096694)
  expect_equal(BIC(g), 80.233422084274)

  # ==================
  # NUll
  S <- t(Z) %*% Z / n
  g <- gips(S, 2)
  expect_warning(
    expect_equal(
      AIC(g), NULL
    ),
    class = "likelihood_does_not_exists"
  )
  expect_warning(
    expect_equal(
      BIC(g), NULL
    ),
    class = "likelihood_does_not_exists"
  )

  # ==================
  # Inf
  g <- gips(diag(0, 4), n)
  expect_warning(expect_equal(AIC(g), Inf), class = "singular_matrix")
  expect_warning(expect_equal(BIC(g), Inf), class = "singular_matrix")
})

test_that("as.character.gips() work", {
  A <- matrix(rnorm(4 * 4), nrow = 4)
  S <- t(A) %*% A
  g <- gips(S, 14, perm = "(123)")
  expect_equal(as.character(g), "(1,2,3)")
})


test_that("print.gips() will show the original perm for BF", {
  g <- gips(S[1:4,1:4], number_of_observations, perm = "(123)")
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  expect_output(print(g_map), "than the \\(1,2,3\\)")
})

test_that("print.summary.gips() will not compare with original the unoptimized gips that is in id", {
  g <- gips(S[1:4,1:4], number_of_observations, perm = "(123)")
  expect_output(print(summary(g)), "Times more likely than identity permutation", fixed = TRUE)
  
  g <- gips(S[1:4,1:4], number_of_observations, perm = "()")
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  expect_output(print(summary(g_map)), "Times more likely than starting permutation", fixed = TRUE)
  
  pattern <- "Log_posteriori:\\s*(-?\\d+\\.\\d+)\\s\\sThe number of observations"
  expect_output(print(summary(g)), pattern) # The "Times more likely than starting permutation:" is skipped
})
