#' Find the Maximum A Posteriori Estimation
#'
#' Use one of the optimization algorithms to find the permutation that maximizes
#' a posteriori based on observed data. Not all optimization algorithms will
#' always find the MAP, but they try to find a significant value.
#' More information can be found in the
#' 'Possible algorithms to use as optimizers' section below.
#'
#' @section Possible algorithms to use as optimizers:
#'
#' * `"Metropolis_Hastings"`, `"MH"` - to use Metropolis-Hastings algorithm;
#'     [see Wikipedia](https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm).
#'     The algorithm will draw a random transposition in every iteration
#'     and consider changing the current one.
#'     When the `max_iter` is reached, the algorithm will return
#'     the MAP Estimator as the best permutation calculated so far.
#'
#' * `"hill_climbing"`, `"HC"` - to use hill climbing algorithm;
#'     [see Wikipedia](https://en.wikipedia.org/wiki/Hill_climbing).
#'     The algorithm will check all transpositions in every iteration and
#'     go to the one with the biggest a posteriori value.
#'     The optimization ends when all *neighbors* will have a smaller
#'     a posteriori value. If the `max_iter` is reached before the end,
#'     then the warning is shown, and it is recommended to start
#'     the optimization again on the output of the `find_MAP()`.
#'     Remember that there are `p*(p-1)/2` transpositions to be checked
#'     in every iteration. For bigger `p`, this may be costly.
#'
#' * `"brute_force"`, `"BF"`, `"full"` - to use the Brute Force algorithm that
#'     checks the whole permutation space of a given size. This algorithm will
#'     definitely find the Maximum A Posteriori Estimation but is very
#'     computationally expensive for bigger space.
#'
#' @param g Object of a `gips` class
#' @param max_iter Number of iterations for an algorithm to perform.
#'     At least 2. For `optimizer=="MH"` it has to be finite;
#'     for `optimizer=="HC"` it can be infinite;
#'     for `optimizer=="BF"` it is not used.
#' @param return_probabilities A boolean. TRUE can only be provided
#'     for `optimizer=="MH"`. Whether to use Metropolis-Hastings results to
#'     calculate posterior probabilities.
#' @param show_progress_bar A boolean.
#'     Indicate whether or not to show the progress bar.
#'     When `max_iter` is infinite, `show_progress_bar` has to be `FALSE`.
#' @param optimizer The optimizer for the search of the maximum posteriori.
#'   * `"MH"` (the default for unoptimized `g`) - Metropolis-Hastings
#'   * `"HC"` - Hill Climbing
#'   * `"BF"` - Brute Force
#'   * `"continue"` (the default for optimized `g`) - The same as
#'       the `g` was optimized by (see Examples).
#'
#' For more details, see the "Possible algorithms to use as optimizers"
#' section below.
#'
#' @returns Returns an optimized object of a `gips` class.
#'
#' @export
#'
#' @seealso
#' * [gips()] - The constructor of a `gips` class.
#'     The `gips` object is used as the `g` parameter.
#' * [plot.gips()] - Practical plotting function for
#'     visualizing the optimization process.
#' * [summary.gips()] - The function that summarizes the output of optimization.
#' * [log_posteriori_of_gips()] - The function that the optimizers
#'     of `find_MAP()` tries to find the argmax of.
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 6
#' mu <- runif(6, -10, 10) # Assume we don't know the mean
#' sigma_matrix <- matrix(
#'   data = c(
#'     1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
#'     0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
#'     0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
#'     0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
#'     0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
#'     0.8, 0.6, 0.4, 0.6, 0.8, 1.0
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' ) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- cov(Z) # Assume we have to estimate the mean
#'
#' g <- gips(S, number_of_observations)
#'
#' g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "MH")
#' g_map
#'
#' g_map2 <- find_MAP(g_map, max_iter = 10, show_progress_bar = FALSE, optimizer = "continue")
#'
#' if (require("graphics")) {
#'   plot(g_map2, type = "both", logarithmic_x = TRUE)
#' }
#'
#' g_map_BF <- find_MAP(g, show_progress_bar = FALSE, optimizer = "BF")
#' summary(g_map_BF)
find_MAP <- function(g, max_iter = NA, return_probabilities = FALSE,
                     show_progress_bar = TRUE, optimizer = NA) {
  # check the correctness of the g argument
  validate_gips(g)

  possible_optimizers <- c(
    "MH", "Metropolis_Hastings", "HC", "hill_climbing",
    "BF", "brute_force", "full", "continue"
  )

  # check the correctness of the rest of arguments
  if (length(optimizer) > 1) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = paste0(
        "`optimizer` must be the character vector of length 1. Must be one of: c('",
        paste0(possible_optimizers, collapse = "', '"), "')."
      ),
      "x" = paste0(
        "You provided `optimizer == (",
        paste0(optimizer, collapse = ", "), ")`."
      ),
      "i" = "Did You misspelled the optimizer name?"
    ))
  }
  # default optimizer:
  if (is.na(optimizer)) {
    optimizer <- ifelse(is.null(attr(g, "optimization_info")),
      "MH", "continue"
    )

    rlang::inform(c("You used the default value of the 'optimizer' argument in `find_MAP()`.",
      "i" = paste0(
        "The 'optimizer = NA' was automatically changed to 'optimizer = \"",
        optimizer, "\"'."
      )
    ))
  }

  # get a chosen optimizer, even with part of the name:
  chosen_optimizer_number <- pmatch(optimizer, possible_optimizers)

  if (is.na(chosen_optimizer_number)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = paste0(
        "`optimizer` must be one of: c('",
        paste0(possible_optimizers, collapse = "', '"), "')."
      ),
      "x" = paste0("You provided `optimizer == '", optimizer, "'`."),
      "i" = "Did You misspelled the optimizer name?"
    ))
  }

  if (optimizer != possible_optimizers[chosen_optimizer_number]) {
    rlang::inform(c(
      "You provided a shortcut for the optimization method's name:",
      "i" = paste0("You provided `optimizer == '", optimizer, "'`"),
      "i" = paste0(
        "This will be changed to `optimizer == '",
        possible_optimizers[chosen_optimizer_number], "'`"
      )
    ))

    optimizer <- possible_optimizers[chosen_optimizer_number]
  }

  if (!(optimizer %in% c("BF", "brute_force", "full")) &&
    is.na(max_iter)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "`max_iter = NA` can be provided only for `optimizer` one of: c('BF', 'brute_force', 'full'). For any other, `max_iter` must be a whole number, strictly bigger than 1.",
      "x" = paste0("You provided `optimizer == ", optimizer, "` and `max_iter = NA`."),
      "i" = "Did You forgot to set the `max_iter`?",
      "i" = "Did You misspelled the optimizer name?"
    ))
  }

  continue_optimization <- (optimizer == "continue")
  if (continue_optimization) {
    if (is.null(attr(g, "optimization_info"))) {
      rlang::abort(c("There was a problem identified with provided arguments:",
        "i" = "`optimizer == 'continue'` can be provided only with optimized gips object `g`.",
        "x" = "You provided `optimizer == 'continue'`, but the gips object `g` is not optimized.",
        "i" = "Did You provided wrong `gips` object?",
        "i" = "Did You want to call another optimizer like 'MH' or 'HC'?"
      ))
    }

    optimizer <- attr(g, "optimization_info")[["optimization_algorithm_used"]][length(attr(g, "optimization_info")[["optimization_algorithm_used"]])] # this is the last used optimizer
    if (optimizer == "brute_force") {
      rlang::abort(c("There was a problem identified with provided arguments:",
        "i" = "`optimizer == 'continue'` cannot be provided after optimizating with `optimizer == 'brute_force'`, because the whole space was already browsed.",
        "x" = "You provided `optimizer == 'continue'`, but the gips object `g` was optimized with brute_force optimizer. Better permutation will not be found."
      ))
    }
  }

  if (!(optimizer %in% c("MH", "Metropolis_Hastings")) && return_probabilities) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "Probabilities can only be returned with the `optimizer == 'Metropolis_Hastings'`",
      "x" = "You provided both `optimizer != Metropolis_Hastings` and `return_probabilities == TRUE`!",
      "i" = "Did You want to use `optimizer == Metropolis_Hastings` or `return_probabilities == FLASE`?"
    ))
  }

  if (return_probabilities) {
    rlang::check_installed("stringi",
      reason = "to return probabilities in `find_MAP(optimizer = 'Metropolis_Hastings', return_probabilities = TRUE)`; without this package, probabilities cannot be returned"
    )
    if (!rlang::is_installed("stringi")) {
      rlang::warn(c("There was a problem with return_probabilities:",
        "i" = "Package `stringi` is required to successfully call `find_MAP(optimizer = 'Metropolis_Hastings', return_probabilities = TRUE)`.",
        "x" = "You do not have package `stringi` installed.",
        "i" = "Optimization will proceed as `find_MAP(optimizer = 'Metropolis_Hastings', return_probabilities = FALSE)`."
      ))

      return_probabilities <- FALSE
    }
  }

  # inform that user can consider "BF"
  if ((optimizer %in% c("MH", "Metropolis_Hastings")) &&
    (max_iter * 10 >= prod(1:ncol(attr(g, "S")))) &&
    is.finite(max_iter)) { # infinite max_iter is illegal, but additional check will not hurt
    rlang::inform(c(paste0(
      "You called optimization with Metropolis_Hastings algorith with ",
      max_iter, " iterations."
    ),
    "i" = paste0(
      "Consider using `optimizer = 'brute_force'`, because it will use ",
      ncol(attr(g, "S")), "! (factorial) = ", prod(1:ncol(attr(g, "S"))),
      " iterations and will browse all permutations, therefore it will definitely find the maximum posteriori estimator."
    )
    ))
  }

  # extract parameters
  S <- attr(g, "S")
  number_of_observations <- attr(g, "number_of_observations")
  if (continue_optimization) { # the `ifelse()` function cannot be used because the objects are lists
    start_perm <- attr(g, "optimization_info")[["last_perm"]]
  } else {
    start_perm <- g[[1]]
  }
  delta <- attr(g, "delta")
  D_matrix <- attr(g, "D_matrix")
  was_mean_estimated <- attr(g, "was_mean_estimated")
  
  if(was_mean_estimated) { # one degree of freedom is lost; we will return this 1 to number_of_observations after optimization in `combine_gips()`
    edited_number_of_observations <- number_of_observations - 1
  } else {
    edited_number_of_observations <- number_of_observations
  }

  start_time <- Sys.time()

  if (optimizer %in% c("MH", "Metropolis_Hastings")) {
    gips_optimized <- Metropolis_Hastings(
      S = S, number_of_observations = edited_number_of_observations,
      max_iter = max_iter, start_perm = start_perm,
      delta = delta, D_matrix = D_matrix,
      return_probabilities = return_probabilities,
      show_progress_bar = show_progress_bar
    )
  } else if (optimizer %in% c("HC", "hill_climbing")) {
    gips_optimized <- hill_climbing(
      S = S, number_of_observations = edited_number_of_observations,
      max_iter = max_iter, start_perm = start_perm,
      delta = delta, D_matrix = D_matrix,
      show_progress_bar = show_progress_bar
    )
  } else if (optimizer %in% c("BF", "brute_force", "full")) {
    gips_optimized <- brute_force(
      S = S, number_of_observations = edited_number_of_observations,
      delta = delta, D_matrix = D_matrix,
      show_progress_bar = show_progress_bar
    )
  }

  end_time <- Sys.time()
  attr(gips_optimized, "optimization_info")[["optimization_time"]] <- end_time - start_time
  attr(gips_optimized, "optimization_info")[["full_optimization_time"]] <- end_time - start_time
  
  structure_constants <- get_structure_constants(gips_optimized[[1]])
  n0 <- max(structure_constants[["r"]] * structure_constants[["d"]] / structure_constants[["k"]])
  if (n0 > number_of_observations) {
    rlang::warn(c("The found permutation has n0 value bigger than the number of observations of the normal variable.",
      "i" = "The covariance matrix invariant under the found permutation does not have the likelihood properly defined.",
      "i" = "For more in-depth explanation, see the Vignette (TODO)."
    ))
  }
  

  return(combine_gips(g, gips_optimized))
}


Metropolis_Hastings <- function(S, number_of_observations, max_iter, start_perm = NULL,
                                delta = 3, D_matrix = NULL, return_probabilities = FALSE,
                                show_progress_bar = TRUE) {
  if (is.null(start_perm)) {
    start_perm <- permutations::id
  }

  check_correctness_of_arguments(
    S = S, number_of_observations = number_of_observations,
    max_iter = max_iter, start_perm = start_perm,
    delta = delta, D_matrix = D_matrix, was_mean_estimated = FALSE,
    return_probabilities = return_probabilities,
    show_progress_bar = show_progress_bar
  )

  if (!inherits(start_perm, "gips_perm")) {
    start_perm <- gips_perm(start_perm, nrow(S)) # now we know the `S` is a matrix
  }

  if (is.infinite(max_iter)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "`max_iter` in `Metropolis_Hastings()` must be finite.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    ))
  }

  perm_size <- dim(S)[1]
  if (permutations::is.cycle(start_perm)) {
    start_perm <- gips_perm(start_perm, perm_size)
  }
  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size)
  }

  my_goal_function <- function(perm) {
    log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `number_of_observations` if the mean was estimated!
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix
    )
  }

  acceptance <- rep(FALSE, max_iter)
  log_posteriori_values <- rep(0, max_iter)
  visited_perms <- list()
  visited_perms[[1]] <- start_perm

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  }
  log_posteriori_values[1] <- my_goal_function(visited_perms[[1]])

  found_perm <- start_perm
  found_perm_log_posteriori <- log_posteriori_values[1]

  Uniformly_drawn_numbers <- stats::runif(max_iter, min = 0, max = 1)

  # main loop
  for (i in 1:(max_iter - 1)) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }

    e <- runif_transposition(perm_size)
    perm_proposal <- compose_with_transposition(visited_perms[[i]], e)

    goal_function_perm_proposal <- my_goal_function(perm_proposal)
    if (is.nan(goal_function_perm_proposal) || is.infinite(goal_function_perm_proposal)) {
      # See ISSUE#5; We hope the introduction of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(goal_function_perm_proposal), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))

      break()
    }

    # if goal_function_perm_proposal > log_posteriori_values[i], then it is true, because Uniformly_drawn_numbers[i] \in [0,1]
    if (Uniformly_drawn_numbers[i] < exp(goal_function_perm_proposal - log_posteriori_values[i])) { # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2). That means this is Metropolis algorithm, not necessary Metropolis-Hastings.
      visited_perms[[i + 1]] <- perm_proposal
      log_posteriori_values[i + 1] <- goal_function_perm_proposal
      acceptance[i] <- TRUE

      if (found_perm_log_posteriori < log_posteriori_values[i + 1]) {
        found_perm_log_posteriori <- log_posteriori_values[i + 1]
        found_perm <- visited_perms[[i + 1]]
      }
    } else {
      visited_perms[[i + 1]] <- visited_perms[[i]]
      log_posteriori_values[i + 1] <- log_posteriori_values[i] # TODO(Do we really want to forget the calculated values? The algorithm HC works differently)
    }
  }

  if (show_progress_bar) {
    close(progressBar)
  }

  function_calls <- length(log_posteriori_values)

  if (return_probabilities) {
    probabilities <- estimate_probabilities(visited_perms)
  } else {
    probabilities <- NULL
  }
  optimization_info <- list(
    "acceptance_rate" = mean(acceptance),
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = visited_perms,
    "last_perm" = visited_perms[[function_calls]],
    "last_perm_log_posteriori" = log_posteriori_values[function_calls],
    "iterations_performed" = i,
    "optimization_algorithm_used" = "Metropolis_Hastings",
    "post_probabilities" = probabilities,
    "did_converge" = NULL,
    "best_perm_log_posteriori" = found_perm_log_posteriori,
    "optimization_time" = NA,
    "full_optimization_time" = NA
  )


  new_gips(
    list(found_perm), S, number_of_observations,
    delta, D_matrix, was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}


hill_climbing <- function(S, number_of_observations, max_iter = 5,
                          start_perm = NULL,
                          delta = 3, D_matrix = NULL,
                          show_progress_bar = TRUE) {
  if (is.null(start_perm)) {
    start_perm <- permutations::id
  }

  check_correctness_of_arguments(
    S = S, number_of_observations = number_of_observations,
    max_iter = max_iter, start_perm = start_perm,
    delta = delta, D_matrix = D_matrix, was_mean_estimated = FALSE,
    return_probabilities = FALSE,
    show_progress_bar = show_progress_bar
  )

  if (!inherits(start_perm, "gips_perm")) {
    start_perm <- gips_perm(start_perm, nrow(S)) # now we know the `S` is a matrix
  }

  if (show_progress_bar && is.infinite(max_iter)) {
    stop("Progress bar is not yet supported for infinite max_iter. Rerun the algorithm with show_progress_bar=FALSE or finite max_iter. For more information see ISSUE#8.") # See ISSUE#8
  }

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  }

  perm_size <- dim(S)[1]

  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size)
  }


  my_goal_function <- function(perm) {
    log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `number_of_observations` if the mean was estimated!
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix
    )
  }

  goal_function_best_logvalues <- numeric(0)
  log_posteriori_values <- numeric(0)

  # init
  speciments <- list()
  speciments[[1]] <- start_perm
  goal_function_best_logvalues[1] <- my_goal_function(speciments[[1]])
  log_posteriori_values[1] <- goal_function_best_logvalues[1]

  # mail loop
  iteration <- 0
  did_converge <- FALSE
  while (iteration <= max_iter - 1) {
    iteration <- iteration + 1
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, iteration)
    }

    best_neighbour <- NULL
    best_neighbour_value <- -Inf
    for (i in 1:(perm_size - 1)) {
      for (j in (i + 1):perm_size) {
        neighbour <- compose_with_transposition(speciments[[iteration]], c(i, j))
        neighbour_value <- my_goal_function(neighbour)
        log_posteriori_values[length(log_posteriori_values) + 1] <- neighbour_value

        if (neighbour_value > best_neighbour_value) {
          best_neighbour_value <- neighbour_value
          best_neighbour <- neighbour
        }
      }
    }

    if (best_neighbour_value > goal_function_best_logvalues[iteration]) {
      goal_function_best_logvalues[iteration + 1] <- best_neighbour_value
      speciments[[iteration + 1]] <- best_neighbour
    } else {
      did_converge <- TRUE
      break
    }
  }

  if (show_progress_bar) {
    close(progressBar)
  }

  if (!did_converge) {
    rlang::warn(c(paste0("Hill Climbing algorithm did not converge in ", iteration, " iterations!"), # now, iteration == max_iter
      "i" = "We recommend to run the `find_MAP(optimizer = 'continue')` on the acquired output."
    ))
    iteration <- iteration + 1 # the very first was the starting perm
  } else {
    goal_function_best_logvalues <- goal_function_best_logvalues[1:iteration]
    if (show_progress_bar) {
      print(paste0("Algorithm did converge in ", iteration, " iterations"))
    }
  }

  function_calls <- length(log_posteriori_values)

  optimization_info <- list(
    "acceptance_rate" = 1 / choose(perm_size, 2),
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = speciments,
    "last_perm" = speciments[[iteration]],
    "last_perm_log_posteriori" = goal_function_best_logvalues[iteration],
    "iterations_performed" = iteration,
    "optimization_algorithm_used" = "hill_climbing",
    "post_probabilities" = NULL,
    "did_converge" = did_converge,
    "best_perm_log_posteriori" = goal_function_best_logvalues[iteration],
    "optimization_time" = NA,
    "full_optimization_time" = NA
  )
  

  new_gips(
    list(speciments[[iteration]]), S, number_of_observations,
    delta, D_matrix, was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}


brute_force <- function(S, number_of_observations,
                        delta = 3, D_matrix = NULL,
                        show_progress_bar = TRUE) {
  check_correctness_of_arguments(
    S = S, number_of_observations = number_of_observations,
    max_iter = 5, start_perm = permutations::id, # max_iter, was_mean_estimated and start_perm are not important for optimization with brute_force
    delta = delta, D_matrix = D_matrix, was_mean_estimated = FALSE,
    return_probabilities = FALSE,
    show_progress_bar = show_progress_bar
  )

  perm_size <- dim(S)[1]

  if (perm_size > 30) { # TODO(Test and set smaller, reasonable value. 15?)
    rlang::abort(c("Optimizer 'brute_force' cannot browse such a big permutional space.",
      "x" = paste0(
        "You provided a space with size ", perm_size,
        "! (factorial), which has ", prod(1:perm_size),
        " elements."
      ),
      "i" = "Do You want to use other optimizer for such a big space? For example 'Metropolis_Hastings' or 'hill_climbing'?"
    ))
  }

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = prod(1:perm_size), initial = 1)
  }

  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size)
  }

  my_goal_function <- function(perm) {
    log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `number_of_observations` if the mean was estimated!
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix
    )
  }

  # main loop
  all_perms_list <- permutations::allperms(perm_size)
  all_perms_list <- permutations::as.cycle(all_perms_list)
  log_posteriori_values <- sapply(1:length(all_perms_list), function(i) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }
    this_perm <- permutations::cycle(list(all_perms_list[[i]]))
    my_goal_function(this_perm)
  })

  if (show_progress_bar) {
    close(progressBar)
  }

  best_perm <- gips_perm(permutations::cycle(list(all_perms_list[[which.max(log_posteriori_values)]])), perm_size)

  optimization_info <- list(
    "acceptance_rate" = NULL,
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = all_perms_list,
    "last_perm" = NULL,
    "last_perm_log_posteriori" = NULL,
    "iterations_performed" = prod(1:perm_size),
    "optimization_algorithm_used" = "brute_force",
    "post_probabilities" = NULL,
    "did_converge" = TRUE,
    "best_perm_log_posteriori" = log_posteriori_values[which.max(log_posteriori_values)],
    "optimization_time" = NA,
    "full_optimization_time" = NA
  )
  

  new_gips(
    list(best_perm), S, number_of_observations,
    delta, D_matrix, was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}



#' Combining 2 gips objects
#'
#' g2 was optimized with a single optimization method. g1 was potentially non-optimized or optimized once, or optimized multiple times.
#' If g2 was optimized with "brute_force", forget the g1.
#'
#' @noRd
combine_gips <- function(g1, g2) {
  # first, adjust the number of observations:
  attr(g2, "number_of_observations") <- attr(g1, "number_of_observations")
  attr(g2, "was_mean_estimated") <- attr(g1, "was_mean_estimated")
  
  if (is.null(attr(g1, "optimization_info")) ||
    attr(g2, "optimization_info")[["optimization_algorithm_used"]] == "brute_force") { # when brute_force was used, forget the initial optimization
    
    return(g2)
  }

  # g1 is also an effect of optimization.
  optimization_info1 <- attr(g1, "optimization_info")
  optimization_info2 <- attr(g2, "optimization_info")

  n1 <- length(optimization_info1[["log_posteriori_values"]])
  n2 <- length(optimization_info2[["log_posteriori_values"]])

  visited_perms <- c(optimization_info1[["visited_perms"]], optimization_info2[["visited_perms"]]) # WoW, one can use `c()` to combine lists!
  optimization_algorithm_used <- c(optimization_info1[["optimization_algorithm_used"]], optimization_info2[["optimization_algorithm_used"]])

  if (all(optimization_algorithm_used == "Metropolis_Hastings") &&
    !is.null(optimization_info1[["post_probabilities"]]) &&
    !is.null(optimization_info2[["post_probabilities"]])) {
    post_probabilities <- estimate_probabilities(visited_perms) # TODO(This can be combined more optimally, but I (Adam) think this will be rarely done nevertheless.)
  } else {
    post_probabilities <- NULL
  }

  optimization_info_new <- list(
    "acceptance_rate" = (n1 * optimization_info1[["acceptance_rate"]] + n2 * optimization_info2[["acceptance_rate"]]) / (n1 + n2),
    "log_posteriori_values" = c(optimization_info1[["log_posteriori_values"]], optimization_info2[["log_posteriori_values"]]),
    "visited_perms" = visited_perms,
    "last_perm" = optimization_info2[["last_perm"]],
    "last_perm_log_posteriori" = optimization_info2[["last_perm_log_posteriori"]],
    "iterations_performed" = c(optimization_info1[["iterations_performed"]], optimization_info2[["iterations_performed"]]),
    "optimization_algorithm_used" = optimization_algorithm_used,
    "post_probabilities" = post_probabilities,
    "did_converge" = optimization_info2[["did_converge"]],
    "best_perm_log_posteriori" = max(optimization_info1[["best_perm_log_posteriori"]], optimization_info2[["best_perm_log_posteriori"]]),
    "optimization_time" = c(optimization_info1[["optimization_time"]], optimization_info2[["optimization_time"]]),
    "full_optimization_time" = optimization_info1[["full_optimization_time"]] + optimization_info2[["full_optimization_time"]]
  )

  if (optimization_info1[["best_perm_log_posteriori"]] > optimization_info2[["best_perm_log_posteriori"]]) {
    g_out <- g1 # in the continuation the new best was not found
  } else {
    g_out <- g2
  }

  attr(g_out, "optimization_info") <- optimization_info_new

  g_out
}
