#' Find the Maximum A Posteriori Estimation
#'
#' Use one of the optimization algorithms to find the permutation that
#' maximizes a posteriori probability based on observed data.
#' Not all optimization algorithms will always find the MAP, but they try
#' to find a significant value. More information can be found in
#' the "**Possible algorithms to use as optimizers**" section below.
#'
#' `find_MAP` can produce a warning when:
#' * the optimizer "hill_climbing" gets to the end of
#'   its `max_iter` without converging.
#' * the optimizer will find the permutation with smaller `n0` than
#'   `number_of_observations` (for more information on what it means,
#'   see **\eqn{C\sigma} and `n0`** section
#'   in `vignette("Theory", package = "gips")` or in its
#'   [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'
#' @section Possible algorithms to use as optimizers:
#'
#' For a in-depth explanation, see in
#'   `vignette("Optimizers", package = "gips")` or in its
#'   [pkgdown page](https://przechoj.github.io/gips/articles/Optimizers.html).
#'
#' For every algorithm, there are some aliases available.
#'
#' * `"Metropolis_Hastings"`, `"MH"` - use
#'     the **Metropolis-Hastings** algorithm;
#'     [see Wikipedia](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm).
#'     The algorithm will draw a random transposition in every iteration
#'     and consider changing the current state (permutation).
#'     When the `max_iter` is reached, the algorithm will return the best
#'     permutation calculated as the MAP Estimator. This implements
#'     the [*Second approach* from references, section 4.1.2](https://arxiv.org/abs/2004.03503).
#'     This algorithm used in this context is a special case of the
#'     **Simulated Annealing** the user may be more familiar with;
#'     [see Wikipedia](https://en.wikipedia.org/wiki/Simulated_annealing).
#'
#' * `"hill_climbing"`, `"HC"` - use
#'     the **hill climbing** algorithm;
#'     [see Wikipedia](https://en.wikipedia.org/wiki/Hill_climbing).
#'     The algorithm will check all transpositions in every iteration and
#'     go to the one with the biggest a posteriori value.
#'     The optimization ends when all *neighbors* will have a smaller
#'     a posteriori value. If the `max_iter` is reached before the end,
#'     then the warning is shown, and it is recommended to continue
#'     the optimization on the output of the `find_MAP()` with
#'     `optimizer = "continue"`; see examples.
#'     Remember that `p*(p-1)/2` transpositions will be checked
#'     in every iteration. For bigger `p`, this may be costly.
#'
#' * `"brute_force"`, `"BF"`, `"full"` - use
#'     the **Brute Force** algorithm that checks the whole permutation
#'     space of a given size. This algorithm will find
#'     the actual Maximum A Posteriori Estimation, but it is
#'     very computationally expensive for bigger spaces.
#'     We recommend Brute Force only for `p <= 9`.
#'     For the time the Brute Force takes on our machines, see in
#'     `vignette("Optimizers", package = "gips")` or in its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Optimizers.html).
#'
#' @param g Object of a `gips` class.
#' @param max_iter The number of iterations for an algorithm to perform.
#'     At least 2. For `optimizer="MH"`, it has to be finite;
#'     for `optimizer="HC"`, it can be infinite;
#'     for `optimizer="BF"`, it is not used.
#' @param optimizer The optimizer for the search of the maximum posteriori:
#'   * `"MH"` (the default for unoptimized `g`) - Metropolis-Hastings;
#'   * `"HC"` - Hill Climbing;
#'   * `"BF"` - Brute Force;
#'   * `"continue"` (the default for optimized `g`) - The same as
#'       the `g` was optimized by (see Examples).
#'
#' See the **Possible algorithms to use as optimizers**
#' section below for more details.
#' @param show_progress_bar A boolean.
#'     Indicate whether or not to show the progress bar:
#'   * When `max_iter` is infinite, `show_progress_bar` has to be `FALSE`;
#'   * When `return_probabilities=TRUE`, then
#'       shows an additional progress bar for the time
#'       when the probabilities are calculated.
#' @param save_all_perms A boolean. `TRUE` indicates saving
#'     a list of all permutations visited during optimization.
#'     This can be useful sometimes but need a lot more RAM.
#' @param return_probabilities A boolean. `TRUE` can only be provided
#'     only when `save_all_perms = TRUE`. For:
#'   * `optimizer="MH"` - use Metropolis-Hastings results to
#'       estimate posterior probabilities;
#'   * `optimizer="BF"` - use brute force results to
#'       calculate exact posterior probabilities.
#'
#' These additional calculations are costly, so a second progress bar
#'     is shown (when `show_progress_bar = TRUE`).
#'
#' To examine probabilities after optimization,
#'     call [get_probabilities_from_gips()].
#'
#' @returns Returns an optimized object of a `gips` class.
#'
#' @export
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [arXiv link](https://arxiv.org/abs/2004.03503);
#' \doi{10.1214/22-AOS2174}
#'
#' @seealso
#' * [gips()] - The constructor of a `gips` class.
#'     The `gips` object is used as the `g` parameter of `find_MAP()`.
#' * [plot.gips()] - Practical plotting function for
#'     visualizing the optimization process.
#' * [summary.gips()] - Summarize the output of optimization.
#' * [AIC.gips()], [BIC.gips()] - Get the Information Criterion
#'     of the found model.
#' * [get_probabilities_from_gips()] - When
#'     `find_MAP(return_probabilities = TRUE)` was called,
#'     probabilities can be extracted with this function.
#' * [log_posteriori_of_gips()] - The function that the optimizers
#'     of `find_MAP()` tries to find the argmax of.
#' * [forget_perms()] - When the `gips` object was optimized
#'     with `find_MAP(save_all_perms = TRUE)`, it will be of
#'     considerable size in RAM. `forget_perms()` can make such an object
#'     lighter in memory by forgetting the permutations it considered.
#' * `vignette("Optimizers", package = "gips")` or its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Optimizers.html) -
#'     A place to learn more about
#'     the available optimizers.
#' * `vignette("Theory", package = "gips")` or its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html) -
#'     A place to learn more about
#'     the math behind the `gips` package.
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 5
#' mu <- runif(perm_size, -10, 10) # Assume we don't know the mean
#' sigma_matrix <- matrix(
#'   data = c(
#'     1.0, 0.8, 0.6, 0.6, 0.8,
#'     0.8, 1.0, 0.8, 0.6, 0.6,
#'     0.6, 0.8, 1.0, 0.8, 0.6,
#'     0.6, 0.6, 0.8, 1.0, 0.8,
#'     0.8, 0.6, 0.6, 0.8, 1.0
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' ) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5)
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- cov(Z) # Assume we have to estimate the mean
#'
#' g <- gips(S, number_of_observations)
#'
#' g_map <- find_MAP(g, max_iter = 5, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
#' g_map
#'
#' g_map2 <- find_MAP(g_map, max_iter = 5, show_progress_bar = FALSE, optimizer = "continue")
#'
#' if (require("graphics")) {
#'   plot(g_map2, type = "both", logarithmic_x = TRUE)
#' }
#'
#' g_map_BF <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
#' summary(g_map_BF)
find_MAP <- function(g, max_iter = NA, optimizer = NA,
                     show_progress_bar = TRUE,
                     save_all_perms = FALSE,
                     return_probabilities = FALSE) {
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
    if (optimizer %in% c("BF", "brute_force", "full")) {
      rlang::abort(c("There was a problem identified with provided arguments:",
        "i" = "`optimizer == 'continue'` cannot be provided after optimizating with `optimizer == 'brute_force'`, because the whole space was already browsed.",
        "x" = "You provided `optimizer == 'continue'`, but the gips object `g` was optimized with brute_force optimizer. Better permutation will not be found."
      ))
    }
  }

  if (!(optimizer %in% c("MH", "Metropolis_Hastings", "BF", "brute_force", "full")) && return_probabilities) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "Probabilities can only be returned with the `optimizer == 'Metropolis_Hastings'` or `optimizer == 'brute_force'`",
      "x" = "You provided both `!(optimizer %in% c('Metropolis_Hastings', 'brute_force'))` and `return_probabilities == TRUE`!",
      "i" = "Did You want to use `optimizer == 'Metropolis_Hastings'` or `optimizer == 'brute_force'`, or `return_probabilities == FLASE`?"
    ))
  }

  if (return_probabilities && optimizer %in% c("MH", "Metropolis_Hastings")) {
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
    rlang::inform(c(
      paste0(
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

  if (was_mean_estimated) { # one degree of freedom is lost; we will return this 1 to number_of_observations after optimization in `combine_gips()`
    edited_number_of_observations <- number_of_observations - 1
  } else {
    edited_number_of_observations <- number_of_observations
  }

  start_time <- Sys.time()

  if (optimizer %in% c("MH", "Metropolis_Hastings")) {
    gips_optimized <- Metropolis_Hastings_optimizer(
      S = S, number_of_observations = edited_number_of_observations,
      max_iter = max_iter, start_perm = start_perm,
      delta = delta, D_matrix = D_matrix,
      return_probabilities = return_probabilities,
      save_all_perms = save_all_perms,
      show_progress_bar = show_progress_bar
    )
  } else if (optimizer %in% c("HC", "hill_climbing")) {
    gips_optimized <- hill_climbing_optimizer(
      S = S, number_of_observations = edited_number_of_observations,
      max_iter = max_iter, start_perm = start_perm,
      delta = delta, D_matrix = D_matrix,
      save_all_perms = save_all_perms,
      show_progress_bar = show_progress_bar
    )
  } else if (optimizer %in% c("BF", "brute_force", "full")) {
    gips_optimized <- brute_force_optimizer(
      S = S, number_of_observations = edited_number_of_observations,
      delta = delta, D_matrix = D_matrix,
      return_probabilities = return_probabilities,
      save_all_perms = save_all_perms,
      show_progress_bar = show_progress_bar
    )
  }

  end_time <- Sys.time()
  attr(gips_optimized, "optimization_info")[["optimization_time"]] <- end_time - start_time
  attr(gips_optimized, "optimization_info")[["whole_optimization_time"]] <- end_time - start_time

  structure_constants <- get_structure_constants(gips_optimized[[1]])
  n0 <- max(structure_constants[["r"]] * structure_constants[["d"]] / structure_constants[["k"]])
  if (attr(g, "was_mean_estimated")) { # correction for estimating the mean
    n0 <- n0 + 1
  }
  if (n0 > number_of_observations) {
    rlang::warn(c(
      paste0(
        "The found permutation has n0 = ", n0,
        " which is bigger than the number_of_observations = ",
        number_of_observations, "."
      ),
      "i" = "The covariance matrix invariant under the found permutation does not have the likelihood properly defined.",
      "i" = "For a more in-depth explanation, see the 'Project Matrix - Equation (6)' section in `vignette('Theory', package = 'gips')` or its pkgdown page: https://przechoj.github.io/gips/articles/Theory.html."
    ))
  }


  return(combine_gips(g, gips_optimized))
}


Metropolis_Hastings_optimizer <- function(S,
    number_of_observations, max_iter, start_perm = NULL,
    delta = 3, D_matrix = NULL, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = TRUE) {
  if (is.null(start_perm)) {
    start_perm <- permutations::id
  }

  check_correctness_of_arguments(
    S = S, number_of_observations = number_of_observations,
    max_iter = max_iter, start_perm = start_perm,
    delta = delta, D_matrix = D_matrix, was_mean_estimated = FALSE,
    return_probabilities = return_probabilities,
    save_all_perms = save_all_perms,
    show_progress_bar = show_progress_bar
  )

  if (!inherits(start_perm, "gips_perm")) {
    start_perm <- gips_perm(start_perm, nrow(S)) # now we know the `S` is a matrix
  }

  if (is.infinite(max_iter)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "`max_iter` in `Metropolis_Hastings_optimizer` must be finite.",
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

  my_goal_function <- function(perm, i) {
    out_val <- log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `number_of_observations` if the mean was estimated!
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix
    )

    if (is.nan(out_val) || is.infinite(out_val)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::abort(c(
        "gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(out_val), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrix. If it is not the case for You, please get in touch with us on ISSUE#5.",
        "x" = paste0("The Metropolis Hastings algorithm was stopped after ", i, " iterations.")
      ))
    }

    out_val
  }

  acceptance <- rep(FALSE, max_iter)
  log_posteriori_values <- rep(0, max_iter)
  if (save_all_perms) {
    visited_perms <- list()
    visited_perms[[1]] <- start_perm
  } else {
    visited_perms <- NA
  }
  current_perm <- start_perm

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  }
  log_posteriori_values[1] <- my_goal_function(current_perm, 0)

  found_perm <- start_perm
  found_perm_log_posteriori <- log_posteriori_values[1]

  Uniformly_drawn_numbers <- stats::runif(max_iter, min = 0, max = 1)

  # main loop
  for (i in 1:(max_iter - 1)) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }

    e <- runif_transposition(perm_size)
    perm_proposal <- compose_with_transposition(current_perm, e)

    goal_function_perm_proposal <- my_goal_function(perm_proposal, i)

    # if goal_function_perm_proposal > log_posteriori_values[i], then it is true, because Uniformly_drawn_numbers[i] \in [0,1]
    if (Uniformly_drawn_numbers[i] < exp(goal_function_perm_proposal - log_posteriori_values[i])) { # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2). That means this is Metropolis algorithm, not necessary Metropolis-Hastings.
      current_perm <- perm_proposal
      if (save_all_perms) {
        visited_perms[[i + 1]] <- current_perm
      }
      log_posteriori_values[i + 1] <- goal_function_perm_proposal
      acceptance[i] <- TRUE

      if (found_perm_log_posteriori < log_posteriori_values[i + 1]) {
        found_perm_log_posteriori <- log_posteriori_values[i + 1]
        found_perm <- current_perm
      }
    } else {
      if (save_all_perms) {
        visited_perms[[i + 1]] <- current_perm
      }
      log_posteriori_values[i + 1] <- log_posteriori_values[i] # TODO(Do we really want to forget the calculated values? The algorithm HC works differently)
    }
  }

  if (show_progress_bar) {
    close(progressBar)
  }

  function_calls <- length(log_posteriori_values)

  # visited_perms are already either a list of things or a `NULL` object

  if (return_probabilities) {
    probabilities <- estimate_probabilities(visited_perms, show_progress_bar)
  } else {
    probabilities <- NULL
  }

  optimization_info <- list(
    "acceptance_rate" = mean(acceptance),
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = visited_perms,
    "start_perm" = start_perm,
    "last_perm" = current_perm,
    "last_perm_log_posteriori" = log_posteriori_values[function_calls],
    "iterations_performed" = i,
    "optimization_algorithm_used" = "Metropolis_Hastings",
    "post_probabilities" = probabilities,
    "did_converge" = NULL,
    "best_perm_log_posteriori" = found_perm_log_posteriori,
    "optimization_time" = NA,
    "whole_optimization_time" = NA
  )


  new_gips(
    list(found_perm), S, number_of_observations,
    delta, D_matrix,
    was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}


hill_climbing_optimizer <- function(S,
    number_of_observations, max_iter = 5,
    start_perm = NULL, delta = 3, D_matrix = NULL,
    save_all_perms = FALSE, show_progress_bar = TRUE) {
  if (is.null(start_perm)) {
    start_perm <- permutations::id
  }

  check_correctness_of_arguments(
    S = S, number_of_observations = number_of_observations,
    max_iter = max_iter, start_perm = start_perm,
    delta = delta, D_matrix = D_matrix, was_mean_estimated = FALSE,
    return_probabilities = FALSE, save_all_perms = save_all_perms,
    show_progress_bar = show_progress_bar
  )

  if (!inherits(start_perm, "gips_perm")) {
    start_perm <- gips_perm(start_perm, nrow(S)) # now we know the `S` is a matrix
  }

  if (show_progress_bar && is.infinite(max_iter)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "x" = "You tried to run `find_MAP(show_progress_bar=TRUE, max_iter=Inf)`.",
      "i" = "Progress bar is not yet supported for infinite max_iter.",
      "i" = "Do You want to use `show_progress_bar=FALSE` or a finite `max_iter`?",
      "i" = "For more information on progress bar see ISSUE#8."
    ))
  }

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  }

  perm_size <- dim(S)[1]

  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size)
  }


  my_goal_function <- function(perm, i) {
    out_val <- log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `number_of_observations` if the mean was estimated!
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix
    )

    if (is.nan(out_val) || is.infinite(out_val)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::abort(c(
        "gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(out_val), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrix. If it is not the case for You, please get in touch with us on ISSUE#5.",
        "x" = paste0("The Hill Climbing algorithm was stopped after ", i, " iterations.")
      ))
    }

    out_val
  }

  goal_function_best_logvalues <- numeric(0)
  log_posteriori_values <- numeric(0)

  # init
  if (save_all_perms) {
    visited_perms <- list()
    visited_perms[[1]] <- start_perm
  } else {
    visited_perms <- NA
  }
  current_perm <- start_perm

  goal_function_best_logvalues[1] <- my_goal_function(current_perm, 0)
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
        neighbour <- compose_with_transposition(current_perm, c(i, j))
        neighbour_value <- my_goal_function(neighbour, iteration)
        log_posteriori_values[length(log_posteriori_values) + 1] <- neighbour_value

        if (neighbour_value > best_neighbour_value) {
          best_neighbour_value <- neighbour_value
          best_neighbour <- neighbour
        }
      }
    }

    if (best_neighbour_value > goal_function_best_logvalues[iteration]) {
      goal_function_best_logvalues[iteration + 1] <- best_neighbour_value
      current_perm <- best_neighbour
      if (save_all_perms) {
        visited_perms[[iteration + 1]] <- best_neighbour
      }
    } else {
      did_converge <- TRUE
      break
    }
  }

  last_perm <- current_perm

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
    "visited_perms" = visited_perms,
    "start_perm" = start_perm,
    "last_perm" = last_perm,
    "last_perm_log_posteriori" = goal_function_best_logvalues[iteration],
    "iterations_performed" = iteration,
    "optimization_algorithm_used" = "hill_climbing",
    "post_probabilities" = NULL,
    "did_converge" = did_converge,
    "best_perm_log_posteriori" = goal_function_best_logvalues[iteration],
    "optimization_time" = NA,
    "whole_optimization_time" = NA
  )


  new_gips(
    list(last_perm), S, number_of_observations,
    delta, D_matrix,
    was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}


brute_force_optimizer <- function(
    S,
    number_of_observations,
    delta = 3, D_matrix = NULL,
    return_probabilities = return_probabilities,
    save_all_perms = FALSE, show_progress_bar = TRUE) {
  check_correctness_of_arguments(
    S = S, number_of_observations = number_of_observations,
    max_iter = 5, start_perm = permutations::id, # max_iter, was_mean_estimated and start_perm are not important for optimization with brute_force
    delta = delta, D_matrix = D_matrix, was_mean_estimated = FALSE,
    return_probabilities = return_probabilities, save_all_perms = save_all_perms,
    show_progress_bar = show_progress_bar
  )

  perm_size <- dim(S)[1]

  if (perm_size > 18) {
    rlang::abort(c("Optimizer 'brute_force' cannot browse such a big permutional space.",
      "x" = paste0(
        "You provided a space with size ", perm_size,
        "! (factorial), which has ", prod(1:perm_size),
        " elements."
      ),
      "i" = "Do You want to use other optimizer for such a big space? For example 'Metropolis_Hastings' or 'hill_climbing'?"
    ))
  }

  if (perm_size > 9) { # I don't know how to test this without running the optimization...
    rlang::warn(c("Optimizer 'brute_force' will take very long time to browse such a big permutional space.",
      "x" = paste0(
        "You provided a space with size ", perm_size,
        "! (factorial), which has ", prod(1:perm_size),
        " elements."
      ),
      "i" = "Do You want to use other optimizer for such a big space? For example 'Metropolis_Hastings' or 'hill_climbing'?"
    ))
  }

  iterations_to_perform <-
    if ((3 <= perm_size) && (perm_size <= 9)) {
      # Only the generators are interesting for us:
      # perm_group_generators are calculated only for up to perm_size = 9
      # See ISSUE#21 for more information
      OEIS_A051625[perm_size]
    } else {
      prod(1:perm_size)
    }

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = iterations_to_perform, initial = 1)
  }

  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size)
  }

  my_goal_function <- function(perm, i) {
    out_val <- log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `number_of_observations` if the mean was estimated!
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix
    )

    if (is.nan(out_val) || is.infinite(out_val)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::abort(c(
        "gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(out_val), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrix. If it is not the case for You, please get in touch with us on ISSUE#5.",
        "x" = paste0("The Brute Force algorithm was stopped after ", i, " iterations.")
      ))
    }

    out_val
  }

  # main loop
  all_perms_list <- permutations::allperms(perm_size)
  all_perms_list <- permutations::as.cycle(all_perms_list)
  if ((3 <= perm_size) && (perm_size <= 9)) {
    # Only the generators are interesting for us:
    # perm_group_generators are calculated only for up to perm_size = 9
    # See ISSUE#21 for more information
    all_perms_list <- all_perms_list[perm_group_generators_list[[perm_size - 2]]]
  }
  log_posteriori_values <- sapply(1:length(all_perms_list), function(i) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }
    this_perm <- permutations::cycle(list(all_perms_list[[i]]))
    my_goal_function(this_perm, i)
  })

  if (show_progress_bar) {
    close(progressBar)
  }

  if (return_probabilities) { # calculate exact probabilities
    probabilities <- calculate_probabilities(all_perms_list, log_posteriori_values, show_progress_bar)
  } else {
    probabilities <- NULL
  }

  best_perm <- gips_perm(permutations::cycle(list(all_perms_list[[which.max(log_posteriori_values)]])), perm_size)

  if (save_all_perms) {
    visited_perms <- all_perms_list
  } else {
    visited_perms <- NA
  }

  optimization_info <- list(
    "acceptance_rate" = NULL,
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = visited_perms,
    "start_perm" = permutations::id,
    "last_perm" = NULL,
    "last_perm_log_posteriori" = NULL,
    "iterations_performed" = iterations_to_perform,
    "optimization_algorithm_used" = "brute_force",
    "post_probabilities" = probabilities,
    "did_converge" = TRUE,
    "best_perm_log_posteriori" = log_posteriori_values[which.max(log_posteriori_values)],
    "optimization_time" = NA,
    "whole_optimization_time" = NA
  )


  new_gips(
    list(best_perm), S, number_of_observations,
    delta, D_matrix,
    was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}



#' Combining 2 gips objects
#'
#' g2 was optimized with a single optimization method. g1 was potentially non-optimized or optimized once, or optimized multiple times.
#' If g2 was optimized with "brute_force", forget the g1.
#'
#' @noRd
combine_gips <- function(g1, g2, show_progress_bar = FALSE) {
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

  if (all(is.na(optimization_info1[["visited_perms"]])) || all(is.na(optimization_info2[["visited_perms"]]))) {
    if (!all(is.na(optimization_info1[["visited_perms"]])) || !all(is.na(optimization_info2[["visited_perms"]]))) {
      rlang::warn("You wanted to save visited_perms on one of the optimized `gips` objects but forget it for the other. This is not possible, so both will be forgotten.")
      optimization_info2[["post_probabilities"]] <- NULL
      optimization_info1[["post_probabilities"]] <- NULL
    }
    visited_perms <- NA
  } else {
    visited_perms <- c(optimization_info1[["visited_perms"]], optimization_info2[["visited_perms"]]) # WoW, one can use `c()` to combine lists!
  }
  optimization_algorithm_used <- c(optimization_info1[["optimization_algorithm_used"]], optimization_info2[["optimization_algorithm_used"]])

  if (all(optimization_algorithm_used == "Metropolis_Hastings") &&
    !is.null(optimization_info2[["post_probabilities"]])) {
    post_probabilities <- estimate_probabilities(visited_perms, show_progress_bar) # TODO(This can be combined more optimally when !is.null(optimization_info1[["post_probabilities"]]). It is significant, because those calculations are like the same speed as the MH itself. However, I (Adam) think this will be rarely done nevertheless.)
  } else {
    post_probabilities <- NULL
  }

  optimization_info_new <- list(
    "acceptance_rate" = (n1 * optimization_info1[["acceptance_rate"]] + n2 * optimization_info2[["acceptance_rate"]]) / (n1 + n2),
    "log_posteriori_values" = c(optimization_info1[["log_posteriori_values"]], optimization_info2[["log_posteriori_values"]]),
    "visited_perms" = visited_perms,
    "start_perm" = optimization_info1[["start_perm"]],
    "last_perm" = optimization_info2[["last_perm"]],
    "last_perm_log_posteriori" = optimization_info2[["last_perm_log_posteriori"]],
    "iterations_performed" = c(optimization_info1[["iterations_performed"]], optimization_info2[["iterations_performed"]]),
    "optimization_algorithm_used" = optimization_algorithm_used,
    "post_probabilities" = post_probabilities,
    "did_converge" = optimization_info2[["did_converge"]],
    "best_perm_log_posteriori" = max(optimization_info1[["best_perm_log_posteriori"]], optimization_info2[["best_perm_log_posteriori"]]),
    "optimization_time" = c(optimization_info1[["optimization_time"]], optimization_info2[["optimization_time"]]),
    "whole_optimization_time" = optimization_info1[["whole_optimization_time"]] + optimization_info2[["whole_optimization_time"]]
  )

  if (optimization_info1[["best_perm_log_posteriori"]] > optimization_info2[["best_perm_log_posteriori"]]) {
    g_out <- g1 # in the continuation the new best was not found
  } else {
    g_out <- g2
  }

  attr(g_out, "optimization_info") <- optimization_info_new

  g_out
}
