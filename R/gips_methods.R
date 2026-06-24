#' Printing `gips` object
#'
#' Printing function for a `gips` class.
#'
#' @param x An object of a `gips` class.
#' @param digits The number of decimal places for
#'     the posterior probability. It can be negative.
#'     By default, `Inf`. It is passed to [base::round()].
#' @param compare_to_original A logical. Whether to print how many
#'     times more likely the current permutation is compared to:
#' * the identity permutation `()` (for unoptimized `gips` object);
#' * the starting permutation (for optimized `gips` object).
#' @param log_value A logical. Whether to print the logarithmic value.
#'     Default to `FALSE`.
#' @param oneline A logical. Whether to print in
#'     one or multiple lines. Default to `FALSE`.
#' @param ... The additional arguments passed to [base::cat()].
#'
#' @seealso
#' * [find_MAP()] - The function that makes
#'     an optimized `gips` object out of the unoptimized one.
#' * [compare_posteriories_of_perms()] - The function that prints
#'     the compared posteriories between any two permutations,
#'     not only compared to the starting one or id.
#'
#' @returns An invisible `NULL`.
#' @export
#'
#' @examples
#' S <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
#' g <- gips(S, 10, perm = "(12)")
#' print(g, digits = 4, oneline = TRUE)
print.gips <- function(x, digits = 3, compare_to_original = TRUE,
                       log_value = FALSE, oneline = FALSE, ...) {
  validate_gips(x)

  printing_text <- paste0("The permutation ", as.character(x[[1]]))

  if (is.null(attr(x, "optimization_info"))) { # it is unoptimized gips object
    log_posteriori <- log_posteriori_of_gips(x)
    if (is.nan(log_posteriori) || is.infinite(log_posteriori)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(log_posteriori), "NaN", "Inf"), " occurred!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrix. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))
    }

    if (compare_to_original) {
      x_id <- gips(
        S = attr(x, "S"),
        number_of_observations = attr(x, "number_of_observations"),
        delta = attr(x, "delta"), D_matrix = attr(x, "D_matrix"),
        was_mean_estimated = attr(x, "was_mean_estimated"), perm = ""
      )
      log_posteriori_id <- log_posteriori_of_gips(x_id)

      printing_text <- c(
        printing_text,
        paste0(
          "is ", convert_log_diff_to_str(log_posteriori - log_posteriori_id, digits),
          " times more likely than the () permutation"
        )
      )
    }
  } else { # it is optimized gips object
    log_posteriori <- attr(x, "optimization_info")[["best_perm_log_posteriori"]]
    x_original <- gips(
      S = attr(x, "S"),
      number_of_observations = attr(x, "number_of_observations"),
      delta = attr(x, "delta"), D_matrix = attr(x, "D_matrix"),
      was_mean_estimated = attr(x, "was_mean_estimated"), perm = attr(x, "optimization_info")[["original_perm"]]
    )
    log_posteriori_original <- log_posteriori_of_gips(x_original)
    
    if (is.nan(log_posteriori) || is.infinite(log_posteriori)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(log_posteriori), "NaN", "Inf"), " occurred!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrix. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))
    }

    printing_text <- c(printing_text, paste0(
      "was found after ",
      length(attr(x, "optimization_info")[["log_posteriori_values"]]),
      " posteriori calculations"
    ))

    if (compare_to_original) {
      printing_text <- c(printing_text, paste0(
        "is ", convert_log_diff_to_str(log_posteriori - log_posteriori_original, digits),
        " times more likely than the ",
        as.character(x_original), " permutation"
      ))
    }
  }

  if (log_value) {
    printing_text <- c(
      printing_text,
      paste0(
        "has log posteriori ",
        round(log_posteriori, digits = digits)
      )
    )
  }

  # The first line will end with ":", all following lines will end with ";".
  cat(
    paste0(c(
      printing_text[1],
      paste0(printing_text[-1],
        collapse = ifelse(oneline, "; ", ";\n - ")
      )
    ), collapse = ifelse(oneline, ": ", ":\n - ")),
    ".\n",
    sep = "", ...
  )

  invisible(NULL)
}


#' Convert the log difference to the appropriate string.
#' If bigger than 10 millions, use the scientific notation
#'
#' @param digits Number of digits after comma
#'
#' @examples
#' convert_log_diff_to_str(1009.5, 3) == "2.632e+438"
#' convert_log_diff_to_str(16.1, 3) == "9820670.922"
#' convert_log_diff_to_str(16.2, 3) == "1.085e+7"
#' convert_log_diff_to_str(-7.677, 3) == "4.634e-4"
#'
#' @noRd
convert_log_diff_to_str <- function(log_diff, digits) {
  if (is.infinite(log_diff)) {
    return(ifelse(log_diff > 0, "Inf", "-Inf"))
  }
  if (log_diff == 0) {
    return("1")
  }

  times_more_likely <- round(
    exp(log_diff),
    digits = digits
  )

  log10_diff <- log_diff * log10(exp(1))

  ifelse(0 < times_more_likely && times_more_likely < 10000000,
    as.character(times_more_likely),
    paste0(
      round(10^(log10_diff - floor(log10_diff)), digits = digits),
      "e",
      ifelse(log_diff > 0, "+", ""), # If log_diff < 0, then floor(log10_diff) will have "-" in front
      floor(log10_diff)
    )
  )
}

# Based on `stats::summary.lm()`
#' Summarizing the gips object
#'
#' `summary` method for `gips` class.
#'
#' @param object An object of class `gips`. Usually, a result of a [find_MAP()].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function `summary.gips()` computes and returns a list of summary
#'     statistics of the given `gips` object. Those are:
#' * For unoptimized `gips` object:
#'   1. `optimized` - `FALSE`.
#'   2. `start_permutation` - the permutation this `gips` represents.
#'   3. `start_permutation_log_posteriori` - the log of the A Posteriori
#'       value the start permutation has.
#'   4. `times_more_likely_than_id` - how many times more likely
#'       the `start_permutation` is over the identity permutation, `()`.
#'       It can be less than 1, meaning the identity permutation
#'       is more likely. Remember that this number can become
#'       very large and overflow to `Inf` or small and underflow to 0.
#'   5. `log_times_more_likely_than_id` - log of `times_more_likely_than_id`.
#'   6. `likelihood_ratio_test_statistics`, `likelihood_ratio_test_p_value` - 
#'       statistics and p-value of Likelihood Ratio test, where
#'       the H_0 is that the data was drawn from the normal distribution
#'       with Covariance matrix invariant under the given permutation.
#'       The p-value is calculated from the asymptotic distribution.
#'       Note that this is sensibly defined only for \eqn{n \ge p}.
#'   7. `n0` - the minimum number of observations needed for
#'       the covariance matrix's maximum likelihood estimator
#'       (corresponding to a MAP) to exist. See **\eqn{C\sigma} and `n0`**
#'       section in `vignette("Theory", package = "gips")` or in its
#'       [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'   8. `S_matrix` - the underlying matrix.
#'       This matrix will be used in calculations of
#'       the posteriori value in [log_posteriori_of_gips()].
#'   9. `number_of_observations` - the number of observations that
#'       were observed for the `S_matrix` to be calculated.
#'       This value will be used in calculations of
#'       the posteriori value in [log_posteriori_of_gips()].
#'   10. `was_mean_estimated` - given by the user while creating the `gips` object:
#'       * `TRUE` means the `S` parameter was the output of [stats::cov()] function;
#'       * `FALSE` means the `S` parameter was calculated with
#'           `S = t(X) %*% X / number_of_observations`.
#'   11. `delta`, `D_matrix` - the hyperparameters of the Bayesian method.
#'       See the **Hyperparameters** section of [gips()] documentation.
#'   12. `n_parameters` - number of free parameters in the covariance matrix.
#'   13. `AIC`, `BIC` - output of [AIC.gips()] and [BIC.gips()] functions.
#' * For optimized `gips` object:
#'   1. `optimized` - `TRUE`.
#'   2. `found_permutation` - the permutation this `gips` represents.
#'       The visited permutation with the biggest A Posteriori value.
#'   3. `found_permutation_log_posteriori` - the log of the A Posteriori
#'       value the found permutation has.
#'   4. `start_permutation` - the original permutation this `gips`
#'       represented before optimization. It is the first visited permutation.
#'   5. `start_permutation_log_posteriori` - the log of the A Posteriori
#'       value the start permutation has.
#'   6. `times_more_likely_than_start` - how many more likely
#'       the `found_permutation` is over the `start_permutation`.
#'       It cannot be a number less than 1.
#'       Remember that this number can big and overflow to `Inf`.
#'   7. `log_times_more_likely_than_start` - log of
#'       `times_more_likely_than_start`.
#'   8. `likelihood_ratio_test_statistics`, `likelihood_ratio_test_p_value` - 
#'       statistics and p-value of Likelihood Ratio test, where
#'       the H_0 is that the data was drawn from the normal distribution
#'       with Covariance matrix invariant under `found_permutation`.
#'       The p-value is calculated from the asymptotic distribution.
#'       Note that this is sensibly defined only for \eqn{n \ge p}.
#'   9. `n0` - the minimal number of observations needed for the existence of
#'       the maximum likelihood estimator (corresponding to a MAP) of
#'       the covariance matrix (see **\eqn{C\sigma} and `n0`**
#'       section in `vignette("Theory", package = "gips")` or in its
#'       [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html)).
#'   10. `S_matrix` - the underlying matrix.
#'       This matrix will be used in calculations of
#'       the posteriori value in [log_posteriori_of_gips()].
#'   11. `number_of_observations` - the number of observations that
#'       were observed for the `S_matrix` to be calculated.
#'       This value will be used in calculations of
#'       the posteriori value in [log_posteriori_of_gips()].
#'   12. `was_mean_estimated` - given by the user while creating the `gips` object:
#'       * `TRUE` means the `S` parameter was output of the [stats::cov()] function;
#'       * `FALSE` means the `S` parameter was calculated with
#'           `S = t(X) %*% X / number_of_observations`.
#'   13. `delta`, `D_matrix` - the hyperparameters of the Bayesian method.
#'       See the **Hyperparameters** section of [gips()] documentation.
#'   14. `n_parameters` - number of free parameters in the covariance matrix.
#'   15. `AIC`, `BIC` - output of [AIC.gips()] and [BIC.gips()] functions.
#'   16. `optimization_algorithm_used` - all used optimization algorithms
#'       in order (one could start optimization with "MH", and then
#'       do an "HC").
#'   17. `did_converge` - a boolean, did the last used algorithm converge.
#'   18. `number_of_log_posteriori_calls` - how many times was
#'       the [log_posteriori_of_gips()] function called during
#'       the optimization.
#'   19. `whole_optimization_time` - how long was the optimization process;
#'       the sum of all optimization times (when there were multiple).
#'   20. `log_posteriori_calls_after_best` - how many times was
#'       the [log_posteriori_of_gips()] function called after
#'       the `found_permutation`; in other words, how long ago
#'       could the optimization be stopped and have the same result.
#'       If this value is small, consider running [find_MAP()]
#'       again with `optimizer = "continue"`.
#'       For `optimizer = "BF"`, it is `NULL`.
#'   21. `acceptance_rate` - only interesting for `optimizer = "MH"`.
#'       How often was the algorithm accepting the change of permutation
#'       in an iteration.
#' @export
#' 
#' @importFrom stats pchisq
#'
#' @seealso
#' * [find_MAP()] - Usually, the `summary.gips()`
#'     is called on the output of `find_MAP()`.
#' * [log_posteriori_of_gips()] - Calculate
#'     the likelihood of a permutation.
#' * [AIC.gips()], [BIC.gips()] - Calculate
#'     Akaike's or Bayesian Information Criterion
#' * [project_matrix()] - Project the known
#'     matrix on the found permutations space.
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 6
#' mu <- runif(6, -10, 10) # Assume we don't know the mean
#' sigma_matrix <- matrix(
#'   data = c(
#'     1.1, 0.8, 0.6, 0.4, 0.6, 0.8,
#'     0.8, 1.1, 0.8, 0.6, 0.4, 0.6,
#'     0.6, 0.8, 1.1, 0.8, 0.6, 0.4,
#'     0.4, 0.6, 0.8, 1.1, 0.8, 0.6,
#'     0.6, 0.4, 0.6, 0.8, 1.1, 0.8,
#'     0.8, 0.6, 0.4, 0.6, 0.8, 1.1
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' ) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- cov(Z) # Assume we have to estimate the mean
#'
#' g <- gips(S, number_of_observations)
#' unclass(summary(g))
#'
#' g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
#' unclass(summary(g_map))
#'
#' g_map2 <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "hill_climbing")
#' summary(g_map2)
summary.gips <- function(object, ...) {
  # validate_gips(object) # validation is done in `log_posteriori_of_gips()`
  permutation_log_posteriori <- log_posteriori_of_gips(object)

  tmp <- get_n0_and_edited_number_of_observations_from_gips(object)
  n0 <- tmp[1]
  edited_number_of_observations <- tmp[2]
  
  n_parameters <- sum(get_structure_constants(object[[1]])[["dim_omega"]])
  
  # Likelihood-Ratio test:
  if (edited_number_of_observations < n0 || !is.positive.definite.matrix(attr(object, "S"))) {
    likelihood_ratio_test_statistics <- NULL
    likelihood_ratio_test_p_value <- NULL
  } else {
    likelihood_ratio_test_statistics <- edited_number_of_observations*(determinant(project_matrix(attr(object, "S"), object[[1]]))$modulus - determinant(attr(object, "S"))$modulus)
    attributes(likelihood_ratio_test_statistics) <- NULL
    p <- attr(object[[1]], "size")
    df_chisq <- p*(p+1)/2 - n_parameters
    if (df_chisq == 0) {
      likelihood_ratio_test_p_value <- NULL
    } else {
      # when likelihood_ratio_test_statistics is close to 0, the H_0
      likelihood_ratio_test_p_value <- 1 - pchisq(likelihood_ratio_test_statistics, df_chisq)
    }
  }

  if (is.null(attr(object, "optimization_info"))) {
    log_posteriori_id <- log_posteriori_of_perm(
      perm_proposal = "", S = attr(object, "S"),
      number_of_observations = edited_number_of_observations,
      delta = attr(object, "delta"), D_matrix = attr(object, "D_matrix")
    )

    summary_list <- list(
      optimized = FALSE,
      start_permutation = object[[1]],
      start_permutation_log_posteriori = permutation_log_posteriori,
      times_more_likely_than_id = exp(permutation_log_posteriori - log_posteriori_id),
      log_times_more_likely_than_id = permutation_log_posteriori - log_posteriori_id,
      likelihood_ratio_test_statistics = likelihood_ratio_test_statistics,
      likelihood_ratio_test_p_value = likelihood_ratio_test_p_value,
      n0 = n0,
      S_matrix = attr(object, "S"),
      number_of_observations = attr(object, "number_of_observations"),
      was_mean_estimated = attr(object, "was_mean_estimated"),
      delta = attr(object, "delta"),
      D_matrix = attr(object, "D_matrix"),
      n_parameters = n_parameters,
      AIC = suppressWarnings(AIC(object, classes = c("singular_matrix", "likelihood_does_not_exists"))), # warning for NA and NULL
      BIC = suppressWarnings(BIC(object, classes = c("singular_matrix", "likelihood_does_not_exists"))) # warning for NA and NULL
    )
  } else {
    optimization_info <- attr(object, "optimization_info")

    if (optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force") {
      when_was_best <- which(abs(optimization_info[["log_posteriori_values"]] - permutation_log_posteriori) < 0.0000001) # close enough; this is the first generator of the group
      log_posteriori_calls_after_best <- length(optimization_info[["log_posteriori_values"]]) - when_was_best[1]
      start_permutation <- optimization_info[["start_perm"]]
      start_permutation_log_posteriori <- optimization_info[["log_posteriori_values"]][1]
    } else {
      # for brute_force when_was_best is useless.
      # Also, the `optimization_info[["visited_perms"]]` is a list, but
        # its elements are not of class `gips_perm`, because it was done with
        # `optimization_info[["visited_perms"]] <- permutations::allperms()`
      # Real original permutation is saved in optimization_info[["original_perm"]]
      when_was_best <- NULL
      log_posteriori_calls_after_best <- NULL
      start_permutation <- optimization_info[["original_perm"]]
      gips_start <- gips(
        S = attr(object, "S"),
        number_of_observations = attr(object, "number_of_observations"),
        delta = attr(object, "delta"),
        D_matrix = attr(object, "D_matrix"),
        was_mean_estimated = attr(object, "was_mean_estimated"),
        perm = start_permutation
      )
      start_permutation_log_posteriori <- log_posteriori_of_gips(gips_start)
    }
    
    summary_list <- list(
      optimized = TRUE,
      found_permutation = object[[1]],
      found_permutation_log_posteriori = permutation_log_posteriori,
      start_permutation = start_permutation,
      start_permutation_log_posteriori = start_permutation_log_posteriori,
      times_more_likely_than_start = exp(permutation_log_posteriori - start_permutation_log_posteriori),
      log_times_more_likely_than_start = permutation_log_posteriori - start_permutation_log_posteriori,
      likelihood_ratio_test_statistics = likelihood_ratio_test_statistics,
      likelihood_ratio_test_p_value = likelihood_ratio_test_p_value,
      n0 = n0,
      S_matrix = attr(object, "S"),
      number_of_observations = attr(object, "number_of_observations"),
      was_mean_estimated = attr(object, "was_mean_estimated"),
      delta = attr(object, "delta"),
      D_matrix = attr(object, "D_matrix"),
      n_parameters = n_parameters,
      AIC = suppressWarnings(AIC(object), classes = c("singular_matrix", "likelihood_does_not_exists")), # warning for NA and NULL
      BIC = suppressWarnings(BIC(object), classes = c("singular_matrix", "likelihood_does_not_exists")), # warning for NA and NULL
      optimization_algorithm_used = optimization_info[["optimization_algorithm_used"]],
      did_converge = optimization_info[["did_converge"]],
      number_of_log_posteriori_calls = length(optimization_info[["log_posteriori_values"]]),
      whole_optimization_time = optimization_info[["whole_optimization_time"]],
      log_posteriori_calls_after_best = log_posteriori_calls_after_best,
      acceptance_rate = optimization_info[["acceptance_rate"]]
    )
  }

  structure(summary_list,
    class = c("summary.gips")
  )
}

# Based on `sloop::s3_get_method(print.summary.lm)`
#' @method print summary.gips
#' @param x An object of class `summary.gips` to be printed
#' @describeIn summary.gips Printing method for class `summary.gips`.
#'     Prints every interesting information in a form pleasant for humans.
#' @returns The function `print.summary.gips()` returns an invisible `NULL`.
#' @export
#'
#' @examples
#' # ================================================================================
#' S <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
#' g <- gips(S, 10)
#' print(summary(g))
print.summary.gips <- function(x, ...) {
  cat(
    ifelse(x[["optimized"]],
      "The optimized",
      "The unoptimized"
    ),
    " `gips` object.\n\nPermutation:\n ",
    ifelse(x[["optimized"]],
      as.character(x[["found_permutation"]]),
      as.character(x[["start_permutation"]])
    ),
    "\n\nLog_posteriori:\n ",
    ifelse(x[["optimized"]],
      x[["found_permutation_log_posteriori"]],
      x[["start_permutation_log_posteriori"]]
    ),
    ifelse(!x[["optimized"]] && as.character(x[["start_permutation"]]) == "()",
      "",
      paste0(
        "\n\nTimes more likely than ",
        ifelse(x[["optimized"]],
          paste0(
            "starting permutation:\n ",
            convert_log_diff_to_str(x[["log_times_more_likely_than_start"]], 3)
          ),
          paste0(
            "identity permutation:\n ",
            convert_log_diff_to_str(x[["log_times_more_likely_than_id"]], 3)
          )
        )
      )
    ),
    ifelse(is.null(x[["likelihood_ratio_test_statistics"]]),
      ifelse(is.positive.definite.matrix(x[["S_matrix"]]),
        "\n\ndet(S) == 0, so Likelihood-Ratio test cannot be performed",
        "\n\nn0 > number_of_observations, so Likelihood-Ratio test cannot be performed"
      ),
      ifelse(is.null(x[["likelihood_ratio_test_p_value"]]),
        "\n\nThe current permutation is id, so Likelihood-Ratio test cannot be performed (there is nothing to compare)",
        paste0("\n\nThe p-value of Likelihood-Ratio test:\n ", format(x[["likelihood_ratio_test_p_value"]], digits = 4)))
      ),
    "\n\nThe number of observations:\n ", x[["number_of_observations"]],
    "\n\n", ifelse(x[["was_mean_estimated"]],
      paste0(
        "The mean in the `S` matrix was estimated.\nTherefore, one degree of freedom was lost.\nThere are ",
        x[["number_of_observations"]] - 1, " degrees of freedom left."
      ),
      paste0(
        "The mean in the `S` matrix was not estimated.\nTherefore, all degrees of freedom were preserved (",
        x[["number_of_observations"]], ")."
      )
    ),
    "\n\nn0:\n ", x[["n0"]],
    "\n\nThe number of observations is ",
    ifelse(x[["n0"]] > x[["number_of_observations"]],
      "smaller",
      ifelse(x[["n0"]] == x[["number_of_observations"]],
        "equal",
        "bigger"
      )
    ),
    " than n0 for this permutation,\nso the gips model based on the found permutation does ",
    ifelse(x[["n0"]] > x[["number_of_observations"]],
      "not ", ""
    ), "exist.",
    "\n\nThe number of free parameters in the covariance matrix:\n ", x[["n_parameters"]],
    "\n\nBIC:\n ", ifelse(x[["n0"]] > x[["number_of_observations"]],
      "The number of observations is smaller than n0 for this permutation,\n so the gips model based on the found permutation does not exist.", x[["BIC"]]
    ),
    "\n\nAIC:\n ", ifelse(x[["n0"]] > x[["number_of_observations"]],
      "The number of observations is smaller than n0 for this permutation,\n so the gips model based on the found permutation does not exist.", x[["AIC"]]
    ),
    sep = ""
  )

  if (x[["optimized"]]) {
    cat("\n\n", paste0(rep("-", getOption("width")), collapse = ""), # the line
      sep = ""
    )

    if (length(x[["optimization_algorithm_used"]]) == 1) {
      cat("\nOptimization algorithm:\n ", x[["optimization_algorithm_used"]],
        ifelse(x[["optimization_algorithm_used"]] == "hill_climbing",
          ifelse(x[["did_converge"]],
            " did converge",
            " did not converge"
          ),
          ""
        ),
        sep = ""
      )
    } else {
      # multiple optimizations
      cat("\nOptimization algorithms:\n ", paste0(x[["optimization_algorithm_used"]], collapse = ", "),
        ifelse(x[["optimization_algorithm_used"]][length(x[["optimization_algorithm_used"]])] == "hill_climbing",
          ifelse(x[["did_converge"]],
            "\n The last hill_climbing did converge.",
            "\n The last hill_climbing did not converge."
          ),
          ""
        ),
        sep = ""
      )
    }

    cat("\n\nThe number of log_posteriori calls:\n ",
      x[["number_of_log_posteriori_calls"]],
      "\n\nOptimization time:\n ", unclass(x[["whole_optimization_time"]]),
      " ", attr(x[["whole_optimization_time"]], "units"),
      sep = ""
    )

    if (all(x[["optimization_algorithm_used"]] == "Metropolis_Hastings")) {
      cat(paste0("\n\nAcceptance rate:\n ", x[["acceptance_rate"]]), sep = "")
    }
    if (x[["optimization_algorithm_used"]][length(x[["optimization_algorithm_used"]])] != "brute_force") {
      cat("\n\nLog_posteriori calls after the found permutation:\n ",
        x[["log_posteriori_calls_after_best"]],
        sep = ""
      )
    }
  }

  cat("\n")

  invisible(NULL)
}

#' Internal
#' @return a vector of length 2 with n0 and edited_number_of_observations
#' @noRd
get_n0_and_edited_number_of_observations_from_gips <- function(g) {
  validate_gips(g)
  
  n0 <- get_n0_from_perm(g[[1]], attr(g, "was_mean_estimated"))
  
  edited_number_of_observations <- attr(g, "number_of_observations")
  
  if (attr(g, "was_mean_estimated")) { # correction for estimating the mean
    edited_number_of_observations <- edited_number_of_observations - 1
  }
  
  c(n0, edited_number_of_observations)
}

#' Internal
#' @return (integer) n0
#' @noRd
get_n0_from_perm <- function(g_perm, was_mean_estimated) {
  structure_constants <- get_structure_constants(g_perm)
  n0 <- max(structure_constants[["r"]] * structure_constants[["d"]] / structure_constants[["k"]])

  if (was_mean_estimated) { # correction for estimating the mean
    n0 <- n0 + 1
  }

  c(n0)
}

#' Extract the Log-Likelihood for `gips` class
#'
#' Calculates Log-Likelihood of the sample based on the `gips` object.
#'
#' This will always be largest for `perm = "()"` (provided that `p <= n`).
#'
#' If the found permutation still requires more parameters than `n`,
#'     the likelihood does not exist; thus the function returns `NULL`.
#'
#' If the `projected_cov` (output of [project_matrix()])
#'     is close to singular, the `NA` is returned.
#'
#' @param object An object of class `gips`. Usually, a result of a [find_MAP()].
#' @param ... Further arguments will be ignored.
#'
#' @section Existence of likelihood:
#' We only consider the non-degenerate multivariate normal model.
#' In the `gips` context, such a model exists only when
#' the number of observations is bigger or equal to `n0`. To get `n0`
#' for the `gips` object `g`, call `summary(g)$n0`.
#' 
#' See examples where the `g_n_too_small` had too small
#' `number_of_observations` to have likelihood. After the optimization,
#' the likelihood did exist.
#'
#' For more information, refer to **\eqn{C_\sigma} and `n0`** section in
#' `vignette("Theory", package = "gips")` or its
#' [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'
#' @section Calculation details:
#' For more details and used formulas, see
#' the **Information Criterion - AIC and BIC** section in
#' `vignette("Theory", package = "gips")` or its
#' [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'
#' @importFrom stats logLik
#'
#' @returns Log-Likelihood of the sample. Object of class `logLik`.
#'
#' Possible failure situations:
#' * When the multivariate normal model does not exist
#'     (`number_of_observations < n0`), it returns `NULL`.
#' * When the multivariate normal model cannot be reasonably approximated
#'     (output of [project_matrix()] is singular), it returns `-Inf`.
#'
#' In both failure situations, it shows a warning.
#' More information can be found in the **Existence of likelihood** section below.
#'
#' @export
#'
#' @seealso
#' * [logLik()] - Generic function this [logLik.gips()] extends.
#' * [find_MAP()] - Usually, the `logLik.gips()`
#'     is called on the output of `find_MAP()`.
#' * [AIC.gips()], [BIC.gips()] - Often, one is more
#'     interested in an Information Criterion AIC or BIC.
#' * [summary.gips()] - One can get `n0` by calling `summary(g)$n0`.
#'     To see why one may be interested in `n0`,
#'     see the **Existence of likelihood** section above.
#' * [project_matrix()] - Project the known matrix
#'     onto the found permutations space.
#'     It is mentioned in the **Calculation details** section above.
#'
#' @examples
#' S <- matrix(c(
#'   5.15, 2.05, 3.60, 1.99,
#'   2.05, 5.09, 2.03, 3.57,
#'   3.60, 2.03, 5.21, 1.97,
#'   1.99, 3.57, 1.97, 5.13
#' ), nrow = 4)
#' g <- gips(S, 5)
#' logLik(g) # -32.67048
#' # For perm = "()", which is default, there is p + choose(p, 2) degrees of freedom
#'
#' g_map <- find_MAP(g, optimizer = "brute_force")
#' logLik(g_map) # -32.6722 # this will always be smaller than `logLik(gips(S, n, perm = ""))`
#'
#' g_n_too_small <- gips(S, number_of_observations = 4)
#' logLik(g_n_too_small) # NULL # the likelihood does not exists
#' summary(g_n_too_small)$n0 # 5, but we set number_of_observations = 4, which is smaller
#' 
#' g_MAP <- find_MAP(g_n_too_small)
#' logLik(g_MAP) # -24.94048, this is no longer NULL
#' summary(g_MAP)$n0 # 2
logLik.gips <- function(object, ...) {
  validate_gips(object)

  original_cov <- attributes(object)[["S"]]
  projected_cov <- project_matrix(original_cov, object[[1]])
  p <- ncol(original_cov)
  n <- attr(object, "number_of_observations")

  tmp <- get_n0_and_edited_number_of_observations_from_gips(object)
  n0 <- tmp[1]
  edited_number_of_observations <- tmp[2]

  if (n < n0) { # Likelihood is not defined in that setting
    rlang::warn(c(
      "The likelihood is not defined for this `gips`.",
      "x" = paste0(
        "The n = ", n,
        " is smaller than the minimum required n0 = ", n0,
        ". For more information, see section **Existence of likelihood** in documentation `?logLik.gips` or its [pkgdown page](https://przechoj.github.io/gips/reference/logLik.gips.html)."
      )
    ), class = "likelihood_does_not_exists")

    return(NULL)
  }

  log_det_projected_cov <- determinant(projected_cov, logarithm = TRUE)[["modulus"]]
  attributes(log_det_projected_cov) <- NULL
  if (is.infinite(log_det_projected_cov)) {
    rlang::warn(c(
      "The projected matrix is computationally singular.",
      "x" = "The likelihood for singular matrices cannot be estimated with a satisfying precision.",
      "i" = paste0("Reciprocal condition number = ", rcond(projected_cov), ".")
    ), class = "singular_matrix")

    return(-Inf)
  }

  log_2pi_plus_1 <- 2.837877066409345483560659472811235279722794947275566825634303 # log(2*pi) + 1

  log_L_S <- -edited_number_of_observations * (p * log_2pi_plus_1 + log_det_projected_cov) / 2

  n_parameters <- sum(get_structure_constants(object[[1]])[["dim_omega"]])

  attr(log_L_S, "df") <- n_parameters
  attr(log_L_S, "nobs") <- n # The AIC and BIC will use n, not edited_number_of_observations
  
  class(log_L_S) <- "logLik"
  
  log_L_S
}

#' Akaike's An Information Criterion for `gips` class
#'
#' @section Calculation details:
#' For more details and used formulas, see
#' the **Information Criterion - AIC and BIC** section in
#' `vignette("Theory", package = "gips")` or its
#' [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'
#' @method AIC gips
#'
#' @param object An object of class `gips`. Usually, a result of a [find_MAP()].
#' @param ... Further arguments will be ignored.
#' @param k Numeric, the *penalty* per parameter to be used.
#'     The default `k = 2` is the classical AIC.
#'
#' @returns `AIC.gips()` returns calculated Akaike's An Information Criterion
#'
#' When the multivariate normal model does not exist
#'     (`number_of_observations < n0`), it returns `NULL`.
#' When the multivariate normal model cannot be reasonably approximated
#'     (output of [project_matrix()] is singular), it returns `Inf`.
#'
#' In both failure situations, shows a warning.
#' More information can be found in the **Existence of likelihood**
#' section of [logLik.gips()].
#'
#' @importFrom stats AIC
#'
#' @seealso
#' * [AIC()], [BIC()] - Generic functions
#'     this `AIC.gips()` and `BIC.gips()` extend.
#' * [find_MAP()] - Usually, the `AIC.gips()` and `BIC.gips()`
#'     are called on the output of `find_MAP()`.
#' * [logLik.gips()] - Calculates the log-likelihood for
#'     the `gips` object. An important part of the Information Criteria.
#'
#' @export
#' @examples
#' S <- matrix(c(
#'   5.15, 2.05, 3.10, 1.99,
#'   2.05, 5.09, 2.03, 3.07,
#'   3.10, 2.03, 5.21, 1.97,
#'   1.99, 3.07, 1.97, 5.13
#' ), nrow = 4)
#' g <- gips(S, 14)
#' g_map <- find_MAP(g, optimizer = "brute_force")
#'
#' AIC(g) # 238
#' AIC(g_map) # 224 < 238, so g_map is better than g according to AIC
AIC.gips <- function(object, ..., k = 2) {
  log_likelihood_S <- logLik.gips(object) # in here we will validate object is of class gips

  if (is.null(log_likelihood_S)) {
    return(NULL)
  }

  if (is.infinite(log_likelihood_S)) {
    return(Inf)
  }

  -2 * as.numeric(log_likelihood_S) + k * attr(log_likelihood_S, "df")
}

#' @method BIC gips
#' @describeIn AIC.gips Schwarz's Bayesian Information Criterion
#'
#' @importFrom stats BIC
#'
#' @returns `BIC.gips()` returns calculated
#'     Schwarz's Bayesian Information Criterion.
#'
#' @export
#' @examples
#' # ================================================================================
#' BIC(g) # 244
#' BIC(g_map) # 226 < 244, so g_map is better than g according to BIC
BIC.gips <- function(object, ...) {
  log_likelihood_S <- logLik.gips(object) # in here we will validate object is of class gips

  if (is.null(log_likelihood_S)) {
    return(NULL)
  }

  if (is.infinite(log_likelihood_S)) {
    return(Inf)
  }

  k <- log(attr(log_likelihood_S, "nobs")) # this line is the only difference from `AIC.gips()`
  -2 * as.numeric(log_likelihood_S) + k * attr(log_likelihood_S, "df")
}

#' Extract probabilities for `gips` object optimized with `return_probabilities = TRUE`
#'
#' After the `gips` object was optimized with
#' the `find_MAP(return_probabilities = TRUE)` function, then
#' those calculated probabilities can be extracted with this function.
#'
#' @param g An object of class `gips`.
#'     A result of a `find_MAP(return_probabilities = TRUE)`.
#'
#' @returns A numeric vector of calculated probability values.
#' Names contain the permutations that these probabilities represent.
#' For `gips` object optimized with `find_MAP(return_probabilities = FALSE)`,
#' it returns a `NULL` object.
#' It is sorted according to the probability.
#'
#' @export
#'
#' @seealso
#' * [find_MAP()] - The `get_probabilities_from_gips()`
#'     is called on the output of
#'     `find_MAP(return_probabilities = TRUE, save_all_perms = TRUE)`.
#' * `vignette("Optimizers", package = "gips")` or its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Optimizers.html)) -
#'     A place to learn more about the available optimizers.
#'
#' @examples
#' g <- gips(matrix(c(1, 0.5, 0.5, 1.3), nrow = 2), 13, was_mean_estimated = FALSE)
#' g_map <- find_MAP(g,
#'   optimizer = "BF", show_progress_bar = FALSE,
#'   return_probabilities = TRUE, save_all_perms = TRUE
#' )
#'
#' get_probabilities_from_gips(g_map)
get_probabilities_from_gips <- function(g) {
  validate_gips(g)

  if (is.null(attr(g, "optimization_info"))) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`gips` object has to be optimized with `find_MAP(return_probabilities=TRUE)` to use `get_probabilities_from_gips()` function.",
      "x" = "You did not optimize `g`.",
      "i" = "Did You use the wrong `g` as an argument for this function?",
      "i" = "Did You forget to optimize `g`?"
    ))
  }

  if (is.null(attr(g, "optimization_info")[["post_probabilities"]])) {
    rlang::inform(c(
      "You called `get_probabilities_from_gips(g)` on the `gips` object that does not have saved probabilities.",
      "x" = "`NULL` will be returned",
      "i" = "Did You use the wrong `g` as an argument for this function?",
      "i" = "Did You forget to optimize with `g <- find_MAP(return_probabilities = TRUE)`?",
      "i" = "Did You unintentionally use `g <- forget_perms(g)`?"
    ))
  }

  attr(g, "optimization_info")[["post_probabilities"]]
}


#' Forget the permutations for `gips` object optimized with `save_all_perms = TRUE`
#'
#' Slim the `gips` object by forgetting the visited permutations from `find_MAP(save_all_perms = TRUE)`.
#'
#' For example, `perm_size = 150` and `max_iter = 150000` we checked `forget_perms()` saves ~350 MB of RAM.
#'
#' @param g An object of class `gips`.
#'     A result of a `find_MAP(save_all_perms = TRUE)`.
#'
#' @returns The same object `g` as given,
#'     but without the visited permutation list.
#'
#' @export
#'
#' @seealso
#' * [find_MAP()] - The `forget_perms()` is called on
#'     the output of `find_MAP(save_all_perms = TRUE)`.
#'
#' @examples
#' A <- matrix(rnorm(10 * 10), nrow = 10)
#' S <- t(A) %*% A
#' g <- gips(S, 13, was_mean_estimated = FALSE)
#' g_map <- find_MAP(g,
#'   max_iter = 10, optimizer = "Metropolis_Hastings",
#'   show_progress_bar = FALSE, save_all_perms = TRUE
#' )
#'
#' object.size(g_map) # ~18 KB
#' g_map_slim <- forget_perms(g_map)
#' object.size(g_map_slim) # ~8 KB
forget_perms <- function(g) {
  validate_gips(g)

  optimization_info <- attr(g, "optimization_info")

  if (is.null(optimization_info)) {
    rlang::inform(c(
      "Provided `g` is a `gips` object, but it was not optimized yet.",
      "i" = "Did You provide the wrong `gips` object?"
    ))
  } else if (all(is.na(optimization_info[["visited_perms"]]))) {
    rlang::inform(c(
      "Provided `g` is an optimized `gips` object that already has forgotten all permutations.",
      "i" = "Did You provide the wrong `gips` object?"
    ))
  } else {
    optimization_info[["visited_perms"]] <- I(NA)
    attr(g, "optimization_info") <- optimization_info
  }

  g
}


#' Transform the `gips` object to a character vector
#'
#' Implementation of the S3 method.
#'
#' @inheritParams print.gips
#' @param ... Further arguments (currently ignored).
#'
#' @method as.character gips
#'
#' @returns An object of a `character` type.
#'
#' @seealso
#' * [as.character.gips_perm()] - The underlying `gips_perm` of
#'     the `gips` object is passed to [as.character.gips_perm()].
#' * [permutations::as.character.cycle()] - The underlying permutation of
#'     the `gips` object is passed to [permutations::as.character.cycle()].
#'
#' @export
#'
#' @examples
#' A <- matrix(rnorm(4 * 4), nrow = 4)
#' S <- t(A) %*% A
#' g <- gips(S, 14, perm = "(123)")
#' as.character(g)
as.character.gips <- function(x, ...) {
  validate_gips(x)

  as.character(x[[1]], ...)
}
