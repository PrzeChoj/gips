#' The constructor of a `gips` class.
#'
#' Create a `gips` object.
#' This object will contain initial data and all other information
#' needed to find the most likely invariant permutation.
#' It will not perform optimization; one must call
#' the [find_MAP()] function to do it. See the examples below.
#'
#' @param S A matrix; empirical covariance matrix.
#'     When `Z` is the observed data:
#' * if one does not know the theoretical mean and has to
#'     estimate it with the observed mean, use `S = cov(Z)`,
#'     and leave parameter `was_mean_estimated = TRUE` as default;
#' * if one knows the theoretical mean is 0, use
#'     `S = (t(Z) %*% Z) / number_of_observations`, and set
#'     parameter `was_mean_estimated = FALSE`.
#' @param number_of_observations A number of data points
#'     that `S` is based on.
#' @param delta A number, hyper-parameter of a Bayesian model.
#'     It has to be strictly bigger than 1.
#'     See the **Hyperparameters** section below.
#' @param D_matrix Symmetric, positive-definite matrix of the same size as `S`.
#'     Hyper-parameter of a Bayesian model.
#'     When `NULL`, the (hopefully) reasonable one is derived from the data.
#'     For more details, see the **Hyperparameters** section below.
#' @param was_mean_estimated A boolean.
#' * Set `TRUE` (default) when your `S` parameter is a result of
#'     a [stats::cov()] function.
#' * Set `FALSE` when your `S` parameter is a result of
#'     a `(t(Z) %*% Z) / number_of_observations` calculation.
#' @param perm An optional permutation to be the base for the `gips` object.
#'     It can be of a `gips_perm` or a `permutation` class, or anything
#'     the function [permutations::permutation()] can handle.
#'     It can also be of a `gips` class, but
#'     it will be interpreted as the underlying `gips_perm`.
#'
#' @section Methods for a `gips` class:
#' * [summary.gips()]
#' * [plot.gips()]
#' * [print.gips()]
#' * [logLik.gips()]
#' * [AIC.gips()]
#' * [BIC.gips()]
#' * [as.character.gips()]
#'
#' @section Hyperparameters:
#' We encourage you to try `D_matrix = d * I`, where `I` is a `p` \eqn{\times}
#'     `p` identity matrix and `d > 0` for some different `d`.
#' When `d` is small compared to the data (e.g., `d=0.1 * mean(diag(S))`),
#'     bigger structures will be found.
#' When `d` is big compared to the data (e.g., `d=100 * mean(diag(S))`),
#'     the posterior distribution does not depend on the data.
#'
#' Taking `D_matrix = d * I` is equivalent to setting `S <- S / d`.
#'
#' The default for `D_matrix` is `D_matrix = d * I`, where
#' `d = mean(diag(S))`, which is equivalent to modifying `S`
#' so that the mean value on the diagonal is 1.
#'
#' In the Bayesian model, the prior distribution for
#' the covariance matrix is a generalized case of
#' [Wishart distribution](https://en.wikipedia.org/wiki/Wishart_distribution).
#'
#' For a brief introduction, see the **Bayesian model selection**
#' section in `vignette("Theory", package = "gips")` or in its
#' [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html)).
#'
#' For analysis of the Hyperparameters influence, see **Section 3.2.**
#' of "Learning permutation symmetries with gips in R"
#' by `gips` developers Adam Chojecki, Paweł Morgen, and Bartosz Kołodziejek,
#' [Journal of Statistical Software](https://doi.org/10.18637/jss.v112.i07);
#' \doi{10.18637/jss.v112.i07}.
#'
#' @section Multi-sample:
#' `gips()` supports G groups of samples sharing the same permutation symmetry.
#' Pass a list of G covariance matrices as `S` and a numeric vector of length G
#' as `number_of_observations`:
#'
#' ```r
#' g <- gips(list(S1, S2, S3), c(n1, n2, n3))
#' ```
#'
#' `D_matrix` should then be a list of G positive-definite matrices
#' (defaulting to `diag(mean(diag(S_g)), p)` for each group).
#' `delta` can be a scalar (broadcast to all groups) or a vector of length G
#' (one value per group); it defaults to `rep(3, G)`.
#'
#' @returns `gips()` returns an object of
#'     a `gips` class after the safety checks.
#'
#' @export
#' @seealso
#' * [stats::cov()] - The `S` parameter, as an empirical covariance matrix,
#'     is most of the time a result of the `cov()` function.
#'     For more information, see
#'     [Wikipedia - Estimation of covariance matrices](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices).
#' * [find_MAP()] - The function that finds
#'     the Maximum A Posteriori (MAP) Estimator
#'     for a given `gips` object.
#' * [gips_perm()] - The constructor of a `gips_perm` class.
#'     The `gips_perm` object is used as the base object for
#'     the `gips` object. To be more precise, the base object
#'     for `gips` is a one-element list of a `gips_perm` object.
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 5
#' mu <- runif(5, -10, 10) # Assume we don't know the mean
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
#' g_map <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
#' g_map
#'
#' summary(g_map)
#'
#' if (require("graphics")) {
#'   plot(g_map, type = "both", logarithmic_x = TRUE)
#' }
gips <- function(S, number_of_observations, delta = 3, D_matrix = NULL,
                 was_mean_estimated = TRUE, perm = "") {
  if (inherits(perm, "gips")) {
    validate_gips(perm)
    perm <- perm[[1]]
  }
  if (!inherits(perm, c("gips_perm", "permutation"))) {
    perm <- permutations::permutation(perm)
  }

  if (is.list(S)) {
    # Multi-sample path
    G <- length(S)
    p <- ncol(S[[1]])

    # Normalize delta: a scalar is broadcast to a per-group vector of length G
    if (length(delta) == 1) delta <- rep(delta, G)

    check_gips_mult_arguments(
      Ss = S, numbers_of_observations = number_of_observations,
      delta = delta, D_matrices = D_matrix,
      was_mean_estimated = was_mean_estimated, perm = perm
    )

    if (is.null(D_matrix)) {
      D_matrix <- lapply(S, function(y) diag(x = mean(diag(y)), nrow = ncol(y)))
    }

    gips_perm_object <- if (inherits(perm, "gips_perm")) perm else gips_perm(perm, p)

    return(validate_gips(new_gips(
      list(gips_perm_object), S, number_of_observations,
      delta = delta, D_matrix = D_matrix,
      was_mean_estimated = was_mean_estimated, optimization_info = NULL
    )))
  }

  check_gips_arguments(
    S = S, number_of_observations = number_of_observations,
    delta = delta, D_matrix = D_matrix, was_mean_estimated = was_mean_estimated,
    perm = perm
  )

  if (inherits(perm, "gips_perm")) {
    gips_perm_object <- perm # it is already a `gips_perm`
  } else {
    gips_perm_object <- gips_perm(perm, nrow(S)) # it is of a `cycle` class from permutations package (it was checked in `check_gips_arguments`). Make it 'gips_perm' class
  }


  if (is.null(D_matrix)) {
    D_matrix <- diag(x = mean(diag(S)), nrow = ncol(S))
  }

  validate_gips(new_gips(
    list(gips_perm_object), S, number_of_observations,
    delta = delta, D_matrix = D_matrix,
    was_mean_estimated = was_mean_estimated, optimization_info = NULL
  ))
}


#' @describeIn gips Constructor. It is only intended for low-level use.
#'
#' @param list_of_gips_perm A list with a single element of
#'     a `gips_perm` class. The base object for the `gips` object.
#' @param optimization_info For internal use only. `NULL` or the list with
#'     information about the optimization process.
#'
#' @returns `new_gips()` returns an object of
#'     a `gips` class without safety checks.
#'
#' @export
new_gips <- function(list_of_gips_perm, S, number_of_observations,
                     delta, D_matrix, was_mean_estimated, optimization_info) {
  valid_S <- is.matrix(S) || (is.list(S) && list_of_matrices_check(S))
  valid_n <- is.numeric(number_of_observations) && all(sapply(number_of_observations, is.wholenumber))
  valid_D <- is.matrix(D_matrix) || (is.list(D_matrix) && list_of_matrices_check(D_matrix))

  if (!is.list(list_of_gips_perm) ||
    !inherits(list_of_gips_perm[[1]], "gips_perm") ||
    !valid_S ||
    !valid_n ||
    !is.numeric(delta) ||
    !valid_D ||
    !is.logical(was_mean_estimated) ||
    !(is.null(optimization_info) || is.list(optimization_info))) {
    rlang::abort(c("x" = "`gips` object cannot be created from those arguments."))
  }


  structure(list_of_gips_perm,
    S = S, number_of_observations = number_of_observations,
    delta = delta, D_matrix = D_matrix, was_mean_estimated = was_mean_estimated,
    optimization_info = optimization_info,
    class = c("gips")
  )
}


#' @describeIn gips Validator. It is only intended for low-level use.
#'
#' @param g Object to be checked whether it is a proper object of a `gips` class.
#'
#' @returns `validate_gips()` returns its argument unchanged.
#'     If the argument is not a proper element of a `gips` class,
#'     it produces an error.
#'
#' @export
validate_gips <- function(g) {
  if (!(inherits(g, "gips"))) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "`g` must be of a `gips` class.",
      "x" = paste0(
        "You provided `g` with `class(g) == (",
        paste(class(g), collapse = ", "), ")`."
      )
    ))
  }

  if (!(length(g) == 1)) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `length(g)` must be `1`.",
      "x" = paste0(
        "You provided `g` with `length(g) == ",
        length(g), "`."
      )
    ))
  }
  if (!is.list(g)) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `g` must be a list.",
      "x" = paste0(
        "You provided `g` with `typeof(g) == '",
        typeof(g), "'."
      )
    ))
  }

  perm <- g[[1]]
  S <- attr(g, "S")
  number_of_observations <- attr(g, "number_of_observations")
  delta <- attr(g, "delta")
  D_matrix <- attr(g, "D_matrix")
  was_mean_estimated <- attr(g, "was_mean_estimated")
  optimization_info <- attr(g, "optimization_info")

  if (!inherits(perm, "gips_perm")) {
    rlang::abort(c(
      "There was a problem identified with provided argument:",
      "i" = "The `g[[1]]` must be an object of class `gips_perm`.",
      "x" = paste0(
        "You provided `g[[1]]` with `class(g[[1]]) == (",
        paste(class(perm), collapse = ", "),
        ")`."
      )
    ))
  }
  
  tryCatch(
    {
      validate_gips_perm(perm)
    },
    error = function(cond) {
      rlang::abort(c(
        "There was a problem identified with provided argument:",
        "i" = "The `g[[1]]` must be a valid object of class `gips_perm`.",
        "x" = paste0(
          "You provided `g[[1]]` with class `gips_perm`, ",
          "but it does not pass `validate_gips_perm(g[[1]])`."
        )
      ))
    }
  )

  if (is.list(S)) {
    check_gips_mult_arguments(
      Ss = S, numbers_of_observations = number_of_observations,
      delta = delta, D_matrices = D_matrix,
      was_mean_estimated = was_mean_estimated, perm = perm
    )
  } else {
    check_gips_arguments(
      S = S, number_of_observations = number_of_observations,
      delta = delta, D_matrix = D_matrix, was_mean_estimated = was_mean_estimated,
      perm = perm
    )
  }

  if (!(is.null(optimization_info) || is.list(optimization_info))) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `optimization_info` value must be either a `NULL`, or a list.",
      "x" = paste0(
        "You provided `attr(g, 'optimization_info')` with type ",
        typeof(optimization_info), "."
      )
    ))
  }

  g
}


# Helper for validating logical flag arguments
check_logical_flag <- function(value, name) {
  error_msg <- character(0)
  if (!is.logical(value)) {
    error_msg <- c(
      "i" = paste0("`", name, "` must be a logic value (`TRUE` or `FALSE`)."),
      "x" = paste0(
        "You provided `", name, "` with type ",
        typeof(value), "."
      )
    )
  } else if (is.na(value)) {
    error_msg <- c(
      "i" = paste0("`", name, "` must be a logic value (`TRUE` or `FALSE`)."),
      "x" = paste0("You provided `", name, "` as an `NA`.")
    )
  }
  error_msg
}


#' Validate arguments for gips constructor and validator
#'
#' Internal function to validate arguments passed to gips() and validate_gips().
#'
#' @param S The empirical covariance matrix
#' @param number_of_observations Number of observations
#' @param delta Bayesian hyperparameter
#' @param D_matrix Bayesian hyperparameter matrix
#' @param was_mean_estimated Whether the mean was estimated
#' @param perm Starting or provided permutation
#'
#' @returns Invisibly. Throws an error if validation fails.
#'
#' @keywords internal
#' 
#' @noRd
check_gips_arguments <- function(S, number_of_observations, delta, D_matrix,
                                 was_mean_estimated, perm) {
  if (!is.matrix(S)) {
    rlang::abort(c("There was a problem identified with provided S argument:",
      "i" = "`S` must be a matrix.",
      "x" = paste0(
        "You provided `S` with type ",
        typeof(S), "."
      )
    ))
  }
  abort_text <- character(0)
  additional_info <- 0
  
  if (ncol(S) != nrow(S)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a square matrix.",
      "x" = paste0(
        "You provided `S` as a matrix, but with different sizes: ",
        ncol(S), " and ", nrow(S), "."
      )
    )
  } else if (!is.numeric(S)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a numeric matrix.",
      "x" = paste0(
        "You provided `S` as a matrix, but with non-numeric values. Your provided type ",
        typeof(S), "."
      )
    )
  } else if (!all(abs(S - t(S)) < 0.000001)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a symmetric matrix.",
      "x" = "You provided `S` as a matrix, but a non-symmetric one.",
      "i" = "Is your matrix approximately symmetric? Maybe try setting `S <- (S+t(S))/2`?"
    )
    additional_info <- additional_info + 1
  } else if (!is.positive.semi.definite.matrix(S, tolerance = 1e-06)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be positive semi-definite matrix.",
      "x" = "You provided `S` as a symmetric matrix, but a non-positive-semi-definite one."
    )
  }
  
  if (is.null(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must not be `NULL`.",
      "x" = "Your provided `number_of_observations` is `NULL`."
    )
  } else if (number_of_observations < 1) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must be at least 1.",
      "x" = paste0(
        "You provided `number_of_observations == ",
        number_of_observations, "`."
      )
    )
  } else if (!is.wholenumber(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must be a whole number.",
      "x" = paste0(
        "You provided `number_of_observations == ",
        number_of_observations, "`."
      )
    )
  }
  
  if (is.null(delta)) {
    abort_text <- c(abort_text,
      "i" = "`delta` must not be `NULL`.",
      "x" = "Your provided `delta` is a `NULL`."
    )
  } else if (delta <= 1) {
    abort_text <- c(abort_text,
      "i" = "`delta` must be strictly bigger than 1.",
      "x" = paste0("You provided `delta == ", delta, "`.")
    )
  }
  
  if (!(is.null(D_matrix) || is.matrix(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must either be `NULL` or a matrix.",
      "x" = paste0(
        "You provided `D_matrix` with type ",
        typeof(D_matrix), "."
      )
    )
  } else if (!(is.null(D_matrix) || ncol(D_matrix) == nrow(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must either be `NULL` or a square matrix.",
      "x" = paste0(
        "You provided `D_matrix` as a matrix, but with different sizes: ",
        ncol(D_matrix), " and ", nrow(D_matrix), "."
      )
    )
  } else if (!(is.null(D_matrix) || ncol(S) == ncol(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`S` must be a square matrix with the same shape as a square matrix `D_matrix`.",
      "x" = paste0(
        "You provided `S` with shape ",
        ncol(S), " and ", nrow(S),
        ", but also `D_matrix` with shape ",
        ncol(D_matrix), " and ", nrow(D_matrix), "."
      )
    )
  } else if (any(is.nan(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must not contain any `NaN`s.",
      "x" = "You provided `D_matrix` with `NaN`s!"
    )
  } else if (any(is.infinite(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must not contain any infinite values.",
      "x" = "You provided `D_matrix` with infinite values!"
    )
  }
  
  abort_text <- c(abort_text, check_logical_flag(was_mean_estimated, "was_mean_estimated"))
  
  if (!(permutations::is.cycle(perm) || inherits(perm, "gips_perm"))) {
    abort_text <- c(abort_text,
      "i" = "`perm` must be the output of `gips_perm()` function, or of a `cycle` class from `permutations` package.",
      "x" = paste0(
        "You provided `perm` with `class(perm) == (",
        paste(class(perm), collapse = ", "),
        ")`."
      )
    )
  } else if (!(permutations::is.cycle(perm) || attr(perm, "size") == ncol(S))) {
    abort_text <- c(abort_text,
      "i" = "`perm` must have the `size` attribute equal to the shape of a square matrix `S`",
      "x" = paste0(
        "You provided `perm` with `size == ",
        attr(perm, "size"),
        "`, but the `S` matrix You provided has ",
        ncol(S), " columns."
      )
    )
  }
  
  if (length(abort_text) > 0) {
    abort_text <- c(
      paste0(
        "There were ", (length(abort_text) - additional_info) / 2,
        " problems identified with the provided arguments:"
      ),
      abort_text
    )

    if (length(abort_text) > 11) {
      abort_text <- c(
        abort_text[1:11],
        paste0("... and ", (length(abort_text) - 1) / 2 - 5, " more problems")
      )
    }

    abort_text <- c(
      abort_text,
      ">" = "If You think You've found a bug in a package, please open an ISSUE on https://github.com/PrzeChoj/gips/issues"
    )

    rlang::abort(abort_text)
  }
  
  invisible(NULL)
}


#' Validate arguments for find_MAP and optimization functions
#'
#' Internal function to validate arguments passed to optimization functions.
#'
#' @param S The empirical covariance matrix
#' @param number_of_observations Number of observations
#' @param delta Bayesian hyperparameter
#' @param D_matrix Bayesian hyperparameter matrix
#' @param was_mean_estimated Whether the mean was estimated
#' @param max_iter Maximum number of iterations
#' @param start_perm Starting permutation
#' @param return_probabilities Whether to return probabilities
#' @param save_all_perms Whether to save all permutations
#' @param show_progress_bar Whether to show a progress bar
#'
#' @returns Invisibly. Throws an error if validation fails.
#'
#' @keywords internal
#' 
#' @noRd
check_find_MAP_arguments <- function(S, number_of_observations, max_iter, start_perm,
                                     delta, D_matrix, was_mean_estimated,
                                     return_probabilities, save_all_perms, show_progress_bar) {
  if (!is.matrix(S) && !is.list(S)) {
    rlang::abort(c("There was a problem identified with provided S argument:",
      "i" = "`S` must be a matrix.",
      "x" = paste0(
        "You provided `S` with type ",
        typeof(S), "."
      )
    ))
  }

  # For multi-sample (list S), the detailed validation was done in gips().
  if (is.list(S)) {
    abort_text <- character(0)
    abort_text <- c(abort_text, check_logical_flag(was_mean_estimated, "was_mean_estimated"))
    abort_text <- c(abort_text, check_logical_flag(return_probabilities, "return_probabilities"))
    abort_text <- c(abort_text, check_logical_flag(save_all_perms, "save_all_perms"))
    abort_text <- c(abort_text, check_logical_flag(show_progress_bar, "show_progress_bar"))
    if (length(abort_text) > 0) {
      rlang::abort(c("There were problems identified with the provided arguments:", abort_text))
    }
    if (return_probabilities && !save_all_perms) {
      rlang::abort(c("There was a problem identified with the provided arguments:",
        "i" = "For calculations of probabilities, all perms have to be available after the optimization process.",
        "x" = "You provided `return_probabilities == TRUE` and `save_all_perms == FALSE`!",
        "i" = "Did You want to set `save_all_perms = TRUE`?"
      ))
    }
    return(invisible(NULL))
  }

  abort_text <- character(0)
  additional_info <- 0
  
  if (ncol(S) != nrow(S)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a square matrix.",
      "x" = paste0(
        "You provided `S` as a matrix, but with different sizes: ",
        ncol(S), " and ", nrow(S), "."
      )
    )
  } else if (!is.numeric(S)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a numeric matrix.",
      "x" = paste0(
        "You provided `S` as a matrix, but with non-numeric values. Your provided type ",
        typeof(S), "."
      )
    )
  } else if (!all(abs(S - t(S)) < 0.000001)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a symmetric matrix.",
      "x" = "You provided `S` as a matrix, but a non-symmetric one.",
      "i" = "Is your matrix approximately symmetric? Maybe try setting `S <- (S+t(S))/2`?"
    )
    additional_info <- additional_info + 1
  } else if (!is.positive.semi.definite.matrix(S, tolerance = 1e-06)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be positive semi-definite matrix.",
      "x" = "You provided `S` as a symmetric matrix, but a non-positive-semi-definite one."
    )
  }
  
  if (is.null(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must not be `NULL`.",
      "x" = "Your provided `number_of_observations` is `NULL`."
    )
  } else if (number_of_observations < 1) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must be at least 1.",
      "x" = paste0(
        "You provided `number_of_observations == ",
        number_of_observations, "`."
      )
    )
  } else if (!is.wholenumber(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must be a whole number.",
      "x" = paste0(
        "You provided `number_of_observations == ",
        number_of_observations, "`."
      )
    )
  }
  
  if (!(is.infinite(max_iter) || is.wholenumber(max_iter))) {
    abort_text <- c(abort_text,
      "i" = "`max_iter` must be either infinite (for hill_climbing optimizer) or a whole number.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    )
  } else if (max_iter < 2) {
    abort_text <- c(abort_text,
      "i" = "`max_iter` must be at least 2.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    )
  }
  
  if (!(permutations::is.cycle(start_perm) || inherits(start_perm, "gips_perm"))) {
    abort_text <- c(abort_text,
      "i" = "`start_perm` must be the output of `gips_perm()` function, or of a `cycle` class from `permutations` package.",
      "x" = paste0(
        "You provided `start_perm` with `class(start_perm) == (",
        paste(class(start_perm), collapse = ", "),
        ")`."
      )
    )
  } else if (!(permutations::is.cycle(start_perm) || attr(start_perm, "size") == ncol(S))) {
    abort_text <- c(abort_text,
      "i" = "`start_perm` must have the `size` attribute equal to the shape of a square matrix `S`",
      "x" = paste0(
        "You provided `start_perm` with `size == ",
        attr(start_perm, "size"),
        "`, but the `S` matrix You provided has ",
        ncol(S), " columns."
      )
    )
  }
  
  if (is.null(delta)) {
    abort_text <- c(abort_text,
      "i" = "`delta` must not be `NULL`.",
      "x" = "Your provided `delta` is a `NULL`."
    )
  } else if (delta <= 1) {
    abort_text <- c(abort_text,
      "i" = "`delta` must be strictly bigger than 1.",
      "x" = paste0("You provided `delta == ", delta, "`.")
    )
  }
  
  if (!(is.null(D_matrix) || is.matrix(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must either be `NULL` or a matrix.",
      "x" = paste0(
        "You provided `D_matrix` with type ",
        typeof(D_matrix), "."
      )
    )
  } else if (!(is.null(D_matrix) || ncol(D_matrix) == nrow(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must either be `NULL` or a square matrix.",
      "x" = paste0(
        "You provided `D_matrix` as a matrix, but with different sizes: ",
        ncol(D_matrix), " and ", nrow(D_matrix), "."
      )
    )
  } else if (!(is.null(D_matrix) || ncol(S) == ncol(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`S` must be a square matrix with the same shape as a square matrix `D_matrix`.",
      "x" = paste0(
        "You provided `S` with shape ",
        ncol(S), " and ", nrow(S),
        ", but also `D_matrix` with shape ",
        ncol(D_matrix), " and ", nrow(D_matrix), "."
      )
    )
  } else if (any(is.nan(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must not contain any `NaN`s.",
      "x" = "You provided `D_matrix` with `NaN`s!"
    )
  } else if (any(is.infinite(D_matrix))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrix` must not contain any infinite values.",
      "x" = "You provided `D_matrix` with infinite values!"
    )
  }
  
  abort_text <- c(abort_text, check_logical_flag(was_mean_estimated, "was_mean_estimated"))
  abort_text <- c(abort_text, check_logical_flag(return_probabilities, "return_probabilities"))
  abort_text <- c(abort_text, check_logical_flag(save_all_perms, "save_all_perms"))
  abort_text <- c(abort_text, check_logical_flag(show_progress_bar, "show_progress_bar"))
  
  if (length(abort_text) > 0) {
    abort_text <- c(
      paste0(
        "There were ", (length(abort_text) - additional_info) / 2,
        " problems identified with the provided arguments:"
      ),
      abort_text
    )

    if (length(abort_text) > 11) {
      abort_text <- c(
        abort_text[1:11],
        paste0("... and ", (length(abort_text) - 1) / 2 - 5, " more problems")
      )
    }

    abort_text <- c(
      abort_text,
      ">" = "If You think You've found a bug in a package, please open an ISSUE on https://github.com/PrzeChoj/gips/issues"
    )

    rlang::abort(abort_text)
  }

  if (return_probabilities && !save_all_perms) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "For calculations of probabilities, all perms have to be available after the optimization process.",
      "x" = "You provided `return_probabilities == TRUE` and `save_all_perms == FALSE`!",
      "i" = "Did You want to set `save_all_perms = TRUE`?",
      "i" = "Did You want to set `return_probabilities = FALSE`?",
      "!" = paste0(
        "Remember that setting `return_probabilities == TRUE` can be computationally costly",
        ifelse(show_progress_bar, " and second progress bar will be shown.", ".")
      )
    ))
  }
  
  invisible(NULL)
}


# ---- Multi-sample helpers ---------------------------------------------------

#' Check whether x is a list of matrices
#'
#' @noRd
list_of_matrices_check <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  }
  for (y in x) {
    if (!is.matrix(y)) {
      return(FALSE)
    }
  }
  return(TRUE)
}


#' Check that Ss, D_matrices, and numbers_of_observations are compatible
#'
#' All three must have the same length G, and each D_matrix element must be
#' the same dimension as the corresponding S element.
#'
#' @noRd
SDN_compatibility_check <- function(Ss, D_matrices, numbers_of_observations) {
  n <- length(numbers_of_observations)
  if (length(Ss) != n || length(D_matrices) != n) {
    return(FALSE)
  }
  if (!all(dim(D_matrices[[1]]) == dim(Ss[[1]]))) {
    return(FALSE)
  }
  return(TRUE)
}


#' Validate arguments for the multi-sample gips constructor and validator
#'
#' @noRd
check_gips_mult_arguments <- function(Ss, numbers_of_observations, delta,
                                       D_matrices, was_mean_estimated, perm) {
  if (!list_of_matrices_check(Ss)) {
    rlang::abort(c(
      "There was a problem identified with the provided S argument:",
      "i" = "For multi-sample `gips`, `S` must be a list of matrices.",
      "x" = "Not all elements of `S` are valid matrices."
    ))
  }

  G <- length(Ss)
  p <- ncol(Ss[[1]])

  if (!all(sapply(Ss, ncol) == p) || !all(sapply(Ss, nrow) == p)) {
    rlang::abort(c(
      "There was a problem identified with the provided S argument:",
      "i" = "All matrices in `S` must have the same dimensions.",
      "x" = "Not all matrices in `S` have the same number of rows and columns."
    ))
  }

  if (!is.numeric(numbers_of_observations) || length(numbers_of_observations) != G) {
    rlang::abort(c(
      "There was a problem identified with the provided arguments:",
      "i" = "`number_of_observations` must have one whole-number entry per matrix in `S`.",
      "x" = paste0(
        "You provided `S` with ", G, " matrices but `number_of_observations` has length ",
        length(numbers_of_observations), "."
      )
    ))
  }

  # delta must be either a scalar or a vector of length G; validate length
  if (!is.numeric(delta) || !(length(delta) == 1 || length(delta) == G)) {
    rlang::abort(c(
      "There was a problem identified with the provided delta argument:",
      "i" = "`delta` must be a single number or a numeric vector with one entry per matrix in `S`.",
      "x" = paste0(
        "You provided `S` with ", G, " matrices but `delta` has length ",
        length(delta), "."
      )
    ))
  }
  # Broadcast scalar to vector before per-group validation
  if (length(delta) == 1) delta <- rep(delta, G)

  if (!is.null(D_matrices)) {
    if (!list_of_matrices_check(D_matrices)) {
      rlang::abort(c(
        "There was a problem identified with the provided D_matrix argument:",
        "i" = "For multi-sample `gips`, `D_matrix` must be a list of matrices.",
        "x" = "Not all elements of `D_matrix` are valid matrices."
      ))
    }
    if (!SDN_compatibility_check(Ss, D_matrices, numbers_of_observations)) {
      rlang::abort(c(
        "There was a problem identified with the provided arguments:",
        "i" = "The lists `S`, `D_matrix`, and `number_of_observations` must all have the same length and consistent dimensions.",
        "x" = "Incompatible lengths or dimensions detected."
      ))
    }
  }

  # Validate each group using the single-sample validator
  for (i in seq_len(G)) {
    check_gips_arguments(
      S = Ss[[i]],
      number_of_observations = numbers_of_observations[i],
      delta = delta[i],
      D_matrix = if (is.null(D_matrices)) NULL else D_matrices[[i]],
      was_mean_estimated = was_mean_estimated,
      perm = perm
    )
  }

  invisible(NULL)
}
