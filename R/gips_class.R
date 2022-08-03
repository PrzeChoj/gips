#' Constructor of the `gips` class.
#'
#' Create the `gips` object.
#' This object will consists data and all other information needed to find the invariant group.
#' The optimization itself will not be performed. To do it, one have to call the \code{\link{find_gips}} function. See examples below.
#'
#' @param S matrix, estimated covariance matrix. When Z is observed data: `S = (t(Z) %*% Z) / number_of_observations`, if one know the theoretical mean is 0; # TODO(What if one have to estimate the theoretical mean with the empirical mean)
#' @param number_of_observations number of data points that `S` is based on.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of the same size as `S`. When NULL, the identity matrix is taken.
#' @param perm optional permutation to be the base for `gips` object. Can be of the class `gips_perm` or `permutation` or anything the function `permutations::permutation()` can take.
#'
#' @return Object of class gips.
#'
#' @export
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 6
#' mu <- numeric(perm_size)
#' # sigma is a matrix invariant under permutation (1,2,3,4,5,6)
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
#' )
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- (t(Z) %*% Z) / number_of_observations # the theoretical mean is 0
#'
#' g <- gips(S, number_of_observations)
#'
#' g_opt <- find_gips(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "MH")
#' g_opt
#'
#' if (require("graphics")) {
#'   plot(g_opt, type = "both", logarithmic_x = TRUE)
#' }
gips <- function(S, number_of_observations, delta = 3, D_matrix = NULL,
                 perm = "") {
  if (!inherits(perm, c("gips_perm", "permutation"))) {
    perm <- permutations::permutation(perm)
  }

  check_correctness_of_arguments( # max_iter, return_probabilities and show_progress_bar are to be checked here, but some value has to be passed
    S = S, number_of_observations = number_of_observations,
    max_iter = 2, start_perm = perm,
    delta = delta, D_matrix = D_matrix,
    return_probabilities = FALSE, show_progress_bar = FALSE
  )

  if (inherits(perm, "gips_perm")) {
    gips_perm_object <- perm # it is already a `gips_perm`
  } else {
    gips_perm_object <- gips_perm(perm, nrow(S)) # it is of class `cycle` from permutations package (it was checked in `check_correctness_of_arguments()`. Make it 'gips_perm'
  }


  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = ncol(S))
  }

  validate_gips(new_gips(
    list(gips_perm_object), S, number_of_observations,
    delta = delta, D_matrix = D_matrix,
    optimization_info = NULL
  ))
}


# TODO(Documentation)
#' @describeIn gips Constructor
#'
#' Only intended for low-level use.
#'
#' @param list_of_gips_perm list with a single element of class `gips_perm`. The base object for the `gips` class.
#' @param optimization_info NULL or the list with information about the optimization process.
#'
#' @export
new_gips <- function(list_of_gips_perm, S, number_of_observations, delta,
                     D_matrix, optimization_info) {
  if (!is.list(list_of_gips_perm) ||
    !inherits(list_of_gips_perm[[1]], "gips_perm") ||
    !is.matrix(S) ||
    !is.wholenumber(number_of_observations) ||
    !is.numeric(delta) ||
    !is.matrix(D_matrix) ||
    !(is.null(optimization_info) || is.list(optimization_info))) {
    rlang::abort(c("x" = "`gips` object cannot be created from those arguments."))
  }


  structure(list_of_gips_perm,
    S = S, number_of_observations = number_of_observations,
    delta = delta, D_matrix = D_matrix,
    optimization_info = optimization_info,
    class = c("gips")
  )
}


# TODO(Documentation; I think @describeIn is intended for methods...)
#' @describeIn gips Validator
#'
#' Only intended for low-level use.
#'
#' @param g Element to be checked if it is of class `gips`.
#'
#' @export
validate_gips <- function(g) {
  if (!(inherits(g, "gips"))) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "`g` must be of class `gips`.",
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
  optimization_info <- attr(g, "optimization_info")

  if (!inherits(perm, "gips_perm")) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `g[[1]]` must be an object of a `gips_perm` class.",
      "x" = paste0(
        "You provided `g[[1]]` with `class(g[[1]]) == (",
        paste(class(perm), collapse = ", "),
        ")`."
      )
    ))
  } else {
    tryCatch(
      {
        validate_gips_perm(perm)
      },
      error = function(cond) {
        rlang::abort(c("There was a problem identified with provided argument:",
          "i" = "The `g[[1]]` must be an object of a `gips_perm` class.",
          "x" = paste0(
            "You provided `g[[1]]` with `class(g[[1]]) == 'gips_perm'`, but your g[[1]] does not pass `validate_gips_perm(g[[1]])`."
          )
        ))
      }
    )
  }

  check_correctness_of_arguments( # max_iter, return_probabilities and show_progress_bar are to be checked here, but some value has to be passed
    S = S, number_of_observations = number_of_observations,
    max_iter = 2, start_perm = perm,
    delta = delta, D_matrix = D_matrix,
    return_probabilities = FALSE, show_progress_bar = FALSE
  )

  if (!(is.null(optimization_info) || is.list(optimization_info))) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `optimization_info` value must be either a NULL, or a list.",
      "x" = paste0(
        "You provided `attr(g, 'optimization_info')` with type ",
        typeof(optimization_info), "."
      )
    ))
  }

  if (is.list(optimization_info)) { # Validate the `optimization_info` after the optimization
    legal_fields <- c("acceptance_rate", "log_likelihood_values", "visited_perms", "last_perm", "last_perm_log_likelihood", "iterations_performed", "optimization_algorithm_used", "post_probabilities", "did_converge", "best_perm_log_likelihood", "optimization_time", "full_optimization_time")

    lacking_fields <- setdiff(legal_fields, names(optimization_info))
    illegal_fields <- setdiff(names(optimization_info), legal_fields)

    abort_text <- character(0)

    if (!(length(lacking_fields) == 0)) {
      abort_text <- c("x" = paste0(
        "Your `attr(g, 'optimization_info')` lacks the following fields: ",
        paste(lacking_fields, collapse = ", "), "."
      ))
    }
    if (!(length(illegal_fields) == 0)) {
      abort_text <- c(abort_text,
        "x" = paste0(
          "Your `attr(g, 'optimization_info')` has the following, unexpected fields: ",
          paste(illegal_fields, collapse = ", "), "."
        )
      )
    }

    # abort the validation
    if (length(abort_text) > 0) {
      rlang::abort(c("There was a problem with the 'optimization_info' attribute.",
        "i" = paste0(
          "After optimiation, `attr(g, 'optimization_info')` must be a list of 11 elements with names: ",
          paste(legal_fields, collapse = ", "), "."
        ),
        "x" = paste0("You have a list of ", length(names(optimization_info)), " elements."),
        abort_text,
        "i" = "Did You accidentally edited `attr(g, 'optimization_info')` by yourself?",
        "i" = "Did You accidentally set one of `attr(g, 'optimization_info')` elements to `NULL` or `NA`?"
      ))
    }

    # All the fields as named as they should be. Check if their content are as expected:
    abort_text <- character(0)
    additional_info <- 0 # for calculation of the number of problems
    if (!((is.numeric(optimization_info[["acceptance_rate"]]) &&
      (length(optimization_info[["acceptance_rate"]]) == 1) &&
      optimization_info[["acceptance_rate"]] >= 0 &&
      optimization_info[["acceptance_rate"]] <= 1) ||
      optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] == "brute_force")) { # when brute_force, acceptance_rate is NULL
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['acceptance_rate']]` must be a number in range [0, 1].",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['acceptance_rate']] == (",
          paste(optimization_info[["acceptance_rate"]], collapse = ", "),
          ")`."
        )
      )
    }
    if (!(optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force" ||
      is.null(optimization_info[["acceptance_rate"]]))) {
      abort_text <- c(abort_text,
        "i" = "When brute force algorithm was used for optimization, `attr(g, 'optimization_info')[['acceptance_rate']]` must be a NULL.",
        "x" = paste0(
          "You have used brute force algorithm, but `attr(g, 'optimization_info')[['acceptance_rate']] == (",
          paste(optimization_info[["acceptance_rate"]], collapse = ", "),
          ")`."
        )
      )
    }
    if (!(is.numeric(optimization_info[["log_likelihood_values"]]))) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['log_likelihood_values']]` must be a vector of numbers.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['log_likelihood_values']] == ",
          typeof(optimization_info[["log_likelihood_values"]]),
          "`."
        )
      )
    }
    if (optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force") {
      if (!(is.list(optimization_info[["visited_perms"]]))) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['visited_perms']]` must be a list.",
          "x" = paste0(
            "You have `attr(g, 'optimization_info')[['visited_perms']]` of type ",
            typeof(optimization_info[["visited_perms"]]),
            "."
          )
        )
      } else if (!(inherits(optimization_info[["visited_perms"]][[1]], "gips_perm"))) { # It only checks for the first one, because checking for every would be too expensive
        abort_text <- c(abort_text,
          "i" = "Elements of `attr(g, 'optimization_info')[['visited_perms']]` must be of class `gips_perm`.",
          "x" = paste0(
            "You have `class(attr(g, 'optimization_info')[['visited_perms']][[1]]) == (",
            paste(class(optimization_info[["visited_perms"]][[1]]), collapse = ", "),
            ")`."
          )
        )
      } else if (!(identical(optimization_info[["last_perm"]], optimization_info[["visited_perms"]][[length(optimization_info[["visited_perms"]])]]))) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['last_perm']]` must be the last element of `attr(g, 'optimization_info')[['visited_perms']]` list.",
          "x" = paste0("You have `attr(g, 'optimization_info')[['last_perm']]` different from `attr(g, 'optimization_info')[['visited_perms']][[length(attr(g, 'optimization_info')[['visited_perms']])]]`.")
        )
      }

      if (inherits(optimization_info[["last_perm"]], "gips_perm")) {
        abort_text <- tryCatch(
          {
            validate_gips_perm(optimization_info[["last_perm"]])

            # optimization_info[["last_perm"]] is proper gips_perm object
            last_perm_gips <- gips(S, number_of_observations,
              delta = delta, D_matrix = D_matrix,
              perm = optimization_info[["last_perm"]]
            )

            if (!(optimization_info[["last_perm_log_likelihood"]] == log_likelihood_of_gips(last_perm_gips))) {
              abort_text <- c(abort_text,
                "i" = "`attr(g, 'optimization_info')[['last_perm_log_likelihood']]` must be the log_likelihood of `optimization_info[['last_perm']]`.",
                "x" = paste0(
                  "You have `attr(g, 'optimization_info')[['last_perm_log_likelihood']] == ",
                  optimization_info[["last_perm_log_likelihood"]],
                  "`, but `log_likelihood_of_gips(gips(attr(g, 'S'), attr(g, 'number_of_observations'), delta=attr(g, 'delta'), D_matrix=attr(g, 'D_matrix'), perm=attr(g, 'optimization_info')[['last_perm']])) == ",
                  log_likelihood_of_gips(last_perm_gips), "`."
                )
              )
            }

            abort_text # if optimization_info[["last_perm"]] passes the validation, return original text
          },
          error = function(cond) { # error can only be thrown in `validate_gips_perm()`, so ass appropriate note to the `abort_text`:
            c(abort_text,
              "i" = "The `attr(g, 'optimization_info')[['last_perm']]` must be an object of a `gips_perm` class.",
              "x" = paste0(
                "You provided `attr(g, 'optimization_info')[['last_perm']]` with `class(attr(g, 'optimization_info')[['last_perm']]) == 'gips_perm'`, but your attr(g, 'optimization_info')[['last_perm']] does not pass `validate_gips_perm(attr(g, 'optimization_info')[['last_perm']])`."
              )
            )
          }
        )
      } else {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['last_perm']]` must be an object of class 'gips_perm'.",
          "x" = paste0(
            "You have `attr(g, 'optimization_info')[['last_perm']]` of class ('",
            paste(class(optimization_info[["last_perm"]]), collapse = "', '"), "')."
          )
        )
      }
    } else { # for brute_force, the visited_perms are of class "cycle"
      if (!(is.list(optimization_info[["visited_perms"]]))) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['visited_perms']]` must be a list.",
          "x" = paste0(
            "You have `attr(g, 'optimization_info')[['visited_perms']]` of type ",
            typeof(optimization_info[["visited_perms"]]),
            "."
          )
        )
      } else if (!(inherits(optimization_info[["visited_perms"]][[1]], "list"))) { # It only checks for the first one, because checking for every would be too expensive
        abort_text <- c(abort_text,
          "i" = "After optimization with brute force algorithm, elements of `attr(g, 'optimization_info')[['visited_perms']]` must be of class `list`.",
          "x" = paste0(
            "You have `class(attr(g, 'optimization_info')[['visited_perms']][[1]]) == (",
            paste(class(optimization_info[["visited_perms"]][[1]]), collapse = ", "),
            ")`."
          )
        )
      } else if (!(is.null(optimization_info[["last_perm"]]))) {
        abort_text <- c(abort_text,
          "i" = "After optimization with brute force algorithm, `attr(g, 'optimization_info')[['last_perm']]` must be a NULL.",
          "x" = paste0("You have `attr(g, 'optimization_info')[['last_perm']]` of type ", typeof(optimization_info[["last_perm"]]), ".")
        )
      }
    }

    if (!(all(is.wholenumber(optimization_info[["iterations_performed"]])))) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['iterations_performed']]` must be a vector of whole numbers.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['iterations_performed']] == (",
          paste(optimization_info[["iterations_performed"]], collapse = ", "),
          ")`."
        )
      )
    } else if (!(sum(optimization_info[["iterations_performed"]]) <= length(optimization_info[["log_likelihood_values"]]))) {
      abort_text <- c(abort_text,
        "i" = "In every iteration at least one value of log_likelihood is calculated.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['iterations_performed']] == ",
          sum(optimization_info[["iterations_performed"]]),
          "`, which is more than `length(attr(g, 'optimization_info')[['log_likelihood_values']]) == ",
          length(optimization_info[["log_likelihood_values"]]), "`."
        )
      )
    }
    if (!all(optimization_info[["optimization_algorithm_used"]] %in% c("Metropolis_Hastings", "best_growth", "brute_force"))) { # Even if MH was used, it would produce the text "Metropolis_Hastings"
      abort_text <- c(abort_text,
        "i" = "The available optimization algorithms are 'Metropolis_Hastings', 'best_growth' and 'brute_force'.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_algorithm_used']] == (",
          paste(optimization_info[["optimization_algorithm_used"]], collapse = ", "),
          ")`."
        )
      )
    } else if (all(optimization_info[["optimization_algorithm_used"]] != "Metropolis_Hastings") &&
      !is.null(optimization_info[["post_probabilities"]])) {
      abort_text <- c(abort_text,
        "i" = "`post_probabilities` can olny be obtained with 'Metropolis_Hastings' optimization method.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_algorithm_used']] == ",
          optimization_info[["optimization_algorithm_used"]],
          "` and the `attr(g, 'optimization_info')[['post_probabilities']]` is not NULL, but is of type `",
          typeof(optimization_info[["post_probabilities"]]), "`."
        )
      )
    } else if (!(is.null(optimization_info[["post_probabilities"]]) ||
      length(optimization_info[["post_probabilities"]]) <= length(optimization_info[["visited_perms"]]))) {
      abort_text <- c(abort_text,
        "i" = "Every element of `attr(g, 'optimization_info')[['post_probabilities']]` was taken from a visided permutation, so it is in `attr(g, 'optimization_info')[['visited_perms']]`.",
        "x" = paste0(
          "You have `length(attr(g, 'optimization_info')[['visited_perms']]) == ",
          length(optimization_info[["post_probabilities"]]),
          "`, but `length(attr(g, 'optimization_info')[['post_probabilities']]) == ",
          length(optimization_info[["visited_perms"]]),
          "` which are not equal."
        )
      )
    } else if (!(is.null(optimization_info[["post_probabilities"]]) ||
      (all(optimization_info[["post_probabilities"]] <= 1) &&
        all(optimization_info[["post_probabilities"]] > 0) &&
        (sum(optimization_info[["post_probabilities"]]) < 1.001) && # Allow small error
        (sum(optimization_info[["post_probabilities"]]) > 0.999)))) {
      abort_text <- c(abort_text,
        "i" = "The vector of `attr(g, 'optimization_info')[['post_probabilities']]` must have properties of probability. All elements in range [0, 1] and sums to 1. What is more, every element of `attr(g, 'optimization_info')[['post_probabilities']]` was visided, so has to have post_probability bigger than 0.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['post_probabilities']]` in a range [",
          min(optimization_info[["post_probabilities"]]), ",",
          max(optimization_info[["post_probabilities"]]),
          "] and with the sum ",
          sum(optimization_info[["post_probabilities"]]),
          "."
        )
      )
    }
    if ((!(optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] %in% c("best_growth", "brute_force"))) && # The last optimization_algorithm_used has to be best_growth or brute_force to make the convergence
      !is.null(optimization_info[["did_converge"]])) {
      abort_text <- c(abort_text,
        "i" = "`did_converge` can only be obtained with 'best_growth' or 'brute_force' optimization method.",
        "x" = paste0(
          "The last optimization method You used was `attr(g, 'optimization_info')[['optimization_algorithm_used']][length(attr(g, 'optimization_info')[['optimization_algorithm_used']])] == ",
          optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])],
          "` and the `attr(g, 'optimization_info')[['did_converge']]` is not NULL, but is of type ",
          typeof(optimization_info[["did_converge"]]), "."
        )
      )
    } else if ((optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] == "best_growth") &&
      !is.logical(optimization_info[["did_converge"]])) {
      abort_text <- c(abort_text,
        "i" = "When 'best_growth' optimization method, the `did_converge` must be TRUE or FALSE.",
        "x" = paste0(
          "The last optimization method You used was `attr(g, 'optimization_info')[['optimization_algorithm_used']][length(attr(g, 'optimization_info')[['optimization_algorithm_used']])] == ",
          optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])],
          "` and the `attr(g, 'optimization_info')[['did_converge']]` is not of type logical, but it is of type ",
          typeof(optimization_info[["did_converge"]]), "."
        )
      )
    } else if ((optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] == "best_growth") &&
      is.na(optimization_info[["did_converge"]])) {
      abort_text <- c(abort_text,
        "i" = "When 'best_growth' optimization method, the `did_converge` must be TRUE or FALSE.",
        "x" = paste0(
          "The last optimization method You used was `attr(g, 'optimization_info')[['optimization_algorithm_used']][length(optimization_info[['optimization_algorithm_used']])] == ",
          optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])],
          "` and the `attr(g, 'optimization_info')[['did_converge']]` is of type logical, but it is a NA."
        )
      )
    }
    best_perm_gips <- gips(S, number_of_observations, delta = delta, D_matrix = D_matrix, perm = perm) # this perm is g[[1]]
    if (!(optimization_info[["best_perm_log_likelihood"]] == log_likelihood_of_gips(best_perm_gips))) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['best_perm_log_likelihood']]` must be the log_likelihood of the base object, `g[[1]]`.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['best_perm_log_likelihood']] == ",
          optimization_info[["best_perm_log_likelihood"]],
          "`, but `log_likelihood_of_gips(gips(attr(g, 'S'), attr(g, 'number_of_observations'), delta=attr(g, 'delta'), D_matrix=attr(g, 'D_matrix'), perm=g[[1]])) == ",
          log_likelihood_of_gips(best_perm_gips), "`."
        )
      )
    }
    if (any(is.na(optimization_info[["optimization_time"]]))) {
      additional_info <- additional_info + 2
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['optimization_time']]` is initially set to NA, but that state of the gips object should not be available to the user.",
        "x" = "You have `is.na(attr(g, 'optimization_info')[['optimization_time']]) == TRUE`.",
        "i" = "Did You used the inner optimizers like `gips:::Metropolis_Hastings()` or `gips:::best_growth()` in stead of the exported function `gips::find_gips()`?",
        "i" = "Did You modified the `find_gips()` function?"
      )
    } else if (!inherits(optimization_info[["optimization_time"]], "difftime")) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['optimization_time']]` has to be of a class 'difftime'.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_time']]` of a class (",
          paste0(class(optimization_info[["optimization_time"]]), collapse = ", "), ")."
        )
      )
    } else if (any(optimization_info[["optimization_time"]] <= 0)) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['optimization_time']]` has to be a time difference bigger than 0.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_time']] == ",
          optimization_info[["optimization_time"]], "`."
        )
      )
    }
    if (is.na(optimization_info[["full_optimization_time"]])) {
      additional_info <- additional_info + 2
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['full_optimization_time']]` is initially set to NA, but that state of the gips object should not be available to the user.",
        "x" = "You have `is.na(attr(g, 'optimization_info')[['full_optimization_time']]) == TRUE`.",
        "i" = "Did You used the inner optimizers like `gips:::Metropolis_Hastings()` or `gips:::best_growth()` in stead of the exported function `gips::find_gips()`?",
        "i" = "Did You modified the `find_gips()` function?"
      )
    } else if (!inherits(optimization_info[["full_optimization_time"]], "difftime")) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['full_optimization_time']]` has to be of a class 'difftime'.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['full_optimization_time']]` of a class (",
          paste0(class(optimization_info[["full_optimization_time"]]), collapse = ", "), ")."
        )
      )
    } else if (optimization_info[["full_optimization_time"]] <= 0) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['full_optimization_time']]` has to be a time difference bigger than 0.",
        "x" = paste0(
          "You have `attr(g, 'full_optimization_info')[['optimization_time']] == ",
          optimization_info[["full_optimization_time"]], "`."
        )
      )
    }

    # TODO(Validate that there is the same number of algorithms used and length of other things dependent on it)


    if (length(abort_text) > 0) {
      if (!is.wholenumber((length(abort_text) - additional_info) / 2)) {
        rlang::inform(paste0(
          "You found a small bug in gips package. We calculated there was ",
          (length(abort_text) - additional_info) / 2,
          " problems, but it is not a whole number. Please inform us about that bug by opening the ISSUE on https://github.com/PrzeChoj/gips/issues"
        ))
      }
      abort_text <- c(
        paste0(
          "There were ", (length(abort_text) - additional_info) / 2,
          " problems identified with `attr(g, 'optimization_info')`:"
        ),
        abort_text
      )

      if (length(abort_text) > 11) {
        abort_text <- c(abort_text[1:11],
          "x" = paste0("... and ", (length(abort_text) - 1) / 2 - 5, " more problems")
        )
      }

      abort_text <- c(abort_text,
        "i" = "Did You accidentally edited `attr(g, 'optimization_info')` by yourself?",
        ">" = "If You think You've found a bug in a package, please open the ISSUE on https://github.com/PrzeChoj/gips/issues"
      )

      rlang::abort(abort_text)
    }
  }

  g
}


check_correctness_of_arguments <- function(S, number_of_observations, max_iter,
                                           start_perm, delta, D_matrix,
                                           return_probabilities, show_progress_bar) {
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
  additional_info <- 0 # for calculation of the number of problems
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
  } else if (sum(S == t(S)) != (nrow(S)^2)) { # this would mean the matrix is not symmetric
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a symmetric matrix.",
      "x" = "You provided `S` as a matrix, but a non-symmetric one.",
      "i" = "Is your matrix approximatelly symmetric? Maybe try setting `S <- (S+t(S))/2`?"
    )
    additional_info <- additional_info + 1 # for calculation of the number of problems
  } else if (!is.positive.semi.definite.matrix(S, tolerance = 1e-06)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be positive semi-definite matrix.",
      "x" = "You provided `S` as a symmetric matrix, but a non-positive-semi-definite one."
    )
  }
  if (is.null(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must not be `NULL`.",
      "x" = "Your provided `number_of_observations` is NULL."
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
      "i" = "`max_iter` must be either infinite (for best_growth optimizer) or a whole number.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    )
  } else if (max_iter < 2) { # TODO(Make it work for max_iter == 1)
    abort_text <- c(abort_text,
      "i" = "`max_iter` must be at least 2.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    )
  }
  if (!(permutations::is.cycle(start_perm) || inherits(start_perm, "gips_perm"))) {
    abort_text <- c(abort_text,
      "i" = "`start_perm` must be the output of `gips_perm()` function, or of class `cycle` form `permutations` package.", # this is not true, but it is close enough
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
      "x" = "Your provided `delta` is a NULL."
    )
  } else if (delta <= 2) {
    abort_text <- c(abort_text,
      "i" = "`delta` must be strictly bigger than 2.",
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
  }
  if (!is.logical(return_probabilities)) {
    abort_text <- c(abort_text,
      "i" = "`return_probabilities` must be a logic value (TRUE or FALSE).",
      "x" = paste0(
        "You provided `return_probabilities` with type ",
        typeof(return_probabilities), "."
      )
    )
  }
  if (!is.logical(show_progress_bar)) {
    abort_text <- c(abort_text,
      "i" = "`show_progress_bar` must be a logic value (TRUE or FALSE).",
      "x" = paste0(
        "You provided `show_progress_bar` with type ",
        typeof(show_progress_bar), "."
      )
    )
  }

  if (length(abort_text) > 0) {
    abort_text <- c(
      paste0(
        "There were ", (length(abort_text) - additional_info) / 2,
        " problems identified with provided arguments:"
      ),
      abort_text
    )

    if (length(abort_text) > 11) {
      abort_text <- c(
        abort_text[1:11],
        paste0("... and ", (length(abort_text) - 1) / 2 - 5, " more problems"),
        ">" = "If You think You've found a bug in a package, please open the ISSUE on https://github.com/PrzeChoj/gips/issues"
      )
    }

    rlang::abort(abort_text)
  }
}



#' Printing `gips` object
#'
#' Printing function for `gips` class.
#'
#' @param x object of class `gips`.
#' @param log_value logical. Weather to print the exp of a value of a \code{\link{log_likelihood_of_gips}} or leave it in logarithmic form.
#' @param digits Number of digits after the comma for likelihood to be presented. Can be negative. Be default, `Inf`. It is passed to \code{\link{round}}.
#' @param ... additional arguments passed to \code{\link{cat}}.
#'
#' @return Invisible NULL.
#' @export
print.gips <- function(x, log_value = TRUE, digits = Inf, ...) {
  # TODO(it is not likelihood, but sth proportional to it. See #ISSUE11)
  validate_gips(x)

  if (is.null(attr(x, "optimization_info"))) { # it is unoptimized gips object
    log_likelihood <- log_likelihood_of_perm(
      perm_proposal = x[[1]], S = attr(x, "S"),
      number_of_observations = attr(x, "number_of_observations"),
      delta = attr(x, "delta"), D_matrix = attr(x, "D_matrix")
    )
    if (is.nan(log_likelihood) || is.infinite(log_likelihood)) {
      # See ISSUE#5; We hope the introduction of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The likelihood value of ", ifelse(is.nan(log_likelihood), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))
    }
    cat(paste0(
      "The permutation ", x[[1]],
      ifelse(log_value,
        paste0(
          " has log likelihood ",
          round(log_likelihood, digits = digits)
        ),
        paste0(
          " has likelihood ",
          round(exp(log_likelihood), digits = digits)
        )
      ), "."
    ), ...)
  } else { # it is optimized gips object
    log_likelihood <- attr(x, "optimization_info")[["best_perm_log_likelihood"]]
    if (is.nan(log_likelihood) || is.infinite(log_likelihood)) {
      # See ISSUE#5; We hope the introduction of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The likelihood value of ", ifelse(is.nan(log_likelihood), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))
    }
    cat(paste0(
      "The permutation ", x[[1]],
      ifelse(log_value,
        paste0(
          " has log likelihood ",
          round(log_likelihood, digits = digits)
        ),
        paste0(
          " has likelihood ",
          round(exp(log_likelihood), digits = digits)
        )
      ), " which was found after ",
      length(attr(x, "optimization_info")[["log_likelihood_values"]]),
      " log_likelihood calculations."
    ), ...)
  }
}




#' Plot optimized matrix or optimization convergence
#'
#' Plot method for `gips` object.
#'
#' @param x Object of class gips. Has to first be optimized with \code{\link{find_gips}}.
#' @param type Character. A type of a plot. one of \code{c("heatmap", "all", "best", "both")}. For "heatmap", plots a heatmap of the `S` matrix inside the `gips` object that was projected on the permutation in the `gips` object. For "all", plots the line of the likelihoods for all visited state. For "best", plots the line of the biggest likelihoods up to the moment For "both", both lines from "all" and "best" are plotted. Default value is `NA`, which will be changed to "heatmap" for non-optimized `gips` objects, and to "both" for optimized ones. Using the default produces a warning. For the `type = "heatmap"`, all other arguments are ignored.
#' @param logarithmic_y boolean.
#' @param logarithmic_x boolean.
#' @param color Vector of colors to be used to plot lines.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend boolean.
#' @param ylim Limits of y axis. When \code{NULL}, the minimum and maximum of the \code{\link{log_likelihood_of_gips}} is taken.
#' @param xlim Limits of x axis. When \code{NULL}, the whole optimization process is shown.
#' @param ... Additional arguments passed to \code{\link{heatmap}} or other various elements of the plot.
#'
#' @return Invisible NULL.
#' @export
plot.gips <- function(x, type = NA,
                      logarithmic_y = TRUE, logarithmic_x = FALSE,
                      color = NULL,
                      title_text = "Convergence plot",
                      xlabel = NULL, ylabel = NULL,
                      show_legend = TRUE,
                      ylim = NULL, xlim = NULL, ...) {
  # TODO(It is not likelihood, but sth proportional to it. See #ISSUE11)
  # TODO(For "MH", those are NOT "All calculated likelihoods", but those that MH was in. Change the legend, or the output of `gips(type="MH")`)

  # checking the correctness of the arguments:
  if (!requireNamespace("graphics", quietly = TRUE)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "Package 'graphics' must be installed to use this function.",
      "x" = "Package 'graphics' seems to be unavailable."
    ))
  }

  validate_gips(x)

  if (is.na(type)) {
    type <- ifelse(is.null(attr(x, "optimization_info")),
      "heatmap",
      "both"
    )

    rlang::inform(c("You used the default value of the 'type' argument in `plot()` for gips object.",
      "i" = paste0(
        "The `type = NA` was automatically changed to `type = '",
        type, "'`."
      )
    ))
  }

  if (!(type %in% c("heatmap", "all", "best", "both"))) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "`type` must be one of: c('heatmap', 'all', 'best', 'both').",
      "x" = paste0("You provided `type == ", type, "`."),
      "i" = "Did You misspell the 'type' argument?"
    ))
  }

  if ((type != "heatmap") && is.null(attr(x, "optimization_info"))) {
    rlang::abort(c("There was a problem identified with provided arguments:",
      "i" = "For non-optimized `gips` objects only the `type = 'heatmap'` can be used.",
      "x" = paste0(
        "You did not optimized `x` and provided `type = '",
        type, "'`."
      ),
      "i" = paste0(
        "Did You want to call `x <- find_gips(g)` and then `plot(x, type = '",
        type, "')`?"
      ),
      "i" = "Did You want to use `type = 'heatmap'`?"
    ))
  }

  # plotting:
  if (type == "heatmap") {
    rlang::check_installed(c("dplyr", "tidyr", "tibble", "ggplot2"),
      reason = "to use `plot.gips(type = 'heatmap')`; without those packages, the `stats::heatmap()` will be used"
    )
    my_projected_matrix <- gips::project_matrix(attr(x, "S"), x[[1]])
    if (rlang::is_installed(c("dplyr", "tidyr", "tibble", "ggplot2"))) {
      p <- ncol(my_projected_matrix)

      if (is.null(colnames(my_projected_matrix))) {
        colnames(my_projected_matrix) <- paste0(seq(1, p))
      }
      if (is.null(rownames(my_projected_matrix))) {
        rownames(my_projected_matrix) <- paste0(seq(1, p))
      }

      # Address the R CMD check's "no visible binding for global variable":
      col_id <- covariance_value <- row_id <- NULL

      # Life would be easier with pipes (%>%)
      my_transformed_matrix <- tibble::rownames_to_column(
        as.data.frame(my_projected_matrix),
        "row_id"
      )
      my_transformed_matrix <- tidyr::pivot_longer(my_transformed_matrix,
        -c(row_id),
        names_to = "col_id",
        values_to = "covariance_value"
      )
      my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
        col_id = as.numeric(col_id)
      )
      my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
        row_id = as.numeric(row_id)
      )
      g_plot <- ggplot2::ggplot(
        my_transformed_matrix,
        ggplot2::aes(x = col_id, y = row_id, fill = covariance_value)
      ) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::scale_x_continuous(breaks = 1:p) +
        ggplot2::scale_y_reverse(breaks = 1:p) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = paste0("Covariance matrix projected on permutation ", x[[1]], "."),
          x = "", y = ""
        )

      return(g_plot)
    } else { # use the basic plot in R, package `graphics`
      if (is.null(color)) { # Setting col = NA or col = NULL turns off the whole plot.
        stats::heatmap(my_projected_matrix,
          symm = TRUE,
          Rowv = NA, Colv = NA, ...
        )
      } else {
        stats::heatmap(my_projected_matrix,
          symm = TRUE,
          Rowv = NA, Colv = NA, col = color, ...
        )
      }
    }
  }
  if (type %in% c("all", "best", "both")) {
    if (is.null(ylabel)) {
      ylabel <- ifelse(logarithmic_y,
        "log likelihood",
        "likelihood"
      )
    }
    if (is.null(xlabel)) {
      xlabel <- ifelse(logarithmic_x,
        "log10 of number of function calls",
        "number of function calls"
      )
    }
    if (is.null(color)) {
      if (type == "both") {
        color <- c("blue", "red")
      } else {
        color <- "blue"
      }
    }
    if (logarithmic_y) {
      y_values_from <- attr(x, "optimization_info")[["log_likelihood_values"]] # values of log_likelihood are logarithmic by default
    } else {
      y_values_from <- exp(attr(x, "optimization_info")[["log_likelihood_values"]])
    }

    y_values_max <- cummax(y_values_from)
    y_values_all <- y_values_from

    num_of_steps <- length(y_values_max)

    if (is.null(xlim)) {
      xlim <- c(1, num_of_steps)
    }

    if (is.null(ylim)) {
      ylim_plot <- c(min(y_values_from), y_values_max[num_of_steps])
      if (type == "best") {
        ylim_plot[1] <- y_values_from[1] # for the "best" type this is the smallest point of the graph
      }
    } else {
      ylim_plot <- ylim
    }

    # make the plot stairs-like
    x_points <- c(1, rep(2:num_of_steps, each = 2))

    if (logarithmic_x) {
      x_points <- log10(x_points)
      xlim <- log10(xlim)
    }

    graphics::plot.new()
    graphics::plot.window(xlim, ylim_plot)

    if (type != "best") {
      # make the plot stairs-like
      y_points <- c(
        rep(y_values_all[1:(length(y_values_all) - 1)], each = 2),
        y_values_all[length(y_values_all)]
      )

      graphics::lines.default(x_points, y_points,
        type = "l", lwd = 3,
        col = color[1], # the first color
        ...
      )
    }
    if (type != "all") {
      # make the plot stairs-like
      y_points <- c(
        rep(y_values_max[1:(length(y_values_max) - 1)], each = 2),
        y_values_max[length(y_values_max)]
      )

      graphics::lines.default(x_points, y_points,
        lwd = 3, lty = 1,
        col = color[length(color)], # the last color
        ...
      )
    }

    graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
    graphics::axis(1, ...)
    graphics::axis(2, ...)
    graphics::box(...)

    if (show_legend) {
      if (type == "both") {
        legend_text <- c(
          "All calculated likelihoods",
          "Maximum likelihoods calculated"
        )
        lty <- c(1, 1)
        lwd <- c(3, 3)
      } else if (type == "all") {
        legend_text <- c("All calculated function values")
        lty <- 1
        lwd <- 3
      } else if (type == "best") {
        legend_text <- c("Maximum function values calculated")
        lty <- 1
        lwd <- 3
      }

      graphics::legend("bottomright",
        inset = .002,
        legend = legend_text,
        col = color,
        lty = lty, lwd = lwd,
        cex = 0.7, box.lty = 0
      )
    }
  }

  invisible(NULL)
}


# TODO(Documentation. Can be based on `stats::summary.lm()`)
#' @export
summary.gips <- function(object, ...) {
  validate_gips(object)
  permutation_log_likelihood <- log_likelihood_of_gips(object)

  structure_constants <- get_structure_constants(object[[1]])
  n0 <- max(structure_constants[["r"]] * structure_constants[["d"]] / structure_constants[["k"]])

  if (is.null(attr(object, "optimization_info"))) {
    summary_list <- list(
      optimized = FALSE,
      start_permutation = object[[1]],
      start_permutation_log_likelihood = permutation_log_likelihood,
      n0 = n0,
      S_matrix = attr(object, "S"),
      number_of_observations = attr(object, "number_of_observations"),
      delta = attr(object, "delta"),
      D_matrix = attr(object, "D_matrix")
    )
  } else {
    optimization_info <- attr(object, "optimization_info")

    if (optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force") {
      start_permutation <- optimization_info[["visited_perms"]][[1]]
      start_permutation_log_likelihood <- log_likelihood_of_perm(
        start_permutation,
        attr(object, "S"),
        attr(object, "number_of_observations"),
        attr(object, "delta"),
        attr(object, "D_matrix")
      )

      when_was_best <- which(optimization_info[["log_likelihood_values"]] == permutation_log_likelihood)
    } else {
      start_permutation <- NULL
      start_permutation_log_likelihood <- NULL
      when_was_best <- NULL
    }


    summary_list <- list(
      optimized = TRUE,
      found_permutation = object[[1]],
      found_permutation_log_likelihood = permutation_log_likelihood,
      start_permutation = start_permutation,
      start_permutation_log_likelihood = start_permutation_log_likelihood,
      n0 = n0,
      S_matrix = attr(object, "S"),
      number_of_observations = attr(object, "number_of_observations"),
      delta = attr(object, "delta"),
      D_matrix = attr(object, "D_matrix"),
      optimization_algorithm_used = optimization_info[["optimization_algorithm_used"]],
      did_converge = optimization_info[["did_converge"]],
      number_of_log_likelihood_calls = length(optimization_info[["log_likelihood_values"]]),
      optimization_algorithm_used = optimization_info[["optimization_algorithm_used"]],
      full_optimization_time = optimization_info[["full_optimization_time"]],
      log_likelihood_calls_after_best = length(optimization_info[["log_likelihood_values"]]) - when_was_best[length(when_was_best)],
      acceptance_rate = optimization_info[["acceptance_rate"]]
    )
  }

  structure(summary_list,
    class = c("summary.gips")
  )
}

# Base on
# `sloop::s3_get_method(print.summary.lm)`
#' @export
print.summary.gips <- function(x, ...) {
  cat(ifelse(x[["optimized"]],
    "The optimized",
    "The unoptimized"
  ),
  " `gips` object.\n\nPermutation:\n ",
  ifelse(x[["optimized"]],
    as.character(x[["found_permutation"]]),
    as.character(x[["start_permutation"]])
  ),
  "\n\nLog_likelihood:\n ",
  ifelse(x[["optimized"]],
    x[["found_permutation_log_likelihood"]],
    x[["start_permutation_log_likelihood"]]
  ),
  "\n\nNumber of observations:\n ", x[["number_of_observations"]],
  "\n\nn0:\n ", x[["n0"]],
  "\n\nNumber of observations is ",
  ifelse(x[["n0"]] > x[["number_of_observations"]],
    "smaller",
    ifelse(x[["n0"]] == x[["number_of_observations"]],
      "equal",
      "bigger"
    )
  ),
  " than n0 for this permutaion, so the gips model based on the found permutation does ",
  ifelse(x[["n0"]] > x[["number_of_observations"]],
    "not ", ""
  ), "exist.",
  sep = ""
  )

  if (x[["optimized"]]) {
    cat("\n\n", paste0(rep("-", getOption("width")), collapse = ""), # the line
      sep = ""
    )

    if (length(x[["optimization_algorithm_used"]]) == 1) {
      cat("\nOptimization algorithm:\n ", x[["optimization_algorithm_used"]],
        ifelse(x[["optimization_algorithm_used"]] == "best_growth",
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
        ifelse(x[["optimization_algorithm_used"]][length(x[["optimization_algorithm_used"]])] == "best_growth",
          ifelse(x[["did_converge"]],
            "\n The last best_growth did converge.",
            "\n The last best_growth did not converge."
          ),
          ""
        ),
        sep = ""
      )
    }

    cat("\n\nNumber of log_likelihood calls:\n ",
      x[["number_of_log_likelihood_calls"]],
      "\n\nOptimization time:\n ", unclass(x[["full_optimization_time"]]),
      " ", attr(x[["full_optimization_time"]], "units"),
      sep = ""
    )

    if (all(x[["optimization_algorithm_used"]] == "Metropolis_Hastings")) {
      cat(paste0("\n\nAcceptance rate:\n ", x[["acceptance_rate"]]), sep = "")
    }
    if (x[["optimization_algorithm_used"]][length(x[["optimization_algorithm_used"]])] != "brute_force") {
      cat("\n\nLog_likelihood calls after the found permutation:\n ",
        x[["log_likelihood_calls_after_best"]],
        sep = ""
      )
    }
  }

  invisible(NULL)
}


# TODO(The base object printed as the gips_perm. The "size" attr omitted)
# TODO(`attr(object, "optimization_info")$visited_perms` can be printed as perms)
# TODO(Now I think it is good as it is with `utils::str.default`)
str.gips <- function(object, ...) {
  validate_gips(object)

  utils::str(object)
}
