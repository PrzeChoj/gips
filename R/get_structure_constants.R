#' Get Structure Constants
#'
#' Finds constants necessary for internal calculations of integrals and
#' eventually the posterior probability in [log_posteriori_of_gips()].
#'
#' Uses [Theorem 5 from references](https://arxiv.org/abs/2004.03503)
#' to calculate the constants.
#'
#' @param perm An object of a `gips_perm` class.
#'     It can also be of a `gips` class, but
#'     it will be interpreted as the underlying `gips_perm`.
#'
#' @returns A list with 5 items:
#'     `r`, `d`, `k`, `L`, `dim_omega` - vectors of constants from
#'     [Theorem 1 from references](https://arxiv.org/abs/2004.03503)
#'     and the beginning of
#'     [section 3.1. from references](https://arxiv.org/abs/2004.03503).
#' @export
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [arXiv link](https://arxiv.org/abs/2004.03503);
#' \doi{10.1214/22-AOS2174}
#'
#' @seealso
#' * [calculate_gamma_function()], [log_posteriori_of_gips()] - The functions
#'     that rely heavily on `get_structure_constants()`.
#'
#' @examples
#' perm <- gips_perm("(1)(2)(3)(4,5)", 5)
#' get_structure_constants(perm)
get_structure_constants <- function(perm) {
  if (inherits(perm, "gips")) {
    validate_gips(perm)
    perm <- perm[[1]]
  }
  if (!(inherits(perm, "gips_perm"))) {
    wrong_argument_abort(
      i = "`perm` must be of a `gips_perm` class.",
      x = paste0(
        "You provided `perm` with `class(perm) == (",
        paste(class(perm), collapse = ", "), ")`."
      )
    )
  }

  validate_gips_perm(perm)
  get_structure_constants_no_checks(perm)
}

# Internal calculations for valid permutations; public inputs must be validated first.
get_structure_constants_no_checks <- function(perm) {
  cycle_lengths <- lengths(perm)

  r <- calculate_r(cycle_lengths)
  L <- length(r)
  d <- calculate_d(L, cycle_lengths)

  k <- d
  dim_omega <- r + r * (r - 1) * d / 2

  list(
    "r" = r,
    "d" = d,
    "k" = k,
    "L" = L,
    "dim_omega" = dim_omega
  )
}

# Internal Euclidean algorithm for a vector of nonnegative integers a and a scalar b > 0.
# Inputs come from validated cycle lengths and generated frequency numerators.
gcd_nonnegative <- function(a, b) {
  b <- rep.int(b, length(a))
  active <- a > 0
  while (any(active)) {
    remainder <- b[active] %% a[active]
    b[active] <- a[active]
    a[active] <- remainder
    active <- a > 0
  }
  b
}

#' Get cycle frequencies as reduced fractions
#'
#' @param cycle_length A positive integer cycle length from a validated permutation.
#' @returns A list with `numerator` and `denominator` vectors for frequencies
#'     between 0 and 1/2, in increasing order.
#' @noRd
get_cycle_frequencies <- function(cycle_length) {
  numerators <- seq.int(0, floor(cycle_length / 2))
  divisors <- gcd_nonnegative(numerators, cycle_length)

  list(
    numerator = numerators / divisors,
    denominator = cycle_length / divisors
  )
}

#' Calculate structure constant r
#' Count cycles containing each distinct frequency, in increasing frequency order.
#' Frequencies are identified by their reduced numerator-denominator pairs.
#'
#' @returns A numeric vector of positive integer multiplicities.
#' @noRd
calculate_r <- function(cycle_lengths) {
  if (length(cycle_lengths) == 0L) {
    rlang::abort("Cannot calculate structure constants for a permutation without cycles.")
  }
  cycle_frequencies <- lapply(cycle_lengths, get_cycle_frequencies)
  numerators <- unlist(lapply(cycle_frequencies, `[[`, "numerator"), use.names = FALSE)
  denominators <- unlist(lapply(cycle_frequencies, `[[`, "denominator"), use.names = FALSE)

  frequencies <- numerators / denominators
  frequency_order <- order(frequencies)
  frequencies <- frequencies[frequency_order]
  numerators <- numerators[frequency_order]
  denominators <- denominators[frequency_order]

  # Group by exact reduced fractions, not by their floating-point values.
  is_new_frequency <- c(TRUE, diff(numerators) != 0 | diff(denominators) != 0)
  if (anyDuplicated(frequencies[is_new_frequency]) > 0L) {
    rlang::abort(c(
      "Cannot reliably order cycle frequencies.",
      "x" = "Distinct fractions have the same floating-point value."
    ), class = "indistinguishable_cycle_frequencies")
  }

  as.double(diff(c(which(is_new_frequency), length(frequencies) + 1)))
}

#' Calculate structure constant d.
#' The first block has frequency 0; the last has frequency 1/2 if any cycle is even.
#'
#' @noRd
calculate_d <- function(r_len, cycle_lengths) {
  d <- rep(2, r_len)
  d[1] <- 1
  # The LCM is even exactly when at least one cycle has even length.
  if (any(cycle_lengths %% 2 == 0)) {
    d[r_len] <- 1
  }
  d
}
