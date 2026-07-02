#' Calculate Gamma function
#'
#' Calculates the integral defined in
#' [Definition 11 from references](https://arxiv.org/abs/2004.03503).
#' Implements
#' [Theorem 8 from references](https://arxiv.org/abs/2004.03503)
#' and uses the
#' [formula (19) from references](https://arxiv.org/abs/2004.03503).
#'
#' @inheritParams get_structure_constants
#' @param lambda A positive real number.
#'
#' @returns The value of the Gamma function of the colored cone
#'     (for the definition of the colored cone, see the **Basic definitions**
#'     section in `vignette("Theory", package = "gips")` or in its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html)).
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [arXiv link](https://arxiv.org/abs/2004.03503);
#' \doi{10.1214/22-AOS2174}
#'
#' @seealso
#' * [get_structure_constants()] - The function useful inside
#'     the `calculate_gamma_function()`.
#' * [log_posteriori_of_gips()] - The function that uses
#'     the values of the gamma function.
#' * `vignette("Theory", package = "gips")` or its
#'     [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html) -
#'     A place to learn more about
#'     the math behind the `gips` package.
#'
#' @examples
#' id_perm <- gips_perm("()", 2)
#' calculate_gamma_function(id_perm, 0.5001) # 44417.54
#' calculate_gamma_function(id_perm, 0.50000001) # 444288280
#' calculate_gamma_function(id_perm, 0.500000000001) # 4.442981e+12
#'
#' oldw <- getOption("warn")
#' options(warn = -1)
#' calculate_gamma_function(id_perm, 0.5) # Inf
#' # Integral diverges; returns Inf and warning
#' options(warn = oldw)
#'
#' @export
calculate_gamma_function <- function(perm, lambda) {
  constants <- get_structure_constants(perm)
  r <- constants[["r"]]
  k <- constants[["k"]]
  d <- constants[["d"]]
  L <- constants[["L"]]
  dim_gamma <- constants[["dim_omega"]]

  if (lambda <= max((r - 1) * d / (2 * k))) {
    rlang::warn(c("Gamma integral is divergent for the given permutation and lambda value.",
      "i" = paste0(
        "Gamma(perm = ", as.character(perm),
        ", lambda = ", lambda, ") = Inf."
      )
    ))
    return(Inf) # the integral does not converge
  }

  A_Gamma <- sum(r * k * log(k))
  B_Gamma <- sum(dim_gamma * log(k)) / 2

  gamma_omega <- sapply(1:L, function(i) {
    calculate_log_gamma_omega(
      k[i] * lambda,
      dim_gamma[i],
      r[i],
      d[i]
    )
  })

  exp(-A_Gamma * lambda + B_Gamma + sum(gamma_omega))
}


#' Calculate the logarithm of a single Gamma omega function
#'
#' Using the formula (12) from the paper
#'
#' @inheritParams calculate_gamma_function
#' @param dim_omega_i Single element from `get_structure_constants`.
#' @param r_i Single element from `get_structure_constants`.
#' @param d_i Single element from `get_structure_constants`.
#'
#' @returns Logarithm of the value of Gamma function.
#'
#' @noRd
calculate_log_gamma_omega <- function(lambda, dim_omega_i, r_i, d_i) {
  if ((lambda + 1) * r_i <= dim_omega_i) { # equivalent to lambda <= dim_omega_i / r_i - 1
    rlang::warn(c("Gamma integral is divergent for the given lambda value and structure constants.",
      "i" = paste0(
        "Gamma(lambda = ", lambda,
        ", dim_omega_i = ", dim_omega_i,
        ", r_i = ", r_i,
        ", d_i = ", d_i,
        ") = Inf."
      )
    ))
    return(Inf) # the integral does not converge
  }

  sum(lgamma((0:(-(r_i - 1))) * d_i / 2 + lambda)) + (dim_omega_i - r_i) / 2 * log(2 * pi)
}


#' Calculate the G part of `log_posteriori_of_gips()`
#'
#' @param structure_constants Constants from `get_structure_constants()`.
#' @param delta Parameter of a method.
#' @param number_of_observations Number of observations.
#'
#' @returns Difference `G(delta + number_of_observations) - G(delta)`.
#' @noRd
calculate_G_part <- function(structure_constants, delta, number_of_observations) {
  L <- structure_constants[["L"]]
  k <- structure_constants[["k"]]
  r <- structure_constants[["r"]]
  d <- structure_constants[["d"]]
  dim_omega <- structure_constants[["dim_omega"]]

  lambda_prior <- k * (delta - 2) / 2 + dim_omega / r
  lambda_posterior <- lambda_prior + k * number_of_observations / 2

  divergent <- (lambda_prior + 1) * r <= dim_omega |
    (lambda_posterior + 1) * r <= dim_omega
  if (any(divergent)) {
    stop(
      paste0(
        "Gamma integral divergence detected while calculating the G part. ",
        "This should not happen for valid `gips` inputs. ",
        "This follows from the facts that `delta > 1`, `k` is `1` or `2` and `r >= 1`. ",
        "If You see this, please open an ISSUE on GitHub to let us know."
      ),
      call. = FALSE
    )
  }

  out <- 0
  for (i in seq_len(L)) {
    offsets <- (0:(-(r[i] - 1))) * d[i] / 2
    gamma_values <- lgamma(c(
      lambda_posterior[i] + offsets,
      lambda_prior[i] + offsets
    ))
    block_indices <- seq_len(r[i])

    out <- out +
      sum(gamma_values[block_indices]) -
      sum(gamma_values[r[i] + block_indices])
  }

  out
}
