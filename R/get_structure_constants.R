#' Get Structure Constants
#'
#' Finds constants necessary for internal calculations of integrals and
#' eventually the posteriori probability in [log_posteriori_of_gips()].
#'
#' Uses [Theorem 5 from references](https://arxiv.org/abs/2004.03503)
#' to calculate the constants.
#'
#' @param perm An object of a `gips_perm` class.
#'     It can also be of a `gips` class, but
#'     it will be interpreted as the underlying `gips_perm`.
#'
#' @returns Returns a list of 5 items:
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

  perm_size <- attr(perm, "size")
  l <- get_cycle_representatives_and_lengths(perm)
  representatives <- l[["representatives"]]
  cycle_lengths <- l[["cycle_lengths"]]
  perm_order <- ifelse(length(cycle_lengths) >= 2,
    numbers::mLCM(cycle_lengths),
    cycle_lengths
  )

  r <- calculate_r(cycle_lengths, perm_order)
  L <- length(r)
  d <- calculate_d(L, perm_order)

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

#' Get cycle representatives and lengths
#'
#' Essentially get iC, pC from paper
#'
#' @inheritParams get_structure_constants
#'
#' @returns A List with 2 items: `representatives` and `cycle_lengths`.
#'
#' @examples
#' perm <- gips_perm(permutations::as.cycle(permutations::as.word(c(4, 3, 6, 5, 1, 2))), 6)
#' get_cycle_representatives_and_lengths(perm)
#'
#' @noRd
get_cycle_representatives_and_lengths <- function(perm) {
  list(
    "representatives" = sapply(perm, function(v) v[1]),
    "cycle_lengths" = sapply(perm, length)
  )
}

#' Calculate structure constant r
#' Utilizes the sparsity for the `r_alfa` to speed up the computation.
#' The vector `r_alfa` is sorted by the value of the `alpha`.
#' Meaning the first element is the `r_alpha` for `alpha` = 0.
#'
#' @returns An integer vector. Structure constant r already without 0 elements.
#' @noRd
calculate_r <- function(cycle_lengths, perm_order) {
  M <- floor(perm_order / 2)
  if (M == 0) {
    # identity function
    return(length(cycle_lengths))
  }
  # for alpha in 0,1,...,floor(perm_order/2)
  # r_{alpha} = #{1:C such that alpha*p_c is a multiple of N}
  # alpha*p_c is a multiple of N iff alpha*p_c %% N == 0
  # Now, N is the Least Common Multiplier of all p_c
  # Which means, that N %% p_c == 0 for every p_c
  # In other words, N / p_c is an integer
  multiples <- round(perm_order / cycle_lengths) # the result of division should be an integer, but floats may interfere

  # Since N/p_c is an integer, alpha*p_c %% N == 0 iff alpha %% (N/p_c) == 0
  # In other words, alpha must be a multiple of (N/p_1, N/p_2,...,N/p_C)
  # (alpha = k*(N/p_1) or k*(N/p_2) or ... or k*(N/p_C) for some integer k (including 0))
  # However, alpha must be at most M, and a valid bound from above for integer k is max(p_1,...,p_C).
  max_order <- max(cycle_lengths)
  
  # Here we create all possible alpha values. 
  # The `multiples` corresponds to N/p_1,...,N/p_C, and `0:max_order` are possible integers k.
  # Use the outer product to get all pairwise multiplications
  alpha_matrix <- multiples %*% t(0:max_order)

  # sort is in ascending order, which means smallest alphas go to start.
  # The end result is as if we iterated over each alpha value (in ascending order),
  # and then deleted entries with 0s.
  possible_alphas <- unique(sort(alpha_matrix[alpha_matrix <= M]))

  # Recalculate the r_alpha vector using its definition directly.
  r_alfa <- sapply(possible_alphas, function(alpha) sum(alpha %% multiples == 0))
  as.double(r_alfa)
}

#' Calculate structure constant d.
#' Utilizing the structure of `r_alfa` vector.
#'
#' @noRd
calculate_d <- function(r_len, perm_order) {
  d <- rep(2, r_len)
  d[1] <- 1
  if (perm_order %% 2 == 0) {
    d[r_len] <- 1
  }
  d
}
