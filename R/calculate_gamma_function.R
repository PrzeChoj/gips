#' Calculate Gamma function
#'
#' Theorem 8 from the paper, using the formula (19) from the paper
#'
#' @param perm an element of class "gips_perm"
#' @param lambda positive real number
#'
#' @export
#'
#' @return Value of Gamma function
#'
#' @examples
#' id_perm <- gips_perm(permutations::id, 2)
#' calculate_gamma_function(id_perm, 0.5001)
#' calculate_gamma_function(id_perm, 0.50000001)
#' calculate_gamma_function(id_perm, 0.500000000001)
#' #calculate_gamma_function(id_perm, 0.5) # integral diverges; returns Inf and warning
calculate_gamma_function <- function(perm, lambda){
  constants <- get_structure_constants(perm)
  r <- constants[['r']]
  k <- constants[['k']]
  d <- constants[['d']]
  L <- constants[['L']]
  dim_gamma <- constants[['dim_omega']]

  if(lambda <= max((r-1)*d/(2*k))){
    warning("Gamma integral does not convarge for a given lambda value")
    return(Inf)  # the integral does not converge
  }

  A_Gamma <- sum(r*k*log(k))
  B_Gamma <- sum(dim_gamma*log(k))/2

  gamma_omega <- sapply(1:L, function(i){calculate_gamma_omega(k[i]*lambda,
                                                               dim_gamma[i],
                                                               r[i],
                                                               d[i])})

  exp(-A_Gamma * lambda + B_Gamma) * prod(gamma_omega)
}


#' Calculate single Gamma omega function
#'
#' Using the formula (12) from the paper
#'
#' @param lambda positive real number
#' @param dim_omega_i single element from `get_structure_constants`
#' @param r_i single element  from `get_structure_constants`
#' @param d_i single element  from `get_structure_constants`
#'
#' @return Logarithm of value of Gamma function
calculate_gamma_omega <- function(lambda, dim_omega_i, r_i, d_i){
  if(lambda <= dim_omega_i/r_i - 1){
    warning("Gamma omega integral does not convarge for a given lambda value")
    return(Inf)  # the integral does not converge
  }

  #prod(gamma((0:(-(r_i-1)))*d_i/2 + lambda)) * (2*pi)^((dim_omega_i - r_i)/2)
  sum(lgamma((0:(-(r_i-1)))*d_i/2 + lambda)) + (dim_omega_i - r_i)/2 * log(2*pi)
}


#' G_function for `log_likelihood_of_perm()`
#'
#' @param delta parameter of a method
#' @param structure_constants constants from `get_structure_constants` function
#'
#' @return Sum of logarithms of elements of `calculate_gamma_omega` from i to L.
#' It is log of a product part of equation (27). For more information, see Issue #3 on `gips`' GitHub
#'
#' @examples
#' perm_size <- 6
#' perm <- permutations::as.cycle(permutations::as.word(c(2,3,1,5,4,6)))
#' gips_perm <- gips_perm(perm, perm_size)
#' structure_constants <- get_structure_constants(gips_perm)
#' gips:::G_function(structure_constants, 3)

G_function <- function(structure_constants, delta=3){

  single_G_i <- sapply(1:structure_constants[['L']], function(i){
    lambda_i <- structure_constants[['k']][i] * (delta-2)/2 + structure_constants[['dim_omega']][i]/structure_constants[['r']][i]

    calculate_gamma_omega(lambda_i, structure_constants[['dim_omega']][i], structure_constants[['r']][i], structure_constants[['d']][i])
  })

  sum(single_G_i)
}












