#' Calculate Gamma function
#'
#' Theorem 8 from the paper, using the formula (19) from the paper
#'
#' @param perm an element of class "cycle"
#' @param lambda positive real number
#' 
#' @export
#'
#' @return Value of Gamma function

calculate_gamma_function <- function(perm, lambda){
  constants <- get_structure_constants(perm)
  r <- constants$r
  k <- constants$k
  d <- constants$d
  L <- constants$L
  
  if(lambda <= max((r-1)*d/(2*k))){
    return(Inf)  # the integral does not converge
  }
  
  A_Gamma <- sum(r*k*log(k))
  dim_gamma <- r + r*(r-1)*d/2
  B_Gamma <- sum(dim_gamma*log(k))/2
  
  gamma_omega <- sapply(1:L, function(i){calculate_gamma_omega(dim_gamma[i],
                                                               k[i]*lambda,
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
#' @return Value of Gamma function
calculate_gamma_omega <- function(lambda, dim_omega_i, r_i, d_i){
  if(lambda <= dim_omega_i/r_i - 1){
    return(Inf)  # the integral does not converge
  }
  
  prod(gamma((0:(r_i-1))*d_i/2 + lambda)) * (2*pi)^((dim_omega_i - r_i)/2)
}


#' G_function for goal_function
#' 
#' @param perm permutation of interest, object of `cycle` class
#' @param delta parameter of a method
#' @param structure_constants constants from `get_structure_constants` function
#' 
#' @return Product of `calculate_gamma_omega` from i to L. It is a product part of equation (27). For more information, see Issue #3 on `gips`' GitHub
#' 
#' @examples
#' perm_size <- 6
#' perm <- permutations::as.cycle(permutations::as.word(c(2,3,1,5,4,6)))
#' structure_constants <- get_structure_constants(perm, perm_size)
#' gips:::G_function(perm, structure_constants, 3)
G_function <- function(perm, structure_constants, delta=3){
  
  single_G_i <- sapply(1:structure_constants[['L']], function(i){
    lambda_i <- structure_constants[['r']][i] * (delta-2)/2 + structure_constants[['dim_omega']][i]/structure_constants[['r']][i]
    
    calculate_gamma_omega(lambda_i, structure_constants[['dim_omega']][i], structure_constants[['r']][i], structure_constants[['d']][i])
  })
  
  prod(single_G_i)
}












