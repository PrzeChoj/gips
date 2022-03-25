#' Calculate Gamma function
#'
#' Theorem 8 from the paper, using the formula (19) from the paper
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
                                                               r[i], d[i])})
  
  exp(-A_Gamma * lambda + B_Gamma) * prod(gamma_omega)
}


#' Calculate single Gamma omega function
#'
#' Using the formula (12) from the paper
#'
#' @return Value of Gamma function

calculate_gamma_omega <- function(dim_omega, lambda, r, d){
  if(lambda <= dim_omega/r - 1){
    return(Inf)  # the integral does not converge
  }
  
  prod(gamma(0:(r-1)*d/2 + lambda)) * (2*pi)^((dim_omega - r)/2)
}





