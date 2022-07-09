#' The log likelihood that the CoV matrix is invariant under permutation.
#' 
#' Calculate the logarithm of function proportional to a posteriori distribution, according to equation (33) and (27). If `Inf` or `NaN` is reached, produces a warning.
#' Keep in mind this is not the log likelihood per se, but rather the funcition that is proportional to the log likelihood.
#' It is the goal function for optimization algorithms.
#'
#' @export
#'
#' @param perm_proposal permutation of interest.
#' @param number_of_observations number of random variables that `S` is based on.
#' @param S matrix, estimated covariance matrix. When Z is observed data: `S = sum(t(Z) %*% Z)/number_of_observations`, if one know the theoretical mean is 0; # TODO(What if one have to estimate the theoretical mean with the empirical mean)
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken.
#'
#' @examples
#' c <- permutations::as.cycle(permutations::as.word(c(2,1)))
#' S1 <- matrix(c(1,0.5,0.5,2), nrow=2, byrow = TRUE)
#' log_likelihood_of_perm(c, 100, S1)
log_likelihood_of_perm <- function(perm_proposal, number_of_observations, S,
                                   delta=3, D_matrix=NULL){
  stopifnot(permutations::is.cycle(perm_proposal) || inherits(perm_proposal, 'gips_perm'),
            is.matrix(S),
            dim(S)[1] == dim(S)[2])
  
  U <- S * number_of_observations  # in the paper there is U everywhere in stead of S, so it is easier to use U matrix in the code
  perm_size <- dim(S)[1]
  
  if(!inherits(perm_proposal, 'gips_perm'))
    perm_proposal <- gips_perm(perm_proposal, perm_size)


  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)  # identity matrix
  }
  stopifnot(is.matrix(D_matrix),
            dim(D_matrix)[1] == dim(D_matrix)[2])

  structure_constants <- get_structure_constants(perm_proposal)

  # Ac_part
  Ac <- sum(structure_constants[['r']]*structure_constants[['k']]*log(structure_constants[['k']]))  # (20)
  Ac_part <- (-number_of_observations/2*Ac)

  # G_part and phi_part
  G_part <- G_function(structure_constants, delta + number_of_observations) -
    G_function(structure_constants, delta)

  # phi_part
  phi_part <- calculate_phi_part(perm_proposal, number_of_observations, U, delta,
                                 D_matrix, structure_constants)

  out <- Ac_part + G_part + phi_part

  if(is.infinite(out)){
    warning("Infinite value of a likelihood")
  }
  if(is.nan(out)){
    warning("NaN value of a likelihood")
  }

  out
}

#' Uniformly random transposition of perm_size elements
#'
#' @param perm_size size from which take transpositions
runif_transposition <- function(perm_size){
  sample(perm_size, 2, replace=FALSE)
}

#' Calculate log phi_part of log_likelihood_of_perm
#'
#' @param structure_constants output of get_structure_constants(perm_proposal, perm_size)
#' Rest of params as in log_likelihood_of_perm
#'
#' @noRd
calculate_phi_part <- function(perm_proposal, number_of_observations, U, delta,
                               D_matrix, structure_constants){

  # projection of matrices on perm_proposal
  equal_indices <- get_equal_indices_by_perm(perm_proposal)
  Dc <- project_matrix(D_matrix, perm_proposal,
                       precomputed_equal_indices=equal_indices)
  Uc <- project_matrix(U, perm_proposal,
                       precomputed_equal_indices=equal_indices)

  # divide by 2 - refer to newest version of the paper
  Dc <- Dc / 2
  Uc <- Uc / 2

  # diagonalization
  # TODO add basis argument? ISSUE#6
  diagonalising_matrix <- prepare_orthogonal_matrix(perm_proposal)
  Dc_diagonalised <- t(diagonalising_matrix) %*% Dc %*% diagonalising_matrix
  DcUc_diagonalised <- t(diagonalising_matrix) %*% (Uc+Dc) %*% diagonalising_matrix

  # block part
  block_ends <- cumsum(structure_constants[['r']] * structure_constants[['d']])
  Dc_block_dets <- calculate_determinants_of_block_matrices(Dc_diagonalised,
                                                            block_ends)
  DcUc_block_dets <- calculate_determinants_of_block_matrices(DcUc_diagonalised,
                                                              block_ends)
  Dc_exponent <- (delta-2)/2 + structure_constants[['dim_omega']] /
    (structure_constants[['r']] * structure_constants[['k']])
  DcUc_exponent <- -(number_of_observations+delta-2)/2 - structure_constants[['dim_omega']] /
    (structure_constants[['r']] * structure_constants[['k']])

  out <- sum(log(Dc_block_dets) * Dc_exponent + log(DcUc_block_dets) * DcUc_exponent)

  out
}

#' Calculate determinants of matrices from block decomposition
#'
#' Block decomposition 1 from paper
#'
#' @param diagonalized_matrix middle matrix from decomposition 1
#' @param block_ends indices of last columns of block matrices.
#' Last element equals size of matrix
#'
#' @return numeric vector
#' @noRd
calculate_determinants_of_block_matrices <- function(diagonalised_matrix,
                                                     block_ends){
  block_starts <- c(0, block_ends[-length(block_ends)]+1)
  sapply(1:length(block_starts), function(i){
    slice <- block_starts[i]:block_ends[i]
    block_matrix <- diagonalised_matrix[slice, slice, drop=FALSE]
    det(block_matrix)
  })
}
