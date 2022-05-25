#' goal function for MH
#'
#' Calculate the logarithm of function proportional to a posteriori distribution, according to equation (33) and (27). If `Inf` or `NaN` is reached, produces a warning.
#'
#' @export
#'
#' @param perm_proposal permutation of interest.
#' @param n_number number of random variables that `U` is based on.
#' @param U matrix that the projection of is wanted.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken.
#'
#' @examples
#' c <- permutations::as.cycle(permutations::as.word(c(2,1)))
#' U1 <- matrix(c(1,0.5,0.5,2), nrow=2,byrow = TRUE)
#' goal_function(c, 100, U1)
goal_function <- function(perm_proposal, n_number, U, delta=3, D_matrix=NULL){
  stopifnot(permutations::is.cycle(perm_proposal),
            is.matrix(U),
            dim(U)[1] == dim(U)[2])
  perm_size <- dim(U)[1]

  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)  # identity matrix
  }
  stopifnot(is.matrix(D_matrix),
            dim(D_matrix)[1] == dim(D_matrix)[2])

  structure_constants <- get_structure_constants(perm_proposal, perm_size)

  # Ac_part
  Ac <- sum(structure_constants[['r']]*structure_constants[['k']]*log(structure_constants[['k']]))  # (20)
  Ac_part <- (-n_number/2*Ac)

  # G_part and phi_part
  G_part <- G_function(perm_proposal, structure_constants, delta + n_number) -
      G_function(perm_proposal, structure_constants, delta)

  # phi_part
  phi_part <- calculate_phi_part(perm_proposal, perm_size, n_number, U, delta,
                                 D_matrix, structure_constants)

  out <- Ac_part + G_part + phi_part

  if(is.infinite(out)){
    warning("Infinite value of a goal function")
  }
  if(is.nan(out)){
    warning("NaN value of a goal function")
  }

  out
}

#' Uniformly random transposition of perm_size elements
#'
#' @param perm_size size from which take transpositions
runif_transposition <- function(perm_size){
  permutations::as.cycle(sample(perm_size, 2, replace=FALSE))
}

#' Calculate log phi_part of goal_function
#'
#' @param structure_constants output of get_structure_constants(perm_proposal, perm_size)
#' Rest of params as in goal_function
#'
#' @noRd
calculate_phi_part <- function(perm_proposal, perm_size, n_number, U, delta,
                               D_matrix, structure_constants){

    # projection of matrices on perm_proposal
    Dc <- project_matrix(D_matrix, perm_proposal, perm_size)
    Uc <- project_matrix(U, perm_proposal, perm_size)

    # divide by 2 - refer to newest version of the paper
    Dc <- Dc / 2
    Uc <- Uc / 2

    # diagonalization
    # TODO add basis argument? ISSUE#6
    diagonalising_matrix <- prepare_orthogonal_matrix(perm_proposal,
                                                      perm_size)
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
    DcUc_exponent <- -(n_number+delta-2)/2 - structure_constants[['dim_omega']] /
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



