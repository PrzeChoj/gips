#' Metropolis-Hastings algorithm
#'
#' Perform the algorithm on the permutations
#'
#' @param U Matrix that the projection of is wanted
#' @param start Start of the algorithm; an element of class "cycle"
#' @param max_iter Number of iterations
#' @param perm_size The dimension of interest. When NULL, the size of U is taken
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken
#'
#' @return list of 3 items: `acceptance_rate`, `goal_function_values`, `points`
#' @export
#'
#' @examples
#' perm_size <- 6
#' mu <- numeric(perm_size)
#' # sigma is a permutation (1,2,3,4,5,6)
#' sigma <- matrix(data = c(1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
#'                          0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
#'                          0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
#'                          0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
#'                          0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
#'                          0.8, 0.6, 0.4, 0.6, 0.8, 1.0),
#'                 nrow=perm_size, byrow = TRUE)
#' N <- 6
#' Z <- MASS::mvrnorm(N, mu = mu, Sigma = sigma )
#' U <- (t(Z) %*% Z)/N
#' start <- permutations::id
#' mh <- MH(U=U, start=start, max_iter=100, perm_size=perm_size,
#'          delta=3, D_matrix=diag(nrow = perm_size))
MH <- function(U, start, max_iter, perm_size=NULL, delta=3, D_matrix=NULL){
  if(is.null(perm_size)){
    perm_size <- dim(U)[1]
  }
  stopifnot(perm_size == dim(U)[1])
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }

  acceptance <- rep(FALSE, max_iter)
  goal_function_values <- rep(0, max_iter)
  points <- list()
  points[[1]] <- start

  goal_function_values[1] <- test_goal_function(points[[1]])
  #goal_function_values[1] <- goal_function(points[[1]])  # TODO(goal_function is work in progress)

  U2 <- stats::runif(max_iter, min = 0, max = 1)

  for (i in 1:(max_iter-1)){
    e <- runif_transposition(perm_size)
    q <- points[[i]] * e

    goal_function_q <- test_goal_function(q)
    #goal_function_q <- goal_function(q)  # TODO(goal_function is work in progress)

    # if goal_function_q > goal_function[i], then it is true
    if(U2[i] < goal_function_q/goal_function_values[i]){ # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2)
      points[[i+1]] <- q
      goal_function_values[i+1] <- goal_function_q
      acceptance[i] <- TRUE
    }
    else{
      points[[i+1]] = points[[i]]
      goal_function_values[i+1] <- goal_function_values[i]
    }
  }

  list("acceptance_rate"=mean(acceptance), "goal_function_values"=goal_function_values, "points"=points)
}



#' goal_function for MH
#'
#' @export
#'
#' @param perm_proposal Permutation of interest
#' @param perm_size size of a permutation
#' @param n Size of a sample
#' @param U Matrix that the projection of is wanted
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken
goal_function <- function(perm_proposal, perm_size, n, U, delta=3, D_matrix=NULL){
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)  # identity matrix
  }

  structure_constants <- get_structure_constants(perm_proposal, perm_size)

  # exp_part
  Ac <- sum(structure_constants[['r']]*structure_constants[['k']]*log(structure_constants[['k']]))  # (20)
  exp_part <- exp(-n/2*Ac)

  # G_part
  G_part <- G_function(perm_proposal, structure_constants, delta + n) / G_function(perm_proposal, structure_constants, delta)

  # projection of matrices on perm_proposal
  Dc <- project_matrix(D_matrix, perm_proposal, perm_size)
  Uc <- project_matrix(U, perm_proposal, perm_size)

  # diagonalisation
  # TODO add basis argument
  diagonalising_matrix <- prepare_orthogonal_matrix(perm_proposal,
                                                    perm_size)
  Dc_diagonalised <- t(diagonalising_matrix) %*% Dc %*% diagonalising_matrix
  DcUc_diagonalised <- t(diagonalising_matrix) %*% (Uc+Dc) %*% diagonalising_matrix

  # det_phi_part
  block_ends <- cumsum(structure_constants[['r']] * structure_constants[['d']])
  Dc_block_dets <- calculate_determinants_of_block_matrices(Dc_diagonalised,
                                                            block_ends)
  DcUc_block_dets <- calculate_determinants_of_block_matrices(DcUc_diagonalised,
                                                              block_ends)
  Dc_exponent <- -(n+delta-2)/2 - structure_constants[['dim_omega']] /
      (structure_constants[['r']] * structure_constants[['k']])
  DcUc_exponent <- (delta-2)/2 + structure_constants[['dim_omega']] /
      (structure_constants[['r']] * structure_constants[['k']])

  det_phi_part <- prod(Dc_block_dets ^ Dc_exponent * DcUc_block_dets ^ DcUc_exponent)

  exp_part * G_part * det_phi_part
}

#' example goal function
#' Used just for testing
#'
#' @param perm_proposal permutation of interest
test_goal_function <- function(perm_proposal){
  permutations::permorder(perm_proposal) + 1
}

#' Uniformly random transposition of perm_size elements
#'
#' @param perm_size size from which take transpositions
runif_transposition <- function(perm_size){
  permutations::as.cycle(sample(perm_size, 2, replace=FALSE))
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


