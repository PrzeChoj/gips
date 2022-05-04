#' Metropolis-Hastings algorithm
#'
#' Perform the algorithm on the permutations
#'
#' @param U matrix that the projection of is wanted.
#' @param n_number number of random variables that `U` is based on.
#' @param max_iter number of iterations.
#' @param start start of the algorithm; an element of class "cycle". When NULL, identity permutation is taken.
#' @param perm_size the dimension of interest. When NULL, the size of U is taken.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken.
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
#' n_number <- 6
#' Z <- MASS::mvrnorm(n_number, mu = mu, Sigma = sigma)
#' U <- (t(Z) %*% Z)/n_number
#' start <- permutations::id
#' mh <- MH(U=U, n_number=10, max_iter=100, start=start, perm_size=perm_size,
#'          delta=3, D_matrix=diag(nrow = perm_size))
MH <- function(U, n_number, max_iter, start=NULL, perm_size=NULL, delta=3, D_matrix=NULL){
  if(is.null(start)){
    start <- permutations::id
  }
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

  #goal_function_values[1] <- test_goal_function(points[[1]])
  goal_function_values[1] <- goal_function(points[[1]],
                                           perm_size, n_number, U,
                                           delta=3, D_matrix=D_matrix)

  U2 <- stats::runif(max_iter, min = 0, max = 1)

  for (i in 1:(max_iter-1)){
    if(i%%100 == 0){print(i)}
    e <- runif_transposition(perm_size)
    perm_proposal <- points[[i]] * e

    #goal_function_perm_proposal <- test_goal_function(perm_proposal)
    goal_function_perm_proposal <- goal_function(perm_proposal,
                                                 perm_size, n_number, U,
                                                 delta=3, D_matrix=D_matrix)
    
    # if goal_function_perm_proposal > goal_function[i], then it is true
    if(U2[i] < goal_function_perm_proposal/goal_function_values[i]){ # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2)
      points[[i+1]] <- perm_proposal
      goal_function_values[i+1] <- goal_function_perm_proposal
      acceptance[i] <- TRUE
    }
    else{
      points[[i+1]] = points[[i]]
      goal_function_values[i+1] <- goal_function_values[i]
    }
  }

  list("acceptance_rate"=mean(acceptance), "goal_function_values"=goal_function_values, "points"=points)
}



#' goal function for MH
#'
#' If infinite value is reached, produces a warning
#'
#' @export
#'
#' @param perm_proposal permutation of interest.
#' @param perm_size size of a permutation.
#' @param n_number number of random variables that `U` is based on.
#' @param U matrix that the projection of is wanted.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken.
#'
#' @examples
#' c <- permutations::as.cycle(permutations::as.word(c(2,1)))
#' U1 <- matrix(c(1,0.5,0.5,2), nrow=2,byrow = TRUE)
#' goal_function(c, 2, 100, U1)
goal_function <- function(perm_proposal, perm_size, n_number, U, delta=3, D_matrix=NULL){
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)  # identity matrix
  }

  structure_constants <- get_structure_constants(perm_proposal, perm_size)

  # exp_part
  Ac <- sum(structure_constants[['r']]*structure_constants[['k']]*log(structure_constants[['k']]))  # (20)
  exp_part <- exp(-n_number/2*Ac)

  # G_part
  G_part <- G_function(perm_proposal, structure_constants, delta + n_number) /
      G_function(perm_proposal, structure_constants, delta)

  phi_part <- calculate_phi_part(perm_proposal, perm_size, n_number, U, delta,
                                 D_matrix, structure_constants)

  out <- exp_part * G_part * phi_part

  if(is.infinite(out)){
    warning("Infinite value of a goal function")
  }

  out
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

#' Calculate phi_part of goal_function
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

    # diagonalisation
    # TODO add basis argument?
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
    
    out <- prod(Dc_block_dets ^ Dc_exponent * DcUc_block_dets ^ DcUc_exponent)
    
    if(is.nan(out)){ # TODO This is temporary solution, see issue#5
      warning("NaN value of a calculate_phi_part function")
      out <- 0
    }
    
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



