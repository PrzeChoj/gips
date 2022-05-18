#' Metropolis-Hastings algorithm
#'
#' Uses Metropolis-Hastings algorithm to find the permutation that maximizes the likelihood of observed data.
#'
#' @param U matrix, sum of outer products of data. `U` = sum(t(Z) %*% Z), where Z is observed data.
#' @param n_number number of data points that `U` is based on.
#' @param max_iter number of iterations for an algorithm to perform.
#' @param start starting permutation for the algorithm; an element of class "cycle". When NULL, identity permutation is taken.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of the same size as `U`. When NULL, the identity matrix is taken.
#' @param show_progress_bar boolean, indicate weather or not show the progress bar.
#'
#' @return list of 3 items: `acceptance_rate`, `goal_function_logvalues`,
#' `points`, `found_point`, `found_point_function_logvalue`
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
#' mh <- MH(U=U, n_number=10, max_iter=100, start=start,
#'          delta=3, D_matrix=diag(nrow=dim(U)[1]))
MH <- function(U, n_number, max_iter, start=NULL,
               delta=3, D_matrix=NULL, show_progress_bar=TRUE){
  if(is.null(start)){
    start <- permutations::id
  }
  stopifnot(dim(U)[2] == dim(U)[1])
  perm_size <- dim(U)[1]
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }

  acceptance <- rep(FALSE, max_iter)
  goal_function_logvalues <- rep(0, max_iter)
  points <- list()
  points[[1]] <- start

  if(show_progress_bar)
    progressBar = utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  goal_function_logvalues[1] <- goal_function(points[[1]],
                                              n_number, U,
                                              delta=delta, D_matrix=D_matrix)

  found_point <- start
  found_point_function_logvalue <- goal_function_logvalues[1]

  U2 <- stats::runif(max_iter, min = 0, max = 1)

  for (i in 1:(max_iter-1)){
    if(show_progress_bar)
      utils::setTxtProgressBar(progressBar, i)

    e <- runif_transposition(perm_size)
    perm_proposal <- permutations::as.cycle(points[[i]] * e)

    goal_function_perm_proposal <- goal_function(perm_proposal,
                                                 n_number, U,
                                                 delta=delta, D_matrix=D_matrix)
    if(is.nan(goal_function_perm_proposal) | is.infinite(goal_function_perm_proposal)){
      #browser()  # needs further investigation. See ISSUE#5
      warning("gips is yet unable to process this U matrix. See ISSUE#5 for more information")
      return(list("acceptance_rate"=mean(acceptance),
                  "goal_function_values"=goal_function_values,
                  "points"=points,
                  "found_point"=found_point,
                  "found_point_function_value"=found_point_function_value))
    }

    # if goal_function_perm_proposal > goal_function_values[i], then it is true
    if(U2[i] < exp(goal_function_perm_proposal-goal_function_logvalues[i])){ # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2)
      points[[i+1]] <- perm_proposal
      goal_function_logvalues[i+1] <- goal_function_perm_proposal
      acceptance[i] <- TRUE

      if(found_point_function_logvalue < goal_function_logvalues[i+1]){
        found_point_function_logvalue <- goal_function_logvalues[i+1]
        found_point <- points[[i+1]]
      }
    }
    else{
      points[[i+1]] = points[[i]]
      goal_function_logvalues[i+1] <- goal_function_logvalues[i]
    }
  }

  if(show_progress_bar)
    close(progressBar)

  list("acceptance_rate"=mean(acceptance),
       "goal_function_logvalues"=goal_function_logvalues,
       "points"=points,
       "found_point"=found_point,
       "found_point_function_logvalue"=found_point_function_logvalue)
}



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
  stopifnot(dim(U)[1] == dim(U)[2])
  perm_size <- dim(U)[1]

  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)  # identity matrix
  }

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



