#' Metropolis-Hastings algorithm
#'
#' Perform the algorithm on the permutations
#'
#' @param U Matrix that the projection of is wanted
#' @param start Start of the algorithm; an element of class "cycle"
#' @param max_iter Number of iterations
#' @param perm_size The dimension of interest. When NULL, the size of U is taken
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken
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
#'          delta=3, D=diag(nrow = perm_size))

MH <- function(U, start, max_iter, perm_size=NULL, delta=3, D=NULL){
  if(is.null(perm_size)){
    perm_size <- dim(U)[1]
  }
  stopifnot(perm_size == dim(U)[1])
  if(is.null(D)){
    D <- diag(nrow = perm_size)
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
#' @param D hyper-parameter of a Bayesian model. Square matrix of size `perm_size`. When NULL, the identity matrix is taken
goal_function <- function(perm_proposal, perm_size, n, U, delta=3, D=NULL){
  if(is.null(D)){
    D <- diag(nrow = perm_size)  # identity matrix
  }
  
  structure_constants <- get_structure_constants(perm_proposal, perm_size)
  
  # exp_part
  Ac <- sum(r*k*log(k))  # (20)
  exp_part <- exp(-n/2*Ac)
  
  # G_part
  # TODO
  G_part <- 7  # G(perm_proposal, delta + n) / G(perm_proposal, delta)
  
  # projection of matrices on perm_proposal
  # TODO
  Dc <- D
  Uc <- U
  
  # det_part
  det_part <- det(Dc+Uc)^(-(n+delta-2)/2) * det(Dc)^((delta-2)/2)  # when D = I identity matrix, then Dc=I, then the second part is 1
  
  # phi_part
  # TODO
  phi_part <- 4  # phi(perm_proposal, D + U) / phi(perm_proposal, D)
  
  exp_part * G_part * det_part * phi_part
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