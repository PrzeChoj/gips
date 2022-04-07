#' Metropolis-Hastings algorithm
#'
#' Perform the algorithm on the permutations
#'
#' @param start Start of the algorithm; an element of class "cycle"
#' @param perm_size The dimension of interest
#' @param max_iter Number of iterations
#'
#' @return list of 3 items: `acceptance_rate`, `goal_function_values`, `points`
#' @export
#' 
#' @examples
#' start <- permutations::id
#' mh <- MH(start = start, 8, 100)

MH <- function(start, perm_size, max_iter){
  acceptance <- rep(FALSE, max_iter)
  goal_function_values <- rep(0, max_iter)
  points <- list()
  points[[1]] <- start
  goal_function_values[1] <- goal_function(points[[1]])
  
  U2 <- stats::runif(max_iter, min = 0, max = 1)
  
  for (i in 1:(max_iter-1)){
    e <- runif_transposition(perm_size)
    q <- points[[i]] * e
      
    goal_function_q <- goal_function(q)
    
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






goal_function <- function(perm){
  permutations::permorder(perm) + 1 # example function
}

#' Uniformly random transposition of perm_size elements
#' 
#' @param perm_size size from which take transpositions
runif_transposition <- function(perm_size){
  permutations::as.cycle(sample(perm_size, 2, replace=FALSE))
}