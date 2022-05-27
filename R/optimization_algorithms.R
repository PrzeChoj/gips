#' Metropolis-Hastings algorithm
#'
#' Uses Metropolis-Hastings algorithm to find the permutation that maximizes the likelihood of observed data.
#'
#' @param U matrix, sum of outer products of data. `U` = sum(t(Z) %*% Z), where Z is observed data.
#' @param n_number number of data points that `U` is based on.
#' @param max_iter number of iterations for an algorithm to perform. At least 2.
#' @param start_perm starting permutation for the algorithm; an element of class "cycle". When NULL, identity permutation is taken.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of the same size as `U`. When NULL, the identity matrix is taken.
#' @param return_probabilities boolean. Whether to use MH results to calculate posterior probabilities.
#' @param show_progress_bar boolean, indicate weather or not show the progress bar.
#'
#' @return object of class gips; list of 8 items: `acceptance_rate`,
#' `goal_function_logvalues`, `points`, `found_point`,
#' `found_point_function_logvalue`, `last_point`,
#' `last_point_function_logvalue`, `iterations_performed`.
#'
#' @export
#'
#' @examples
#' perm_size <- 6
#' mu <- numeric(perm_size)
#' # sigma is a matrix invariant under permutation (1,2,3,4,5,6)
#' sigma_matrix <- matrix(data = c(1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
#'                                 0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
#'                                 0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
#'                                 0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
#'                                 0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
#'                                 0.8, 0.6, 0.4, 0.6, 0.8, 1.0),
#'                        nrow=perm_size, byrow = TRUE)
#' n_number <- 13
#' Z <- MASS::mvrnorm(n_number, mu = mu, Sigma = sigma_matrix)
#' U <- (t(Z) %*% Z)
#' start_perm <- permutations::id
#' mh <- MH(U=U, n_number=n_number, max_iter=10, start_perm=start_perm,
#'          show_progress_bar=FALSE)
#' plot(mh)
MH <- function(U, n_number, max_iter, start_perm=NULL,
               delta=3, D_matrix=NULL, return_probabilities=FALSE,
               show_progress_bar=TRUE){
  if(is.null(start_perm)){
    start_perm <- permutations::id
  }
  stopifnot(permutations::is.cycle(start_perm),
            is.matrix(U),
            dim(U)[1] == dim(U)[2],
            max_iter >= 2) # TODO(Make it work for max_iter == 1)
  perm_size <- dim(U)[1]
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }
  stopifnot(is.matrix(D_matrix),
            dim(D_matrix)[1] == dim(D_matrix)[2])

  acceptance <- rep(FALSE, max_iter)
  goal_function_logvalues <- rep(0, max_iter)
  points <- list()
  points[[1]] <- start_perm

  if(show_progress_bar)
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  goal_function_logvalues[1] <- goal_function(points[[1]],
                                              n_number, U,
                                              delta=delta, D_matrix=D_matrix)

  found_point <- start_perm
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
      #browser()  # See ISSUE#5; We hope the introduction of log calculations
      # will stop this problem.
      warning("gips is yet unable to process this U matrix. We think it can only happen for dim(U)[1] > 500. If it is not the case for you, please get in touch with us on ISSUE#5")

      break()
    }

    # if goal_function_perm_proposal > goal_function_logvalues[i], then it is true
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

  function_calls <- length(goal_function_logvalues)

  out <- list("acceptance_rate"=mean(acceptance),
              "goal_function_logvalues"=goal_function_logvalues,
              "points"=points,
              "found_point"=found_point,
              "found_point_function_logvalue"=found_point_function_logvalue,
              "last_point"=points[[function_calls]],
              "last_point_function_logvalue"=found_point_function_logvalue[function_calls],
              "iterations_performed"=i)

  class(out) <- c("gips", "list")

  if(return_probabilities){
      probabilities <- estimate_probabilities(points)
      out$post_probabilities <- estimate_probabilities(out)
    }
  out
}


#' Best growth algorithm
#'
#' Uses best growth algorithm to find the permutation that maximizes the likelihood of observed data.
#'
#' @param U matrix, sum of outer products of data. `U` = sum(t(Z) %*% Z), where Z is observed data.
#' @param n_number number of data points that `U` is based on.
#' @param max_iter number of iterations for an algorithm to perform. At least 2.
#' @param start_perm starting permutation for the algorithm; an element of class "cycle". When NULL, identity permutation is taken.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of the same size as `U`. When NULL, the identity matrix is taken.
#' @param show_progress_bar boolean, indicate weather or not show the progress bar.
#'
#' @return object of class gips; list of 8 items: `acceptance_rate`,
#' `goal_function_logvalues`, `points`, `found_point`,
#' `found_point_function_logvalue`, `last_point`,
#' `last_point_function_logvalue`, `iterations_performed`.
#'
#' @export
#'
#' @examples
#' perm_size <- 6
#' mu <- numeric(perm_size)
#' # sigma is a matrix invariant under permutation (1,2,3,4,5,6)
#' sigma_matrix <- matrix(data = c(1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
#'                                 0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
#'                                 0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
#'                                 0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
#'                                 0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
#'                                 0.8, 0.6, 0.4, 0.6, 0.8, 1.0),
#'                        nrow=perm_size, byrow = TRUE)
#' n_number <- 13
#' Z <- MASS::mvrnorm(n_number, mu = mu, Sigma = sigma_matrix)
#' U <- (t(Z) %*% Z)
#' start_perm <- permutations::id
#' bg <- best_growth(U=U, n_number=n_number, max_iter=5, start_perm=start_perm,
#'                   show_progress_bar=FALSE) # Algorithm did converge in 4 iterations
#' plot(bg)
best_growth <- function(U, n_number, max_iter=5,
                        start_perm=NULL,
                        delta=3, D_matrix=NULL,
                        show_progress_bar=TRUE){
  if(show_progress_bar)
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)

  stopifnot(dim(U)[1] == dim(U)[2], max_iter > 1)
  perm_size <- dim(U)[1]

  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }
  if(is.null(start_perm)){
    start_perm <- permutations::id
  }

  my_goal_function <- function(perm){
    goal_function(perm, n_number, U,
                  delta=delta, D_matrix=D_matrix)
  }

  f_values <- numeric(max_iter)

  # init
  speciments <- list()
  speciments[[1]] <- start_perm
  f_values[1] <- my_goal_function(speciments[[1]])

  # mail loop
  for(iteration in 2:max_iter){
    if(show_progress_bar)
      utils::setTxtProgressBar(progressBar, iteration)

    best_neighbour <- NULL
    best_neighbour_value <- -Inf
    for(i in 1:(perm_size-1)){
      for(j in (i+1):perm_size){
        neighbour <- permutations::as.cycle(speciments[[iteration - 1]] * permutations::as.cycle(c(i, j)))
        neighbour_value <- my_goal_function(neighbour)

        if(neighbour_value > best_neighbour_value){
          best_neighbour_value <- neighbour_value
          best_neighbour <- neighbour
        }
      }
    }

    if(best_neighbour_value > f_values[iteration - 1]){
      f_values[iteration] <- best_neighbour_value
      speciments[[iteration]] <- best_neighbour
    }else{
      iteration <- iteration-1
      break
    }
  }

  if(show_progress_bar)
    close(progressBar)

  if(iteration == max_iter)
    warning("Algorithm did not converge! Try with bigger max_iter and starting_perm == output$found_perm") # TODO(there will be a function `continue(bg)`; see ISSUE#11)
  else{
    f_values <- f_values[1:iteration]
    print(paste0("Algorithm did converge in ", iteration, " iterations"))
  }


  out <- list("acceptance_rate" = 1/choose(perm_size, 2),
              "goal_function_logvalues" = f_values,
              "points" = speciments,
              "found_point" = speciments[[iteration]],
              "found_point_function_logvalue" = f_values[iteration],
              "last_point" = speciments[[iteration]],
              "last_point_function_logvalue" = f_values[iteration],
              "iterations_performed" = iteration)

  class(out) <- c("gips", "list")

  out
}






