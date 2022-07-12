#' Find the Gaussian model Invariant by Permutation Symmetry
#'
#' Uses one of optimization algorithms to find the permutation that maximizes the likelihood of observed data.
#'
#' @param g object of `gips` class
#' @param max_iter number of iterations for an algorithm to perform. At least 2. For `optimizer=="MH"` has to be finite; for `optimizer=="BG"`, can be infinite.
#' @param return_probabilities boolean. TRUE can only be provided for `optimizer=="MH"`. Whether to use `Metropolis_Hastings()` results to calculate posterior probabilities.
#' @param show_progress_bar boolean. Indicate weather or not show the progress bar.
#' @param optimizer the optimizer for the search of the maximum likelihood. Currently the "MH" - Metropolis-Hastings algorithm, or "BG" - best growth algorithm. See #TODO(reference appropriate documentation pages: Metropolis_Hastings and best_growth)
#'
#' @return object of class gips
#' 
#' @export
#'
#' @examples
#' require(MASS)  # for mvrnorm()
#' 
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
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- (t(Z) %*% Z) / number_of_observations  # the theoretical mean is 0
#' 
#' g <- gips(S, number_of_observations)
#' 
#' g <- find_gips(g, max_iter=10, show_progress_bar=FALSE, optimizer="MH")
#' g
#' 
#' if (require(graphics)) {
#'   plot(g, logarithmic_x=TRUE)
#' }
find_gips <- function(g, max_iter, return_probabilities=FALSE,
                      show_progress_bar=TRUE, optimizer="MH"){
  # check the correctness of the g argument
  validate_gips(g)
  
  # extract parameters
  S <- attr(g, "S")
  number_of_observations <- attr(g, "number_of_observations")
  start_perm <- g[[1]]
  delta <- attr(g, "delta")
  D_matrix <- attr(g, "D_matrix")
  
  # check the correctness of the rest of arguments
  if(!(optimizer %in% c("MH", "Metropolis_Hastings", "BG", "best_growth"))){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "`optimizer` must be one of: c('MH', 'Metropolis_Hastings', 'BG', 'best_growth').",
                   "x" = paste0("You provided `optimizer` == ", optimizer, ".")))
  }
  if(!(optimizer %in% c("MH", "Metropolis_Hastings")) && return_probabilities){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "Probabilities can only be provided with the `optimizer = 'Metropolis_Hastings'",
                   "x" = "You should use either `optimizer == Metropolis_Hastings` or `return_probabilities == FLASE`"))
  }
  
  
  if(optimizer %in% c("MH", "Metropolis_Hastings")){
    return(Metropolis_Hastings(S=S, number_of_observations=number_of_observations,
                               max_iter=max_iter, start_perm=start_perm,
                               delta=delta, D_matrix=D_matrix,
                               return_probabilities=return_probabilities,
                               show_progress_bar=show_progress_bar))
  }
  
  if(optimizer %in% c("BG", "best_growth")){
    return(best_growth(S=S, number_of_observations=number_of_observations,
                       max_iter=max_iter, start_perm=start_perm,
                       delta=delta, D_matrix=D_matrix,
                       show_progress_bar=show_progress_bar))
  }
}


check_correctness_of_arguments <- function(S, number_of_observations, max_iter,
                                           start_perm, delta, D_matrix,
                                           return_probabilities, show_progress_bar){
  abord_text <- character(0)
  if(!is.matrix(S))
    abord_text <- c(abord_text,
                    "i" = "`S` must be a matrix.",
                    "x" = paste0("You provided `S` with type == (",
                                 paste(typeof(S), collapse = ", "),
                                 ")."))
  else if(ncol(S) != nrow(S))
    abord_text <- c(abord_text,
                    "i" = "`S` matrix must be a square matrix.",
                    "x" = paste0("You provided `S` as a matrix, but with different sizes: ",
                                 dim(S)[1], " and ", dim(S)[2], "."))
  if(is.null(number_of_observations))
    abord_text <- c(abord_text,
                    "i" = "`number_of_observations` must not be `NULL`.",
                    "x" = "Your provided `number_of_observations` is NULL.")
  else if(number_of_observations < 1)
    abord_text <- c(abord_text,
                    "i" = "`number_of_observations` must be at least 1.",
                    "x" = paste0("You provided `number_of_observations` == ",
                                 number_of_observations, "."))
  else if(!is.wholenumber(number_of_observations))
    abord_text <- c(abord_text,
                    "i" ="`number_of_observations` must be a whole number.",
                    "x" = paste0("You provided `number_of_observations` == ",
                                 number_of_observations, "."))
  if(!(is.infinite(max_iter) || is.wholenumber(max_iter)))
    abord_text <- c(abord_text,
                    "i" ="`max_iter` must be either infinite (for best_growth optimizer) or a whole number.",
                    "x" = paste0("You provided `max_iter` == ", max_iter, "."))
  else if(max_iter < 2)  # TODO(Make it work for max_iter == 1)
    abord_text <- c(abord_text,
                    "i" = "`max_iter` must be at least 2.",
                    "x" = paste0("You provided `max_iter` == ", max_iter, "."))
  if(!(permutations::is.cycle(start_perm) || inherits(start_perm, 'gips_perm')))
    abord_text <- c(abord_text,
                    "i" = "`start_perm` must be the output of `gips_perm()` function, or of class `cycle` form `permutations` package.",  # this is not true, but it is close enough
                    "x" = paste0("You provided `start_perm` with class == (",
                                 paste(class(start_perm), collapse = ", "),
                                 ")."))
  else if(!(permutations::is.cycle(start_perm) || attr(start_perm, 'size') == ncol(S)))
    abord_text <- c(abord_text,
                    "i" = "`start_perm` must have the `size` attribute equal to the shape of a square matrix `S`",
                    "x" = paste0("You provided `start_perm` with `size` == ",
                                 attr(start_perm, 'size'),
                                 ", but the `S` matrix you provided has ",
                                 ncol(S), " columns."))
  if(is.null(delta))
    abord_text <- c(abord_text,
                    "i" ="`delta` must not be `NULL`.",
                    "x" = "Your provided `delta` is a NULL.")
  else if(delta <= 2)
    abord_text <- c(abord_text,
                    "i" ="`delta` must be strictly bigger than 2.",
                    "x" = paste0("You provided `delta` == ", delta, "."))
  if(!(is.null(D_matrix) || is.matrix(D_matrix)))
    abord_text <- c(abord_text,
                    "i" ="`D_matrix` must either be `NULL` or a matrix.",
                    "x" = paste0("You provided `D_matrix` with type == (",
                                 paste(typeof(D_matrix), collapse = ", "),
                                 ")."))
  else if(!(is.null(D_matrix) || ncol(D_matrix) == nrow(D_matrix)))
    abord_text <- c(abord_text,
                    "i" ="`D_matrix` must either be `NULL` or a square matrix.",
                    "x" = paste0("You provided `D_matrix` as a matrix, but with different sizes: ",
                                 ncol(D_matrix), " and ", nrow(D_matrix), "."))
  else if(!(is.null(D_matrix) || ncol(S) == ncol(D_matrix)))
    abord_text <- c(abord_text,
                    "i" ="`S` must be a square matrix with the same shape as a square matrix `D_matrix`.",
                    "x" = paste0("You provided `S` with shape ",
                                 ncol(S), " and ", nrow(S),
                                 ", but also `D_matrix` with shape ",
                                 ncol(D_matrix), " and ", nrow(D_matrix), "."))
  if(!is.logical(return_probabilities))
    abord_text <- c(abord_text,
                    "i" ="`return_probabilities` must be a logic value (TRUE or FALSE).",
                    "x" = paste0("You provided `return_probabilities` with type == (",
                                 paste(typeof(return_probabilities), collapse = ", "),
                                 ")."))
  if(!is.logical(show_progress_bar))
    abord_text <- c(abord_text,
                    "i" ="`show_progress_bar` must be a logic value (TRUE or FALSE).",
                    "x" = paste0("You provided `show_progress_bar` with type == (",
                                 paste(typeof(show_progress_bar), collapse = ", "),
                                 ")."))
  
  if(length(abord_text) > 0){
    abord_text <- c(paste0("There were ", length(abord_text)/2,
                           " problems identified with provided arguments:"),
                    abord_text)
    
    if(length(abord_text) > 11){
      abord_text <- c(abord_text[1:11],
                      paste0("... and ", (length(abord_text)-1)/2 - 5, " more problems"))
    }
    
    rlang::abort(abord_text)
  }
}



Metropolis_Hastings <- function(S, number_of_observations, max_iter, start_perm=NULL,
                                delta=3, D_matrix=NULL, return_probabilities=FALSE,
                                show_progress_bar=TRUE){
  if(is.null(start_perm)){
    start_perm <- permutations::id
  }
  
  check_correctness_of_arguments(S=S, number_of_observations=number_of_observations,
                                 max_iter=max_iter, start_perm=start_perm,
                                 delta=delta, D_matrix=D_matrix,
                                 return_probabilities=return_probabilities,
                                 show_progress_bar=show_progress_bar)
  
  if(!inherits(start_perm, 'gips_perm')){
    start_perm <- gips_perm(start_perm, nrow(S))  # now we know the `S` is a matrix
  }
  
  if(is.infinite(max_iter)){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "`max_iter` in `Metropolis_Hastings()` must be finite",
                   "x" = paste0("You provided `max_iter` == ", max_iter, ".")))
  }
  
  perm_size <- dim(S)[1]
  if(permutations::is.cycle(start_perm))
      start_perm <- gips_perm(start_perm, perm_size)
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }
  
  my_goal_function <- function(perm){
    log_likelihood_of_perm(perm, S=S,  # We recommend to use the `log_likelihood_of_gips()` function
                           number_of_observations=number_of_observations,
                           delta=delta, D_matrix=D_matrix)
  }

  acceptance <- rep(FALSE, max_iter)
  log_likelihood_values <- rep(0, max_iter)
  visited_perms <- list()
  visited_perms[[1]] <- start_perm

  if(show_progress_bar)
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  log_likelihood_values[1] <- my_goal_function(visited_perms[[1]])

  found_perm <- start_perm
  found_perm_log_likelihood <- log_likelihood_values[1]

  Uniformly_drawn_numbers <- stats::runif(max_iter, min = 0, max = 1)

  # main loop
  for (i in 1:(max_iter-1)){
    if(show_progress_bar)
      utils::setTxtProgressBar(progressBar, i)
    
    e <- runif_transposition(perm_size)
    perm_proposal <- compose_with_transposition(visited_perms[[i]], e)

    goal_function_perm_proposal <- my_goal_function(perm_proposal)
    if(is.nan(goal_function_perm_proposal) | is.infinite(goal_function_perm_proposal)){
      # See ISSUE#5; We hope the introduction of log calculations have stopped this problem.
      warning("gips is yet unable to process this S matrix. We think it can only happen for dim(S)[1] > 500. If it is not the case for you, please get in touch with us on ISSUE#5")

      break()
    }

    # if goal_function_perm_proposal > log_likelihood_values[i], then it is true, because Uniformly_drawn_numbers[i] \in [0,1]
    if(Uniformly_drawn_numbers[i] < exp(goal_function_perm_proposal-log_likelihood_values[i])){ # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2). That means this is Metropolis algorithm, not necessary Metropolis-Hastings.
      visited_perms[[i+1]] <- perm_proposal
      log_likelihood_values[i+1] <- goal_function_perm_proposal
      acceptance[i] <- TRUE

      if(found_perm_log_likelihood < log_likelihood_values[i+1]){
        found_perm_log_likelihood <- log_likelihood_values[i+1]
        found_perm <- visited_perms[[i+1]]
      }
    }
    else{
      visited_perms[[i+1]] = visited_perms[[i]]
      log_likelihood_values[i+1] <- log_likelihood_values[i]
    }
  }

  if(show_progress_bar)
    close(progressBar)

  function_calls <- length(log_likelihood_values)
  
  if(return_probabilities){
    probabilities <- estimate_probabilities(visited_perms)
  }else{
    probabilities <- NULL
  }
  optimization_info <- list("acceptance_rate" = mean(acceptance),
                            "log_likelihood_values" = log_likelihood_values,
                            "visited_perms" = visited_perms,
                            "last_perm" = visited_perms[[function_calls]],
                            "last_perm_log_likelihood" = log_likelihood_values[function_calls],
                            "iterations_performed" = i,
                            "optimization_algorithm_used" = "Metropolis_Hastings",
                            "post_probabilities" = probabilities,
                            "did_converge" = NULL,
                            "best_perm_log_likelihood" = found_perm_log_likelihood)
  
  
  new_gips(list(found_perm), S, number_of_observations, delta,
           D_matrix, optimization_info)
}


best_growth <- function(S, number_of_observations, max_iter=5,
                        start_perm=NULL,
                        delta=3, D_matrix=NULL,
                        show_progress_bar=TRUE){
  if(is.null(start_perm)){
    start_perm <- permutations::id
  }
  
  check_correctness_of_arguments(S=S, number_of_observations=number_of_observations,
                                 max_iter=max_iter, start_perm=start_perm,
                                 delta=delta, D_matrix=D_matrix,
                                 return_probabilities=FALSE,
                                 show_progress_bar=show_progress_bar)
  
  if(!inherits(start_perm, 'gips_perm')){
    start_perm <- gips_perm(start_perm, nrow(S))  # now we know the `S` is a matrix
  }
  
  if(show_progress_bar && is.infinite(max_iter)){
    stop("Progress bar is not yet supported for infinite max_iter. Rerun the algorithm with show_progress_bar=FALSE or finite max_iter. For more information see ISSUE#8.") # See ISSUE#8
  }
  
  if(show_progress_bar)
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)

  perm_size <- dim(S)[1]

  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }
  

  my_goal_function <- function(perm){
    log_likelihood_of_perm(perm, S=S,  # We recommend to use the `log_likelihood_of_gips()` function
                           number_of_observations=number_of_observations,
                           delta=delta, D_matrix=D_matrix)
  }

  goal_function_best_logvalues <- numeric(0)
  log_likelihood_values <- numeric(0)
  
  # init
  speciments <- list()
  speciments[[1]] <- start_perm
  goal_function_best_logvalues[1] <- my_goal_function(speciments[[1]])
  log_likelihood_values[1] <- goal_function_best_logvalues[1]

  # mail loop
  iteration <- 0
  did_converge <- FALSE
  while(iteration <= max_iter-1){
    iteration <- iteration + 1
    if(show_progress_bar)
      utils::setTxtProgressBar(progressBar, iteration)

    best_neighbour <- NULL
    best_neighbour_value <- -Inf
    for(i in 1:(perm_size-1)){
      for(j in (i+1):perm_size){
        neighbour <- compose_with_transposition(speciments[[iteration]], c(i, j))
        neighbour_value <- my_goal_function(neighbour)
        log_likelihood_values[length(log_likelihood_values) + 1] <- neighbour_value

        if(neighbour_value > best_neighbour_value){
          best_neighbour_value <- neighbour_value
          best_neighbour <- neighbour
        }
      }
    }

    if(best_neighbour_value > goal_function_best_logvalues[iteration]){
      goal_function_best_logvalues[iteration + 1] <- best_neighbour_value
      speciments[[iteration + 1]] <- best_neighbour
    }else{
      did_converge <- TRUE
      break
    }
  }

  if(show_progress_bar)
    close(progressBar)

  if(!did_converge){
    rlang::warn(paste0("Algorithm did not converge in ", iteration, # now, iteration == max_iter
                       " iterations! Try one more time with starting_perm = output$found_perm")) # TODO(there will be a function `continue(bg)`; see ISSUE#11)
  }else{
    goal_function_best_logvalues <- goal_function_best_logvalues[1:iteration]
    if(show_progress_bar)
      print(paste0("Algorithm did converge in ", iteration, " iterations"))
  }
  
  function_calls <- length(log_likelihood_values)
  
  optimization_info <- list("acceptance_rate" = 1/choose(perm_size, 2),
                            "log_likelihood_values" = log_likelihood_values,
                            "visited_perms" = speciments,
                            "last_perm" = speciments[[iteration]],
                            "last_perm_log_likelihood" = goal_function_best_logvalues[iteration],
                            "iterations_performed" = iteration,
                            "optimization_algorithm_used" = "best_growth",
                            "post_probabilities" = NULL,
                            "did_converge" = did_converge,
                            "best_perm_log_likelihood" = goal_function_best_logvalues[iteration])
  
  new_gips(list(speciments[[iteration]]), S, number_of_observations, delta,
           D_matrix, optimization_info)
}






