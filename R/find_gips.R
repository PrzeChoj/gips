#' Find the Gaussian model Invariant by Permutation Symmetry
#'
#' Uses one of optimization algorithms to find the permutation that maximizes the likelihood of observed data.
#'
#' @param g object of `gips` class
#' @param max_iter number of iterations for an algorithm to perform. At least 2. For `optimizer=="MH"` has to be finite; for `optimizer=="BG"`, can be infinite; for `optimizer=="BF"` it is not used.
#' @param return_probabilities boolean. TRUE can only be provided for `optimizer=="MH"`. Whether to use `Metropolis_Hastings()` results to calculate posterior probabilities.
#' @param show_progress_bar boolean. Indicate weather or not show the progress bar.
#' @param optimizer the optimizer for the search of the maximum likelihood. Currently the "MH" - Metropolis-Hastings algorithm, or "BG" - best growth algorithm, or "BF" - brute force algorithm, or "continue" to continue the optimization performed on the `g` object (see Examples). By default, NA that is changed into "MH" when `g` is unoptimized and "continue", when `g` is optimized. See #TODO(In "Details" explain: Metropolis_Hastings and best_growth and brute_force)
#'
#' @return object of class gips
#' 
#' @export
#'
#' @examples
#' require("MASS")  # for mvrnorm()
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
#' g_opt <- find_gips(g, max_iter=10, show_progress_bar=FALSE, optimizer="MH")
#' g_opt
#' 
#' g_opt2 <- find_gips(g_opt, max_iter=10, show_progress_bar=FALSE, optimizer="continue")
#' 
#' if (require("graphics")) {
#'   plot(g_opt2, type="both", logarithmic_x=TRUE)
#' }
#' 
#' g_opt_BF <- find_gips(g, show_progress_bar=FALSE, optimizer="BF")
#' summary(g_opt_BF)
find_gips <- function(g, max_iter=NA, return_probabilities=FALSE,
                      show_progress_bar=TRUE, optimizer=NA){
  # check the correctness of the g argument
  validate_gips(g)
  
  # check the correctness of the rest of arguments
  if(length(optimizer) > 1)
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "`optimizer` must be the character vector of length 1. Must be one of: c('MH', 'Metropolis_Hastings', 'BG', 'best_growth', 'continue').",
                   "x" = paste0("You provided `optimizer == (",
                                paste0(optimizer, collapse = ", "), ")`."),
                   "i" = "Did You misspelled the optimizer name?"))
  # default optimizer:
  if(is.na(optimizer)){
    optimizer <- ifelse(is.null(attr(g, "optimization_info")),
                        "MH", "continue")
    
    rlang::inform(c("You used the default value of the 'optimizer' argument in `find_gips()`.",
                    "i" = paste0("The 'optimizer = NA' was automatically changed to 'optimizer = \"",
                                 optimizer, "\"'.")))
  }
  if(!(optimizer %in% c("MH", "Metropolis_Hastings", "BG", "best_growth", "BF", "brute_force", "full", "continue"))){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "`optimizer` must be one of: c('MH', 'Metropolis_Hastings', 'BG', 'best_growth', 'BF', 'brute_force', 'full', 'continue').",
                   "x" = paste0("You provided `optimizer == ", optimizer, "`."),
                   "i" = "Did You misspelled the optimizer name?"))
  }
  if(!(optimizer %in% c("BF", "brute_force", "full")) &&
     is.na(max_iter)){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "`max_iter = NA` can be provided only for `optimizer` one of: c(BF', 'brute_force', 'full'). For any other, `max_iter` must be a whole number bigger than 1.",
                   "x" = paste0("You provided `optimizer == ", optimizer, "` and `max_iter = NA`."),
                   "i" = "Did You forgot to set the `max_iter`?",
                   "i" = "Did You misspelled the optimizer name?"))
  }
  
  continue_optimization <- (optimizer == "continue")
  if(continue_optimization){
    if(is.null(attr(g, "optimization_info"))){
      rlang::abort(c("There was a problem identified with provided arguments:",
                     "i" = "`optimizer == 'continue'` can be provided only with optimized gips object `g`.",
                     "x" = "You provided `optimizer == 'continue'`, but the gips object `g` is not optimized.",
                     "i" = "Did You provided wrong `gips` object?",
                     "i" = "Did You want to call another optimizer like 'MH' or 'BG'?"))
    }
    optimizer <- attr(g, "optimization_info")[["optimization_algorithm_used"]][length(attr(g, "optimization_info")[["optimization_algorithm_used"]])]  # this is the last used optimizer
    
    if(optimizer == "brute_force")
      rlang::abort(c("There was a problem identified with provided arguments:",
                     "i" = "`optimizer == 'continue'` cannot be provided after optimizating with `optimizer == 'brute_force'`, because the whole space was already browsed.",
                     "x" = "You provided `optimizer == 'continue'`, but the gips object `g` was optimized with brute_force optimizer. Better permutation will not be found."))
  }
  
  if(!(optimizer %in% c("MH", "Metropolis_Hastings")) && return_probabilities){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "Probabilities can only be returned with the `optimizer == 'Metropolis_Hastings'`",
                   "x" = "You provided both `optimizer != Metropolis_Hastings` and `return_probabilities == TRUE`!",
                   "i" = "Did You want to use `optimizer == Metropolis_Hastings` or `return_probabilities == FLASE`?"))
  }
  
  # inform, that user can use "BF"
  if((optimizer %in% c("MH", "Metropolis_Hastings")) && (max_iter*10 >= prod(1:ncol(attr(g, "S"))))){
    rlang::inform(c(paste0("You called optimization with Metropolis_Hastings algorith with ",
                           max_iter, " iterations."),
                    "i" = paste0("Consider using `optimizer = 'brute_force'`, because it will use ",
                                 ncol(attr(g, "S")), "! (factorial) = ", prod(1:ncol(attr(g, "S"))),
                                 " iterations and will browse all permutations, therefore it will definitely find the maximum likelihood estimator.")))
  }
  
  # extract parameters
  S <- attr(g, "S")
  number_of_observations <- attr(g, "number_of_observations")
  if(continue_optimization){  # the `ifelse()` function cannot be used because the objects are lists
    start_perm <- attr(g, "optimization_info")[["last_perm"]]
  }else{
    start_perm <- g[[1]]
  }
  delta <- attr(g, "delta")
  D_matrix <- attr(g, "D_matrix")
  
  start_time <- Sys.time()
  
  if(optimizer %in% c("MH", "Metropolis_Hastings")){
    gips_optimized <- Metropolis_Hastings(S=S, number_of_observations=number_of_observations,
                                          max_iter=max_iter, start_perm=start_perm,
                                          delta=delta, D_matrix=D_matrix,
                                          return_probabilities=return_probabilities,
                                          show_progress_bar=show_progress_bar)
  }else if(optimizer %in% c("BG", "best_growth")){
    gips_optimized <- best_growth(S=S, number_of_observations=number_of_observations,
                                  max_iter=max_iter, start_perm=start_perm,
                                  delta=delta, D_matrix=D_matrix,
                                  show_progress_bar=show_progress_bar)
  }else if(optimizer %in% c("BF", "brute_force", "full")){
    gips_optimized <- brute_force(S=S, number_of_observations=number_of_observations,
                                  delta=delta, D_matrix=D_matrix,
                                  show_progress_bar=show_progress_bar)
  }
  
  end_time <- Sys.time()
  attr(gips_optimized, "optimization_info")[["optimization_time"]] <- end_time - start_time
  attr(gips_optimized, "optimization_info")[["full_optimization_time"]] <- end_time - start_time
  
  return(combine_gips(g, gips_optimized))
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
                   "x" = paste0("You provided `max_iter == ", max_iter, "`.")))
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
    if(is.nan(goal_function_perm_proposal) || is.infinite(goal_function_perm_proposal)){
      # See ISSUE#5; We hope the introduction of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
                    "x"=paste0("The likelihood value of ", ifelse(is.nan(goal_function_perm_proposal), "NaN", "Inf"), " occured!"),
                    "i"="We think it can only happen for ncol(S) > 500. If it is not the case for You, please get in touch with us on ISSUE#5."))
      
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
      log_likelihood_values[i+1] <- log_likelihood_values[i]  # TODO(Do we really want to forget the calculated values? The algorithm BG works differently)
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
                            "best_perm_log_likelihood" = found_perm_log_likelihood,
                            "optimization_time" = NA,
                            "full_optimization_time" = NA)
  
  
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
    rlang::warn(c(paste0("Best Growth algorithm did not converge in ", iteration, " iterations!"), # now, iteration == max_iter
                  "i" = "We reccomend to run the `find_gips()` one more time on the conquered output")) # TODO(There will be a function `continue(bg)`; see ISSUE#11)
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
                            "best_perm_log_likelihood" = goal_function_best_logvalues[iteration],
                            "optimization_time" = NA,
                            "full_optimization_time" = NA)
  
  new_gips(list(speciments[[iteration]]), S, number_of_observations, delta,
           D_matrix, optimization_info)
}


brute_force <- function(S, number_of_observations,
                        delta=3, D_matrix=NULL,
                        show_progress_bar=TRUE){
  check_correctness_of_arguments(S=S, number_of_observations=number_of_observations,
                                 max_iter=5, start_perm=permutations::id,  # max_iter and start_perm are not important for optimization with brute_force
                                 delta=delta, D_matrix=D_matrix,
                                 return_probabilities=FALSE,
                                 show_progress_bar=show_progress_bar)
  
  perm_size <- dim(S)[1]
  
  if(perm_size>30){  # TODO(Test and set smaller, reasonable value. 15?)
    rlang::abort(c("Optimizer 'brute_force' cannot browse such a big permutional space.",
                   "x" = paste0("You provided a space with size ", perm_size,
                                "! (factorial), which has ", prod(1:perm_size),
                                " elements."),
                   "i" = "Do You want to use other optimizer for such a big slace? For example 'Metropolis_Hastings' or 'best_growth'?"))
  }
  
  if(show_progress_bar)
    progressBar <- utils::txtProgressBar(min = 0, max = prod(1:perm_size), initial = 1)
  
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = perm_size)
  }
  
  my_goal_function <- function(perm){
    log_likelihood_of_perm(perm, S=S,  # We recommend to use the `log_likelihood_of_gips()` function
                           number_of_observations=number_of_observations,
                           delta=delta, D_matrix=D_matrix)
  }
  
  # main loop
  all_perms_list <- permutations::allperms(perm_size)
  all_perms_list <- permutations::as.cycle(all_perms_list)
  log_likelihood_values <- sapply(1:length(all_perms_list), function(i){
    if(show_progress_bar){
      utils::setTxtProgressBar(progressBar, i)
    }
    this_perm <- permutations::cycle(list(all_perms_list[[i]]))
    my_goal_function(this_perm)
  })
  
  if(show_progress_bar)
    close(progressBar)
  
  best_perm <- gips_perm(permutations::cycle(list(all_perms_list[[which.max(log_likelihood_values)]])), perm_size)
  
  optimization_info <- list("acceptance_rate" = NULL,
                            "log_likelihood_values" = log_likelihood_values,
                            "visited_perms" = all_perms_list,
                            "last_perm" = NULL,
                            "last_perm_log_likelihood" = NULL,
                            "iterations_performed" = prod(1:perm_size),
                            "optimization_algorithm_used" = "brute_force",
                            "post_probabilities" = NULL,
                            "did_converge" = TRUE,
                            "best_perm_log_likelihood" = log_likelihood_values[which.max(log_likelihood_values)],
                            "optimization_time" = NA,
                            "full_optimization_time" = NA)
  
  new_gips(list(best_perm), S, number_of_observations, delta,
           D_matrix, optimization_info)
}



#' Combining 2 gips objects
#' 
#' g2 was optimized with a single optimization method. g1 was potentially non-optimized or optimized once, or optimized multiple times.
#' If g2 was optimized with "brute_force", forget the g1.
#' 
#' @noRd
combine_gips <- function(g1, g2){
  if(is.null(attr(g1, "optimization_info")) ||
     attr(g2, "optimization_info")[["optimization_algorithm_used"]] == "brute_force"){  # when brute_force was used, forget the initial optimization
    return(g2)
  }
  
  # g1 is also an effect of optimization.
  optimization_info1 <- attr(g1, "optimization_info")
  optimization_info2 <- attr(g2, "optimization_info")
  
  n1 <- length(optimization_info1[["log_likelihood_values"]])
  n2 <- length(optimization_info2[["log_likelihood_values"]])
  
  visited_perms <- c(optimization_info1[["visited_perms"]], optimization_info2[["visited_perms"]])  # WoW, one can use `c()` to combine lists!
  optimization_algorithm_used <- c(optimization_info1[["optimization_algorithm_used"]], optimization_info2[["optimization_algorithm_used"]])
  
  if(all(optimization_algorithm_used == "Metropolis_Hastings") &&
     !is.null(optimization_info1[["post_probabilities"]]) &&
     !is.null(optimization_info2[["post_probabilities"]])){
    post_probabilities <- estimate_probabilities(visited_perms)  # TODO(This can be combined more optimally, but I (Adam) think this will be rarely done nevertheless.)
  }else{
    post_probabilities <- NULL
  }
  
  optimization_info_new <- list("acceptance_rate" = (n1*optimization_info1[["acceptance_rate"]] + n2*optimization_info2[["acceptance_rate"]])/(n1+n2),
                                "log_likelihood_values" = c(optimization_info1[["log_likelihood_values"]], optimization_info2[["log_likelihood_values"]]),
                                "visited_perms" = visited_perms,
                                "last_perm" = optimization_info2[["last_perm"]],
                                "last_perm_log_likelihood" = optimization_info2[["last_perm_log_likelihood"]],
                                "iterations_performed" = c(optimization_info1[["iterations_performed"]], optimization_info2[["iterations_performed"]]),
                                "optimization_algorithm_used" = optimization_algorithm_used,
                                "post_probabilities" = post_probabilities,
                                "did_converge" = optimization_info2[["did_converge"]],
                                "best_perm_log_likelihood" = max(optimization_info1[["best_perm_log_likelihood"]], optimization_info2[["best_perm_log_likelihood"]]),
                                "optimization_time" = c(optimization_info1[["optimization_time"]], optimization_info2[["optimization_time"]]),
                                "full_optimization_time" = optimization_info1[["full_optimization_time"]] + optimization_info2[["full_optimization_time"]])
  
  attr(g2, "optimization_info") <- optimization_info_new
  
  g2
}






