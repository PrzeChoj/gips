#' constructor of the `gips` class.
#' 
#' Create the `gips` object.
#' This object will consists data and all other information needed to find the invariant group.
#' The optimization itself will not be performed. To do it, one have to call the TODO() function. See examples below.
#' 
#' @param S matrix, estimated covariance matrix. When Z is observed data: `S = (t(Z) %*% Z) / number_of_observations`, if one know the theoretical mean is 0; # TODO(What if one have to estimate the theoretical mean with the empirical mean)
#' @param number_of_observations number of data points that `S` is based on.
#' @param delta hyper-parameter of a Bayesian model. Has to be bigger than 2.
#' @param D_matrix hyper-parameter of a Bayesian model. Square matrix of the same size as `S`. When NULL, the identity matrix is taken.
#' 
#' @return Object of class gips.
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
#' S <- (t(Z) %*% Z)/number_of_observations  # the theoretical mean is 0
#' 
#' # TODO()
gips <- function(S, number_of_observations, delta=3, D_matrix=NULL){
  check_correctness_of_arguments(S=S, number_of_observations=number_of_observations,
                                 max_iter=2, start_perm=permutations::id,
                                 delta=delta, D_matrix=D_matrix,
                                 return_probabilities=FALSE, show_progress_bar=FALSE)
  
  gips_perm_object <- gips_perm(permutations::id, nrow(S))
  
  if(is.null(D_matrix)){
    D_matrix <- diag(nrow = ncol(S))
  }
  
  validate_gips(new_gips(gips_perm_object, S, number_of_observations,
                         delta=delta, D_matrix=D_matrix,
                         optimization_info=NULL))
}


# TODO(Documentation)
#' @describeIn gips Constructor
#'
#' Only intended for low-level use.
#'
#' @param gips_perm_object of class `gips_perm`. The base object for the `gips` class.
#' @param optimization_info NULL or the list with information about the optimizaion process.
#'
#' @export
new_gips <- function(gips_perm_object, S, number_of_observations, delta,
                     D_matrix, optimization_info){
  stopifnot(inherits(gips_perm_object, "gips_perm"),
            is.matrix(S),
            is.wholenumber(number_of_observations),
            is.numeric(delta),
            is.matrix(D_matrix),
            is.null(optimization_info) || is.list(optimization_info))
  
  structure(gips_perm_object, S=S, number_of_observations=number_of_observations,
            delta=delta, D_matrix=D_matrix, optimization_info=optimization_info,
            class=c("gips", "gips_perm"))
}


# TODO(Documentation)
validate_gips <- function(x){
  perm <- unclass(x)
  S <- attr(x, "S")
  number_of_observations <- attr(x, "number_of_observations")
  delta <- attr(x, "delta")
  D_matrix <- attr(x, "D_matrix")
  optimization_info <- attr(x, "optimization_info")
  
  
  if(!inherits(x, "gips_perm")){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "The `x` must be also of a `gips_perm` class.",
                   "x" = paste0("You provided `x` with class == (",
                                paste(class(x), collapse = ", "),
                                ").")))
  }
  
  class(perm) <- "gips_perm"
  attr(perm, "size") <- attr(x, "size")
  check_correctness_of_arguments(S=S, number_of_observations=number_of_observations,
                                 max_iter=2, start_perm=perm,
                                 delta=delta, D_matrix=D_matrix,
                                 return_probabilities=FALSE, show_progress_bar=FALSE)
  if(!(is.null(optimization_info) || is.list(optimization_info))){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "The `optimization_info` value must be either a NULL, or a list.",
                   "x" = paste0("You provided `x` with type == (",
                                paste(typeof(perm), collapse = ", "),
                                ").")))
  }
  
  # TODO(Validate the `optimization_info` more carefully (when it will be ready))
  
  x
}


# TODO(The base object printed as the gips_perm. The "size" attr omitted)
str.gips <- function(object, ...){
  str(object)
}


#' Printing gips object
#' 
#' Printing function for gips class.
#' @param x object of class gips.
#' @param log_value logical. Weather to print the exp of a value of a \code{\link{log_likelihood_of_perm}} or leave it in logarithmic form.
#' @param ... additional arguments passed to \code{\link{cat}}.
#' 
#' @return Invisible NULL.
#' @export
print.gips <- function(x, log_value = TRUE, ...){
  # TODO(it is not likelihood, but sth proportional to it. See #ISSUE11)
  
  value_part <- ifelse(log_value,
                       paste0(" with log likelihood ",
                              x[["found_perm_log_likelihood"]]),
                       paste0(" with likelihood ",
                              exp(x[["found_perm_log_likelihood"]])))
  cat(paste0("Optimization algorithm ",
             x[["optimization_algorithm_used"]], " after ",
             length(x[["log_likelihood_values"]]),
             " iterations found permutation ",
             x[["found_perm"]],
             value_part),
      ...)
}




#' Plot optimization convergence
#' 
#' Plot method for gips objects.
#' @param x Object of class gips.
#' @param type Character. A type of a plot. Either "all", "best" or "both". For "all" plots likelihood for all visited state. For "best" the biggest likelihood up to the point are plotted. For "both" both lines are plotted.
#' @param logarithmic_y boolean.
#' @param logarithmic_x boolean.
#' @param color Vector of olors to be used to plot lines.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend boolean.
#' @param ylim Limits of y axis. When \code{NULL}, the minimum and maximum of the \code{\link{log_likelihood_of_perm}} is taken.
#' @param ... additional arguments passed to \code{\link{print}}.
#' 
#' @return Invisible NULL.
#' @export
plot.gips <- function(x, type="both",
                      logarithmic_y=TRUE,
                      logarithmic_x=FALSE,
                      color=NULL,
                      title_text="Convergence plot",
                      xlabel=NULL,
                      ylabel=NULL, show_legend=TRUE,
                      ylim=NULL, ...){
  # TODO(it is not likelihood, but sth proportional to it. See #ISSUE11)
  
  if (!requireNamespace("graphics", quietly = TRUE)) {
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "Package \"graphics\" must be installed to use this function.",
                   "x" = "Package \"graphics\" seems to be unavailable."))
  }
  
  if(!(type %in% c("all", "best", "both"))){
    rlang::abort(c("There was a problem identified with provided arguments:",
                   "i" = "`type` must be one of: c('all', 'best', 'both').",
                   "x" = paste0("You provided `type` == ", type, ".")))
  }
  
  if(is.null(ylabel)){
    ylabel <- ifelse(logarithmic_y,
                     "log likelihood",
                     "likelihood")
  }
  if(is.null(xlabel)){
    xlabel <- ifelse(logarithmic_x,
                     "log10 of number of function calls",
                     "number of function calls")
  }
  if(is.null(color)){
    if(type == "both"){
      color <- c("blue", "red")
    }else{
      color <- "blue"
    }
  }
  if(logarithmic_y){
    y_values_from <- x[["log_likelihood_values"]] # values of likelihood are logarithmic by default
  }else{
    y_values_from <- exp(x[["log_likelihood_values"]])
  }

  y_values_max <- cummax(y_values_from)
  y_values_all <- y_values_from
  
  num_of_steps <- length(y_values_max)
  
  xlim <- c(1, num_of_steps)
  if(is.null(ylim)){
    ylim_plot <- c(min(y_values_from), y_values_max[num_of_steps])
    if(type == "best"){
      ylim_plot[1] <- y_values_from[1] # for the "best" type this is the smallest point of the graph
    }
  }
  else
    ylim_plot <- ylim
  
  # make the plot stairs-like
  x_points <- c(1, rep(2:num_of_steps, each = 2))
  
  if(logarithmic_x){
    x_points <- log10(x_points)
    xlim <- log10(xlim)
  }
  
  graphics::plot.new()
  graphics::plot.window(xlim, ylim_plot)
  
  if(type != "best"){
    # make the plot stairs-like
    y_points <- c(rep(y_values_all[1:(length(y_values_all)-1)], each = 2),
                  y_values_all[length(y_values_all)])
    
    graphics::lines.default(x_points, y_points, type = "l",
                            lwd=3, col = color[1], # the first color
                            ...)
  }
  if(type != "all"){
    # make the plot stairs-like
    y_points <- c(rep(y_values_max[1:(length(y_values_max)-1)], each = 2),
                  y_values_max[length(y_values_max)])
    
    graphics::lines.default(x_points, y_points, lwd=3,
                            lty = 1,
                            col = color[length(color)], # the last color
                            ...)
  }
  
  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)
  
  if(show_legend){
    if(type == "both"){
      legend_text <- c("All calculated likelihoods",
                       "Maximum likelihoods calculated")
      lty <- c(1, 1)
      lwd <- c(3, 3)
    }else if(type == "all"){
      legend_text <- c("All calculated function values")
      lty <- 1
      lwd <- 3
    }else if(type == "best"){
      legend_text <- c("Maximum function values calculated")
      lty <- 1
      lwd <- 3
    }
    
    graphics::legend("bottomright", inset=.002,
                     legend = legend_text,
                     col = color, lty = lty, cex = 0.7,
                     box.lty=0, lwd = lwd)
  }
  
  invisible(NULL)
}


# TODO(summary, the n_0 of best perms, the distribution of likelihood values)
summary.gips <- function(object, ...){
  invisible(NULL)
}




