#' Shift vector
#'
#' Move k elements from the start of a vector to its end.
#'
#' @noRd

shift_vector <- function(v, k){
    if(k==0) return(v)
    c(v[-(1:k)], v[1:k])
}

#' Rearrange vector
#'
#' Move elements from the start of a vector to its end, so that the minimal
#' element will be first.
#'
#' @examples
#' v <- c(5,3,2,1,4)
#' rearranged <- rearrange_vector(v)
#' all(rearranged == c(1,4,5,3,2)) # TRUE
#'
#' @noRd

rearrange_vector <- function(v){
    shift_vector(v, which.min(v)-1)
}


#' Is matrix symmetric
#' 
#' We did not use the `matrixcalc::is.positive.semi.definite` function, because
#' here we b=have no checks(because they were done before) and the `tol`
#' argument is taken relative, not absolute. The `tol` argument is 1e-06,
#' the same as in the `MASS::mvrnorm` function.
#' 
#' @noRd
is.positive.semi.definite.matrix <- function (x, tol = 1e-06)
{
  ev <- eigen(x, symmetric = TRUE)$values
  
  return(all(ev >= -tol * abs(ev[1L])))
}