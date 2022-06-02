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
