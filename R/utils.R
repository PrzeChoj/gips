#' Get subcycles
#'
#' `permutations::cycle` returns permutation without cycles of length 1.
#' We need a function, that fixes that
#'
#' @param perm permutations::cycle object
#'
#' @example get_subcycles(permutations::as.cycle(permutations::as.word(c(4,3,2,1,5))), 7)
#'
#' @return list of integer vectors - cycles INCLUDING cycles of length 1
#' @noRd
get_subcycles <- function(perm, perm_size){
    if (permutations::is.id(perm)) {
        return(as.list(1:perm_size))
    }
    cycles <- unclass(perm)[[1]]
    representatives <- permutations::get1(perm)

    # unfortunately, cycles of length 1 are ignored
    fixed_boolean <- permutations::fixed(perm)
    # if "last" elements are fixed, they are not returned.
    # correct for that
    if(length(fixed_boolean) < perm_size){
        fixed_boolean[(length(fixed_boolean) + 1):perm_size] <- TRUE
    }
    fixed_elements <- which(fixed_boolean)

    subcycles <- c(cycles, as.list(fixed_elements))

    # Sort for now - cause it's simple
    # could be replace by merge-like function (used i.e. in mergesort)
    # use `findInterval`
    representatives <- sapply(subcycles, min)
    subcycles[order(representatives)]
}

#' Shift vector
#'
#' Move k elements from the start of a vector to its end.
#'
#' @noRd

shift_vector <- function(v, k){
    if(k==0) return(v)
    c(v[-(1:k)], v[1:k])
}
