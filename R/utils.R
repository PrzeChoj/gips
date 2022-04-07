#' Get subcycles
#'
#' `permutations::cycle` returns permutation without cycles of length 1.
#' We need a function, that fixes that
#'
#' @param perm permutations::cycle object
#' 
#' @example get_subcycles(as.cycle(as.word(c(4,3,2,1,5))), 10)
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
    
    # A: Why can't we just:
    # return(c(cycles, as.list(fixed_elements))) #?
    
    elements_before_fixed <- fixed_elements - 1
    fixed_elements_with_nonfixed_predecessor <-
        fixed_elements[!elements_before_fixed %in% fixed_elements]

    # I know, i know
    for (el in fixed_elements) {
        if (el == 1) {
            representatives <- c(1, representatives)
            cycles <- append(cycles, 1, 0)
            next
        }
        if (el %in% fixed_elements_with_nonfixed_predecessor) {
            previous_representative <-
                permutations::get1(permutations::get_cyc(perm, el - 1))
            index <-
                which(representatives == previous_representative)
            cycles <- append(cycles, el, index)
            representatives <- append(representatives, el, index)
            next
        }
        index <- which(representatives == el - 1)
        representatives <- append(representatives, el, index)
        cycles <- append(cycles, el, index)
    }
    cycles
}
