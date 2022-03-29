#' Get Structure Constants
#'
#' Get Structure Constants from theorem 5 of the paper
#'
#' @param perm an element of class "cycle"
#'
#' @return list of 4 items: `r`, `d`, `k`, `L` - vectors of constants from theorem 1
#' @export
#' 
#' @examples
#' perm <- permutations::as.cycle(permutations::as.word(c(1,2,3,5,4)))
#' get_structure_constants(perm)

get_structure_constants <- function(perm){
    l <- get_cycle_represesentatives_and_lenghts(perm)
    representatives <- l[['representatives']]
    cycle_lenghts <- l[['cycle_lengths']]

    perm_order <- permutations::permorder(perm)
    M <- floor(perm_order / 2)
    # TBC
}

#' Get cycle representatives and lengths
#'
#' Essentially get iC, pC from paper
#'
#' @param perm an element of class "cycle"
#'
#' @return list with 2 items: `representatives` and `cycle_lengths`
#'

get_cycle_represesentatives_and_lenghts <-function(perm){
    representatives <- permutations::get1(perm)
    cycle_lengths <- sapply(perm[[1]], length)

    # unfortunately, above functions ignore cycles of length 1
    fixed_elements <- which(permutations::fixed(perm))
    elements_before_fixed <- fixed_elements - 1
    fixed_elements_with_nonfixed_predecessor <- fixed_elements[!elements_before_fixed %in% fixed_elements]

    # I know, i know
    for(el in fixed_elements){
        if(el == 1){
            representatives <- c(el, representatives)
            cycle_lengths <- c(1, cycle_lengths)
            next
        }
        if(el %in% fixed_elements_with_nonfixed_predecessor){
            previous_representative <- permutations::get1(permutations::get_cyc(el-1))
            index <- which(representatives == previous_representative)
            representatives <- append(representatives, el, index)
            cycle_lengths <- append(cycle_lengths, 1, index)
            next
        }
        index <- which(representatives == el-1)
        representatives <- append(representatives, el, index)
        cycle_lengths <- append(cycle_lengths, 1, index)
    }
    list('representatives' = representatives,
         'cycle_lenghts' = cycle_lengths)
}

