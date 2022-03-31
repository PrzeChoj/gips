#' Get Structure Constants
#'
#' Get Structure Constants from theorem 5 of the paper
#'
#' @param perm an element of class "cycle"
#' @param perm_size size of permutation
#'
#' @return list of 4 items: `r`, `d`, `k`, `L` - vectors of constants from theorem 1
#' @export
#'
#' @examples
#' perm <- permutations::as.cycle(permutations::as.word(c(1,2,3,5,4)))
#' get_structure_constants(perm)

get_structure_constants <- function(perm, perm_size) {
    l <- get_cycle_representatives_and_lengths(perm, perm_size)
    representatives <- l[['representatives']]
    cycle_lenghts <- l[['cycle_lengths']]
    perm_order <- permutations::permorder(perm)

    r <- calculate_r(cycle_lenghts, perm_order)
    d <- calculate_d(perm_order)

    L <- sum(r > 0)
    d <- d[r>0]
    r <- r[r>0]
    k <- d
    list ('r'=r,
          'd'=d,
          'k'=k,
          'L'=L)
}

#' Get cycle representatives and lengths
#'
#' Essentially get iC, pC from paper
#'
#' @param perm an element of class "cycle". Can't be an identity.
#' @param perm_size size of permutation
#'
#' @return list with 2 items: `representatives` and `cycle_lengths`
#' @noRd

get_cycle_representatives_and_lengths <- function(perm, perm_size) {
    if (permutations::is.id(perm)) {
        return(list('representatives'=1:perm_size,
                    'cycle_lengths'=rep(1, perm_size)))
    }
    representatives <- permutations::get1(perm)
    cycle_lengths <- sapply(perm[[1]], length)

    # unfortunately, above functions ignore cycles of length 1
    fixed_boolean <- permutations::fixed(perm)
    # if "last" elements are fixed, they are not returned.
    # correct for that
    if(length(fixed_boolean) < perm_size){
        fixed_boolean[(length(fixed_boolean) + 1):perm_size] <- TRUE
    }
    fixed_elements <- which(fixed_boolean)
    elements_before_fixed <- fixed_elements - 1
    fixed_elements_with_nonfixed_predecessor <-
        fixed_elements[!elements_before_fixed %in% fixed_elements]

    # I know, i know
    for (el in fixed_elements) {
        if (el == 1) {
            representatives <- c(el, representatives)
            cycle_lengths <- c(1, cycle_lengths)
            next
        }
        if (el %in% fixed_elements_with_nonfixed_predecessor) {
            previous_representative <-
                permutations::get1(permutations::get_cyc(perm, el - 1))
            index <-
                which(representatives == previous_representative)
            representatives <- append(representatives, el, index)
            cycle_lengths <- append(cycle_lengths, 1, index)
            next
        }
        index <- which(representatives == el - 1)
        representatives <- append(representatives, el, index)
        cycle_lengths <- append(cycle_lengths, 1, index)
    }
    list('representatives' = representatives,
         'cycle_lengths' = cycle_lengths)
}

#' Calculate structure constant r
#'
#' @return integer vector. Structure constant r WITH elements equal to 0.
#' @noRd
calculate_r <- function(cycle_lengths, perm_order) {
    M <- floor(perm_order / 2)
    if (M == 0) {
        # identity function
        return(length(cycle_lengths))
    }
    # for a in 0,1,...,floor(perm_order/2)
    # r_a = #{1:C such that a*p_c is a multiple of N}
    # AKA a*p_c %% N == 0
    # Corollary: N %% p_c == 0 for each p_c
    raw_alphas <- perm_order / cycle_lengths
    filtered_alphas <- raw_alphas[raw_alphas <= M]
    alpha_count <- table(filtered_alphas)

    r <- rep(0, M)
    r[as.integer(names(alpha_count))] <- as.integer(alpha_count)
    # correct for alpha == 0
    r <- c(length(cycle_lengths), r)
    r
}

#' Calculate structure constant d
#'
#' @noRd
calculate_d <- function(perm_order) {
    M <- floor(perm_order / 2)
    d <- c(1, rep(2, M))
    if (perm_order %% 2 == 0)
        d[M + 1] <- 1
    d
}
