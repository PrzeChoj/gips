#' Permutation object
#'
#' Create permutation objects to be passed to functions of `gips` package.
#'
#' @param x an object created with `permutations` package, or any object that
#' can be coerced using \code{\link[permutations]{permutation}} function.
#' @param size integer. Size of permutation (AKA cardinality of set, on which permutation
#' is defined).
#'
#' @seealso
#' \code{\link[permutations]{permutation}}
#'
#' @examples
#' gperm <- gips_perm(permutations::as.word(c(1,2,3,5,4)), 5)
#' gperm <- gips_perm(permutations::as.cycle('(5,4)'), 5)
#' # note the necessity of `size` parameter
#' gperm <- gips_perm(permutations::as.cycle('(5,4)'), 7)
#' gperm <- gips_perm('(1,2)(5,4)', 7)
#'
#' @export

gips_perm <- function(x, size){
    if(!inherits(x, 'permutation'))
        x <- permutations::permutation(x)
    if(rlang::is_missing(size))
        rlang::abort(c("There was a problem identified with provided argument:",
                       "i" = '`size` argument must be provided.',
                       "x" = "You did not provide the `size` argument."))
    if(!is.wholenumber(size))
        rlang::abort(c("There was a problem identified with provided argument:",
                       "i" = '`size` must be a whole number.',
                       "x" = paste0("You provided `size` == ", size, ".")))
    x <- permutations::as.cycle(x)

    if(length(unclass(x)) > 1){
        rlang::warn("Passing multiple permutations passed to `gips_perm` is not supported. Taking only the first one",
                    "i" = paste0("Passed ", length(unclass(x)), "permutations."))
        x <- x[1]
    }

    if(is.null(permutations::is.id(x))){
        return(new_gips_perm(list(), 0))
    }


    if (permutations::is.id(x)) {
        x <- as.list(1:size)
        return(new_gips_perm(x, size))
    }

    cycles <- unclass(x)[[1]]

    representatives <- permutations::get1(x)

    # unfortunately, cycles of length 1 are ignored
    fixed_boolean <- permutations::fixed(x)
    # if "last" elements are fixed, they are not returned.
    # correct for that
    if(length(fixed_boolean) < size){
        fixed_boolean[(length(fixed_boolean) + 1):size] <- TRUE
    }
    fixed_elements <- which(fixed_boolean)

    subcycles <- c(cycles, as.list(fixed_elements))

    new_gips_perm(subcycles, size)
}

#' @describeIn gips_perm Constructor
#'
#' Only intended for low-level use.
#'
#' @param cycles list of integer vectors. Each vector corresponds to a single cycle
#' of a permutation.
#'
#' @export

new_gips_perm <- function(cycles, size){
    rearranged_cycles <- lapply(cycles, rearrange_vector)
    representatives <- sapply(rearranged_cycles, function(v)v[1])
    reordered_cycles <- rearranged_cycles[order(representatives)]
    structure(reordered_cycles, size=size, class='gips_perm')
}

#' Print gips_perm
#'
#' Implementation of S3 method.
#'
#' @param x `gips_perm`
#' @param ... further arguments passed to \code{\link[permutations]{print.cycle}}
#'
#' @export
print.gips_perm <- function(x, ...){
    x <- permutations::as.cycle(x)
    permutations::print.cycle(x, ...)
}

#' Coerce gips_perm to character vector
#'
#' Implementation of S3 method.
#'
#' @param x `gips_perm`
#' @param ... further arguments passed to \code{\link[base]{as.character}}
#'
#' @export
as.character.gips_perm <- function(x, ...)
    as.character(permutations::as.cycle(x), ...)

#' Compose permutation with transposition
#'
#' @param gips_perm Object of `gips_perm` class
#' @param transposition integer vector of length 2. Transposition in a form of
#' cycle.
#'
#' @return `gips_perm` object. Composition of `gips_perm` parameter and `transposition`.
#'
#' @noRd
#' @examples
#' perm <- permutations::as.cycle('(1,2,3)(4,5)')
#' gperm <- gips_perm(perm, 6)
#' tr <- c(2,3)
#' tr_perm <- permutations::as.cycle(tr)
#'
#' composed <- compose_with_transposition(gperm, tr)
#' composed2 <- perm * tr_perm
#'
#' # composed and composed 2 refer to the same permutation
#'
compose_with_transposition <- function(gips_perm, transposition){
    cycle_1_index <- which(sapply(gips_perm, function(cycle)
        transposition[1] %in% cycle))
    cycle_2_index <- which(sapply(gips_perm, function(cycle)
        transposition[2] %in% cycle))
    cycle_1 <- gips_perm[[cycle_1_index]]
    cycle_2 <- gips_perm[[cycle_2_index]]
    composed_gips_perm <- gips_perm[c(-cycle_1_index, -cycle_2_index)]
    if(cycle_1_index == cycle_2_index){
        # We are breaking cycle into 2 cycles
        shifted_cycle <- shift_vector(cycle_1, which(cycle_1 == transposition[1])-1)
        new_cycle_1 <- shifted_cycle[1:(which(shifted_cycle == transposition[2])-1)]
        new_cycle_2 <- shifted_cycle[(which(shifted_cycle == transposition[2])):length(shifted_cycle)]
        composed_gips_perm <- add_cycle(composed_gips_perm, new_cycle_1)
        composed_gips_perm <- add_cycle(composed_gips_perm, new_cycle_2)
    } else {
        # We are merging 2 cycles
        ind <- which(cycle_1 == transposition[1])
        fragment_1 <- shift_vector(cycle_2, which(cycle_2 == transposition[2])-1)
        fragment_2 <- shift_vector(cycle_1, which(cycle_1 == transposition[1])-1)
        new_cycle <- c(fragment_1, fragment_2)
        composed_gips_perm <- add_cycle(composed_gips_perm, new_cycle)
    }
    new_gips_perm(composed_gips_perm, attr(gips_perm, 'size'))
}

#' Add a new cycle to permutation
#'
#' @param cycles list of integer vectors. Each corresponds to cycles of a permutation
#' @param new_cycle integer vector. None of its elements are present in `cycles`
#'
#' @noRd

add_cycle <- function(cycles, new_cycle){
    # Assume, that cycles are sorted by their min element
    # new_cycle - not necessarily
    new_cycle <- rearrange_vector(new_cycle)
    min_representatives <- sapply(cycles, function(v)v[1])
    insert_index <- findInterval(new_cycle[1], min_representatives)
    append(cycles, list(new_cycle), after=insert_index)
}
