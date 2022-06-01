gips_perm <- function(x, size){
    if(!inherits(x, 'permutation'))
        x <- permutations::permutation(x)
    x <- permutations::as.cycle(x)

    if (permutations::is.id(x)) {
        x <- as.list(1:size)
        return(structure(x, size=size, class='gips_perm'))
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

    # Sort for now - cause it's simple
    # could be replace by merge-like function (used i.e. in mergesort)
    # use `findInterval`
    representatives <- sapply(subcycles, min)
    structure(subcycles[order(representatives)], size=size, class='gips_perm')
}
