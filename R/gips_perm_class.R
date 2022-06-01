gips_perm <- function(x, size){
    if(!inherits(x, 'permutation'))
        x <- permutations::permutation(x)
    if(!is.wholenumber(size))
        rlang::abort('`size` must be an integer.')
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

is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

compose_with_transposition <- function(gips_perm, transposition){
    cycle_1_index <- which(sapply(gips_perm, function(cycle)
        transposition[1] %in% cycle))
    cycle_2_index <- which(sapply(gips_perm, function(cycle)
        transposition[2] %in% cycle))
    cycle_1 <- gips_perm[[cycle_1_index]]
    cycle_2 <- gips_perm[[cycle_2_index]]
    new_gips_perm <- gips_perm[c(-cycle_1_index, -cycle_2_index)]
    if(cycle_1_index == cycle_2_index){
        # We are breaking cycle into 2 cycles
        shifted_cycle <- shift_vector(cycle_1, which(cycle_1 == transposition[1])-1)
        new_cycle_1 <- shifted_cycle[1:(which(shifted_cycle == transposition[2])-1)]
        new_cycle_2 <- shifted_cycle[(which(shifted_cycle == transposition[2])):length(shifted_cycle)]
        new_gips_perm <- add_cycle(new_gips_perm, new_cycle_1)
        new_gips_perm <- add_cycle(new_gips_perm, new_cycle_2)
    } else {
        # We are merging 2 cycles
        ind <- which(cycle_1 == transposition[1])
        fragment_1 <- shift_vector(cycle_2, which(cycle_2 == transposition[2])-1)
        fragment_2 <- shift_vector(cycle_1, which(cycle_1 == transposition[1])-1)
        new_cycle <- c(fragment_1, fragment_2)
        new_gips_perm <- add_cycle(new_gips_perm, new_cycle)
    }
    structure(new_gips_perm, size=attr(gips_perm, 'size'), class='gips_perm')
}

add_cycle <- function(cycles, new_cycle){
    # Assume, that cycles are sorted by their min element
    # new_cycle - not necessarily
    new_cycle <- shift_vector(new_cycle, which.min(new_cycle)-1)
    min_representatives <- sapply(cycles, function(v)v[1])
    insert_index <- findInterval(new_cycle[1], min_representatives)
    append(cycles, list(new_cycle), after=insert_index)
}
