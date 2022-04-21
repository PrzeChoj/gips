#' Project matrix
#'
#' Project matrix on the space of symmetrical matrices invariant
#' by a cyclic group of permutations.
#'
#' @param U matrix to be projected
#' @param perm permutation. Generator of a permutation group
#' @param perm_size size of permutation
#'
#' @return projected matrix
#' @export
#' 
#' @examples project_matrix(U = matrix(rnorm(49), nrow = 7),
#'                          perm = permutations::as.cycle(permutations::as.word(c(4,3,2,1,5))),
#'                          perm_size = 7)
project_matrix <- function(U, perm, perm_size){
    equal_indices_by_perm <- get_equal_indices_by_perm(perm, perm_size)
    mean_values <- sapply(equal_indices_by_perm, function(indices)
        mean(U[indices]))
    means <- rep(mean_values, sapply(equal_indices_by_perm, length))
    projected_matrix <- matrix(nrow=nrow(U), ncol=ncol(U))
    projected_matrix[unlist(equal_indices_by_perm)] <- means
    projected_matrix
}

#' Get indices of elements of perm_size x perm_size matrix, which should be equal
#'
#' @return a list of integer vectors. Each vector contatins SINGLE indices
#' of elements, which should be equal in symmetrical matrix invariant
#' by permutation
#'
#' @examples
#' perm <- permutations::as.cycle(permutations::as.word(c(2,3,1,5,4,6)))
#' matrix_symvariant <- matrix(c(
#' 2, 1, 1, 3, 3, 4,
#' 1, 2, 1, 3, 3, 4,
#' 1, 1, 2, 3, 3, 4,
#' 3, 3, 3, 5, 6, 7,
#' 3, 3, 3, 6, 5, 7,
#' 4, 4, 4, 7, 7, 8
#' ), byrow=TRUE, ncol=6)
#' out <- get_equal_indices_by_perm(perm, 6)
#' all(sapply(out, function(v) all.equal(matrix_symvariant[v]))) # TRUE
#'
#' @noRd
get_equal_indices_by_perm <- function(perm, perm_size){
    # We are essentially looking for cycles of permutation defined on elements
    # of matrix (perm_size^2 -> perm_size^2)
    permuted_indices <- as.integer(permutations::as.word((perm)))
    # correct for last fixed elements
    l <- length(permuted_indices)
    if (l < perm_size){
        permuted_indices[(l+1):perm_size] <- (l+1):perm_size
    }
    original_matrix_elements <- matrix(1:(perm_size^2), nrow=perm_size)
    permuted_matrix_elements <- original_matrix_elements[permuted_indices,][,permuted_indices]
    large_perm <- permutations::as.cycle(permutations::as.word(as.integer(permuted_matrix_elements)))
    subcycles <- get_subcycles(large_perm, perm_size^2)

    # Correct for symmetry
    subcycle_representatives <- sapply(subcycles, function(cyc){
        get_diagonal_representative(cyc, perm_size)
    })
    is_subcycle_symmetrical <- is.na(subcycle_representatives)
    if (all(is_subcycle_symmetrical)){
        return(subcycles)
    }
    nonsymmetric_subcycles <- subcycles[!is_subcycle_symmetrical]
    merge_pairs <- which_subcycles_merge(subcycle_representatives[!is_subcycle_symmetrical],
                                         perm_size)
    merged_subcycles <- lapply(1:nrow(merge_pairs), function(i){
        pair <- merge_pairs[i,]
        c(nonsymmetric_subcycles[[pair[1]]],
          nonsymmetric_subcycles[[pair[2]]])
    })
    append(merged_subcycles, subcycles[is_subcycle_symmetrical])
}

#' Which subcycles should be merged
#'
#' @param subcycle_representatives vector of indices, that point to equal values in
#' a matrix invariant by a permutation. Each element has its symmetrical
#' counterpart also in that vector.
#'
#' @return matrix with two columns, each row are indices of subcycles to be
#' merged. Those pairs are symmetrical with regard to main diagonal.
#' @noRd

which_subcycles_merge <- function(subcycle_representatives, perm_size){
    double_indices_of_representatives <- get_double_from_single_indices(
        subcycle_representatives,
        perm_size
    )
    lower_triangle_representatives <-
        double_indices_of_representatives[double_indices_of_representatives[,1] >
                                                    double_indices_of_representatives[,2],,drop=FALSE]
    upper_triangle_counterparts <- lower_triangle_representatives[,2:1,drop=FALSE]
    representative_pairs <- matrix(c(
        get_single_from_double_indices(lower_triangle_representatives, perm_size),
        get_single_from_double_indices(upper_triangle_counterparts, perm_size)),ncol=2)
    matrix(match(representative_pairs, subcycle_representatives),ncol=2)
}

#' Get representative based on diagonal
#'
#' Select index corresponding to place in matrix, which is
#' a) closest to main diagonal
#' b) if a) equal-> closest to (1,1).
#' Essentially two cycles are symmetrical to each other with regard to main
#' diagonal iff their indices are symmetrical.
#'
#' @param indices integer vector interpreted as SINGLE indices of matrix
#' @param matrix_size number of rows of square matrix
#'
#' @return Either single integer or NA if places corresponding to them are
#' placed in symmetrical way (when there is an index on main diagonal
#' or when selection above does not yield single index)
#' @noRd

get_diagonal_representative <- function(indices, matrix_size) {
    double_indices <-
        get_double_from_single_indices(indices, matrix_size)
    which_diag <- abs(double_indices[, 1] - double_indices[, 2])
    if (min(which_diag) == 0)
        return(NA)
    indices_on_middlest_diag <- which(which_diag == min(which_diag))
    distances_from_lupper_corner <-
        apply(double_indices[indices_on_middlest_diag, , drop = FALSE], 1, min)
    closest_to_lupper_corner <-
        which(distances_from_lupper_corner == min(distances_from_lupper_corner))
    if (length(closest_to_lupper_corner) > 1) {
        return(NA)
    }
    indices[indices_on_middlest_diag][closest_to_lupper_corner]
}

get_double_from_single_indices <- function(indices, matrix_size){
    row_indices <- indices %% matrix_size
    row_indices[row_indices == 0] <- matrix_size
    matrix(c(row_indices,
             ceiling(indices / matrix_size)), ncol=2)
}

get_single_from_double_indices <- function(indices, matrix_size){
    (indices[,2] - 1) * matrix_size + indices[,1]
}


