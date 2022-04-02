#' Project matrix
#'
#' Project matrix on the space of symmetrical matrices invariant
#' by a cyclic group of permutations.
#'
#' @param U matrix to be projected
#' @param perm permutation. Generator of a permutation group
#'
#' @return projected matrix
#' @export

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
#' all(sapply(out, function(v) all.equal(matrix_symvariant[v]))) == TRUE
#'
#' @noRd
get_equal_indices_by_perm <- function(perm, perm_size){
    # We are essentially looking for cycles of permutation defined on elements
    # of matrix (perm_size^2 -> perm_size^2)
    permuted_indices <- as.integer(permutations::as.word((perm)))
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
    is_subcycle_symmetrical <- is.null(subcycle_representatives)
    nonsymmetrical_subcycle_representatives <- get_double_from_single_indices(
        subcycle_representatives[!is_subcycle_symmetrical]
    )
    lower_triangle_representatives <-
       nonsymmetrical_subcycle_representatives[nonsymmetrical_subcycle_representatives[,1] >
                                                   nonsymmetrical_subcycle_representatives[,2]]
    upper_triangle_counterparts <- lower_triangle_representatives[,2:1]
    join_pairs <- matrix(c(get_single_from_double_indices(lower_triangle_representatives, perm_size),
                         get_single_from_double_indices(upper_triangle_counterparts, perm_size)),
                         ncol=2)
    joined_subcycles <- lapply(1:sum(!is_subcycle_symmetrical), function(i){
        pair <- join_pairs[i,]
        c(subcycles[!is_subcycle_symmetrical][[pair[1]]],
          subcycles[!is_subcycle_symmetrical][[pair[2]]])
    })
    append(joined_subcycles, subcycles[is_subcycle_symmetrical])
}

#' Get representative based on diagonal
#'
#' Select index corresponding to place in matrix, which is
#' a) closest to main diagonal
#' b) if a) equal-> closest to (1,1)
#'
#' @param indices integer vector interpreted as SINGLE indices of matrix with
#' matix_size rows columns
#'
#' @return Either single integer or NULL if places corresponding to them are
#' placed in symmetrical way (when there is an index on main diagonal
#' or when selection above does not yield single index)
#' @noRd

get_diagonal_representative <- function(indices, matrix_size){
    double_indices <- get_double_from_single_indices(indices, matrix_size)
    which_diag <- abs(double_indices[,1] - double_indices[,2])
    if (min(which_diag) == 0) return(NULL)
    considered_indices <- which(which_diag == min(which_diag))
    diag_positions <- apply(double_indices[considered_indices,], 2, min)
    representative_index <- which(diag_positions == min(diag_positions))
    if (length(repsesentative_index) > 1){
        return(NULL)
    }
    indices[considered_indices,][representative_index]
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


