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
project_matrix <- function(U, perm, perm_size, precomputed_equal_indices=NULL){
    if(is.null(precomputed_equal_indices))
        equal_indices_by_perm <- get_equal_indices_by_perm(perm, perm_size)
    else
        equal_indices_by_perm <- precomputed_equal_indices
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
#' by permutation `perm`.
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
#' @noRd
get_equal_indices_by_perm <- function(perm, perm_size){
    # We are essentially looking for cycles of permutation defined on elements
    # of matrix (perm_size^2 -> perm_size^2)
    subcycles <- get_subcycles(perm, perm_size)

    # Let's go
    subcycle_indice_pairs <- matrix(c(rep(1:length(subcycles), each=length(subcycles)),
                                      rep(1:length(subcycles), times=length(subcycles))),
                                    ncol=2)
    subcycle_indice_pairs <- subcycle_indice_pairs[subcycle_indice_pairs[,1] <=
                                                       subcycle_indice_pairs[,2],,
                                                   drop=FALSE]

    nested_list <- lapply(1:nrow(subcycle_indice_pairs), function(pair_index){
        i <- subcycle_indice_pairs[pair_index, 1]
        j <- subcycle_indice_pairs[pair_index, 2]
        subcycle_1 <- subcycles[[i]]
        subcycle_2 <- subcycles[[j]]

        matrix_subcycle_length <- numbers::LCM(length(subcycle_1),
                                               length(subcycle_2))
        number_of_matrix_subcycles <- length(subcycle_1)*length(subcycle_2)/
            matrix_subcycle_length

        cycle_indices <- 1:matrix_subcycle_length
        subcycle_1_indices <- cycle_indices %% length(subcycle_1)
        subcycle_1_indices[subcycle_1_indices==0] <- length(subcycle_1)
        subcycle_1_elements <- subcycle_1[subcycle_1_indices]

        subcycle_2_indices <- cycle_indices %% length(subcycle_2)
        subcycle_2_indices[subcycle_2_indices==0] <- length(subcycle_2)

        lapply(1:(number_of_matrix_subcycles),function(k){
            subcycle_2_elements <- subcycle_2[shift(subcycle_2_indices, k-1)]

            double_indices <- matrix(c(
                subcycle_1_elements, subcycle_2_elements,
                subcycle_2_elements, subcycle_1_elements
            ), ncol=2)
            single_indices <- get_single_from_double_indices(double_indices,
                                                             perm_size)

            if(i == j && single_indices[matrix_subcycle_length+1] %in%
               single_indices[1:matrix_subcycle_length]){
                single_indices <- single_indices[1:matrix_subcycle_length]
            }
            shift_vector(single_indices, which.min(single_indices)-1)
        })
    })
    unlist(nested_list, recursive = FALSE)
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

    # contrdiag is a diagonal parallel to main contrdiag
    # main contrdiag is going from right-upper to left-bottom corner
    which_contrdiag <- double_indices[indices_on_middlest_diag, 1] +
        double_indices[indices_on_middlest_diag, 2]

    closest_to_lupper_corner <-
        which(which_contrdiag == min(which_contrdiag))
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


