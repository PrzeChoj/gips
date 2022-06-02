#' Project matrix
#'
#' Project matrix on the space of symmetrical matrices invariant
#' by a cyclic group of permutations.
#'
#' @param U matrix to be projected.
#' @param perm permutation. Generator of a permutation group.
#'             Either of `gips_perm` or `permutations::cycle` class.
#' @param precomputed_equal_indices used in internal calculations in case when the equal indices have already been calculated; If it is not the case, leave this parameter as \code{NULL} and those will be computed
#'
#' @return projected matrix
#' @export
#'
#' @examples
#' gperm <- gips_perm(permutations::as.word(c(4,3,2,1,5)), 7)
#' U <- matrix(rnorm(49), nrow = 7)
#' projected_U <- project_matrix(U, perm = gips_perm)
project_matrix <- function(U, perm, perm_size=NULL, precomputed_equal_indices=NULL){
    if(!is.matrix(U) || nrow(U) != ncol(U))
        rlang::abort('`U` must be a square matrix.')
    if(is.null(precomputed_equal_indices)){
        perm_size <- ncol(U)
        if(!inherits(perm, 'gips_perm')){
            perm <- gips_perm(perm, perm_size)
        } else if(attr(perm, 'size') != ncol(U))
            rlang::abort('Size of `perm` must be equal to number of columns and rows of `U`.')
        equal_indices_by_perm <- get_equal_indices_by_perm(perm)
    }
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
#' @param perm gips_perm
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
get_equal_indices_by_perm <- function(perm){
    perm_size <- attr(perm, 'size')
    # We'll be iterating over pairs of subcycles
    subcycle_indice_pairs <- matrix(c(rep(1:length(perm), each=length(perm)),
                                      rep(1:length(perm), times=length(perm))),
                                    ncol=2)
    subcycle_indice_pairs <- subcycle_indice_pairs[subcycle_indice_pairs[,1] <=
                                                       subcycle_indice_pairs[,2],,
                                                   drop=FALSE]
    # Let's go
    nested_list <- lapply(1:nrow(subcycle_indice_pairs), function(pair_index){
        i <- subcycle_indice_pairs[pair_index, 1]
        j <- subcycle_indice_pairs[pair_index, 2]
        subcycle_1 <- perm[[i]]
        subcycle_2 <- perm[[j]]

        # matrix_subcycle is a subcycle of permutation P defined as
        # P(k,l) = (perm(k), perm(l))
        matrix_subcycle_length <- numbers::LCM(length(subcycle_1),
                                               length(subcycle_2))
        number_of_matrix_subcycles <- length(subcycle_1)*length(subcycle_2)/
            matrix_subcycle_length

        # Instead of operating on subcycle elements, we will be operating
        # on elements' indices
        # I.e. for subcycle s=[3,5,2] we operate on indices [1,2,3]
        # s[1] = 3 etc
        elements_indices <- 1:matrix_subcycle_length
        subcycle_1_indices <- elements_indices %% length(subcycle_1)
        subcycle_1_indices[subcycle_1_indices==0] <- length(subcycle_1)
        subcycle_1_elements <- subcycle_1[subcycle_1_indices]

        subcycle_2_indices <- elements_indices %% length(subcycle_2)
        subcycle_2_indices[subcycle_2_indices==0] <- length(subcycle_2)

        lapply(1:(number_of_matrix_subcycles),function(k){
            subcycle_2_elements <- subcycle_2[shift_vector(subcycle_2_indices,
                                                           k-1)]

            double_indices <- matrix(c(
                subcycle_1_elements, subcycle_2_elements,
                subcycle_2_elements, subcycle_1_elements
            ), ncol=2)
            single_indices <- get_single_from_double_indices(double_indices,
                                                             perm_size)

            # We need to correct for matrix subcycles, that are symmetric
            if(i == j &&
               single_indices[matrix_subcycle_length+1] %in%
               single_indices[1:matrix_subcycle_length]){
                single_indices <- single_indices[1:matrix_subcycle_length]
            }
            shift_vector(single_indices, which.min(single_indices)-1)
        })
    })
    unlist(nested_list, recursive = FALSE)
}

#' Indices utils
#'
#' Elements of `matrix` can be accessed by double indices `M[i,j]`
#' or, when treating matrix as a long vector, single indices `M[k]`.
#' These functions allow to switch from one kind of indices to another.
#'
#' @noRd

get_double_from_single_indices <- function(indices, matrix_size){
    row_indices <- indices %% matrix_size
    row_indices[row_indices == 0] <- matrix_size
    matrix(c(row_indices,
             ceiling(indices / matrix_size)), ncol=2)
}

get_single_from_double_indices <- function(indices, matrix_size){
    (indices[,2] - 1) * matrix_size + indices[,1]
}


