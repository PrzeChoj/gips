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
    projected_matrix <- matrix(
        means[unlist(equal_indices_by_perm)], nrow=nrow(U)
    )# that's not it. WIP
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
    get_subcycles(large_perm, perm_size^2)
}
