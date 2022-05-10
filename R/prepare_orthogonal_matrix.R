#' Prepare orthogonal matrix
#'
#' Calculate orthogonal matrix U_Gamma
#' for decomposition in Theorem 1. To calculate it, we use Theorem 6.
#'
#' @param perm permuations::cycle object
#' @param perm_size size of permutation
#' @param basis matrix with basis vectors in COLUMNS. Identity by default
#' @return matrix p x p with columns from V object elements, sorted according to
#'     Theorem 6
#' @export
prepare_orthogonal_matrix <- function(perm, perm_size, basis=NULL){
    if(is.null(basis))
        basis <- diag(nrow=perm_size)
    l <- get_cycle_representatives_and_lengths(perm, perm_size)
    v_object <- get_v_object(perm, l[['representatives']], l[['cycle_lengths']],
                             basis)
    arrange_v_object(v_object)
}

#' Get V object defined before theorem 6
#'
#' @param basis matrix
#' @return list of length(cycle_lengths) length. ith element is a
#'    matrix p x cycle_lengths[i].
#' @noRd

get_v_object <- function(perm, cycle_representatives, cycle_lengths, basis){
    permuted_representatives <- get_permuted_representatives(perm, cycle_representatives,
                                                             cycle_lengths)
    v_object <- lapply(1:length(cycle_representatives), function(i){
        curr_cycle_length <- cycle_lengths[i]
        relevant_permutation_power_indices <- 1:curr_cycle_length
        curr_permuted_representatives <- permuted_representatives[relevant_permutation_power_indices,i]
        get_v_matrix_for_subcycle(curr_permuted_representatives, basis)
    })
    v_object
}

get_permuted_representatives <- function(perm, cycle_representatives, cycle_lengths){
    # perm ignores last fixed elements
    if(permutations::is.id(perm)){
        out <- matrix(1:max(cycle_representatives), nrow=1)
        return(out)
    }
    perm_size <- permutations::size(perm)
    last_fixed_representatives <- cycle_representatives[cycle_representatives > perm_size]
    other_representatives <- cycle_representatives[cycle_representatives <= perm_size]
    other_permuted_representatives <- as.function(perm ^ (1:max(cycle_lengths)-1))(other_representatives)
    if(!is.matrix(other_permuted_representatives)){
        other_permuted_representatives <- matrix(other_permuted_representatives,
                                                ncol=length(other_representatives))
    }
    if(length(last_fixed_representatives) == 0) return(other_permuted_representatives)
    last_permuted_representatives <- matrix(rep(last_fixed_representatives,
                                                each = nrow(other_permuted_representatives)),
                                            nrow = nrow(other_permuted_representatives))
    cbind(other_permuted_representatives, last_permuted_representatives)
}

get_v_matrix_for_subcycle <- function(permuted_representative, basis){
    cycle_length <- length(permuted_representative)
    k_s <- 1:cycle_length - 1
    chosen_basis_columns <- basis[,permuted_representative, drop=FALSE] # matrix p x curr_cycle_length

    first_element <- apply(chosen_basis_columns, 1, sum) /
        sqrt(cycle_length)
    v_matrix <- matrix(nrow=nrow(basis), ncol=cycle_length)
    v_matrix[,1] <- first_element

    max_beta <- floor((cycle_length - 1) / 2)
    if (max_beta >= 1) {
        betas <- 1:max_beta
        trygonometric_argument <- 2 * pi * outer(k_s, betas) / cycle_length
        even_elements <-
            chosen_basis_columns %*% cos(trygonometric_argument) *
            sqrt(2 / cycle_length)
        odd_elements <-
            chosen_basis_columns %*% sin(trygonometric_argument) *
            sqrt(2 / cycle_length)

        v_matrix[, 2 * betas] <- even_elements
        v_matrix[, 2 * betas + 1] <- odd_elements
    }
    if (cycle_length %% 2 == 0) {
        last_element <- chosen_basis_columns %*% cos(pi * k_s) /
            sqrt(cycle_length)
        v_matrix[, cycle_length] <- last_element
    }
    v_matrix
}

#' Execute multiple permutations on a single object
#' @noRd
permute_vectorized <- function(perms, x){
    cycles <- get_cyc(perms, x)
    sapply(unclass(cycles), function(l){
        cycle <- l[[1]]
        index <- which(cycle == x) + 1 %% length(cycle)
        cycle[index]
    })
}

#' Arrange V object
#'
#' @return matrix p x p with columns from V object elements, sorted according to
#'     Theorem 6
#' @noRd
arrange_v_object <- function(v_object){
    # we put v_k^v earlier than v_k'^c' if
    # i) [k/2]/p_c < [k'/2]/p_c' or
    # ii) [k/2]/p_c == [k'/2]/p_c' and c<c' or
    # iii) [k/2]/p_c == [k'/2]/p_c' and c == c' and k is even and k' is odd

    v_matrix <- do.call(cbind, v_object)

    features_list <- lapply(1:length(v_object), function(i){
        cycle_length <- ncol(v_object[[i]])
        f_1 <- floor(1:cycle_length / 2) / cycle_length
        f_2 <- rep(i, cycle_length)
        f_3 <- 1:cycle_length %% 2
        data.frame('f_1' = f_1,
                   'f_2' = f_2,
                   'f_3' = f_3)
    })

    df <- do.call(rbind, features_list)
    sorting_indices <- order(df$f_1, df$f_2, df$f_3)
    v_matrix[,sorting_indices]
}
