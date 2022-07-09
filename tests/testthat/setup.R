to_perm <- function(v) permutations::as.cycle(permutations::as.word(v))
example_perm <- to_perm(c(2,3,1,5,4,6))
example_perm2 <- to_perm(c(2,3,4,5,1))
example_cycle_representatives <- c(1,4,6)
example_cycle_lengths <- c(3,2,1)
example_basis <- diag(6)
example_structure_constants <- list('r'=c(3,1,1),
                                    'd'=c(1,2,1),
                                    'k'=c(1,2,1),
                                    'L'=3,
                                    'dim_omega'=c(6,1,1))
example_structure_constants2 <- list('r'=c(1,1,1),
                                     'd'=c(1,2,2),
                                     'k'=c(1,2,2),
                                     'L'=3,
                                     'dim_omega'=c(1,1,1))
example_v_object <- list(matrix(c(rep(1/sqrt(3), 3), rep(0, 3),
                                  sqrt(2/3), rep(-sqrt(1/6),2), rep(0,3),
                                  0, 1/sqrt(2), -1/sqrt(2), rep(0, 3)), nrow=6),
                         matrix(c(rep(0,3), 1/sqrt(2), 1/sqrt(2), 0,
                                  rep(0,3), 1/sqrt(2), -1/sqrt(2), 0),nrow=6),
                         matrix(c(rep(0, 5), 1), nrow=6))

example_orth_matrix <- matrix(c(rep(1/sqrt(3), 3), rep(0, 3),
                                rep(0, 3), rep(1/sqrt(2), 2), 0,
                                rep(0, 5), 1,
                                sqrt(2/3), rep(-sqrt(1/6),2), rep(0,3),
                                0, 1/sqrt(2), -1/sqrt(2), rep(0, 3),
                                rep(0,3),1/sqrt(2), -1/sqrt(2), 0), nrow=6)

matrix_invariant_by_example_perm <- matrix(c(
    2, 1, 1, 3, 3, 4,
    1, 2, 1, 3, 3, 4,
    1, 1, 2, 3, 3, 4,
    3, 3, 3, 5, 6, 7,
    3, 3, 3, 6, 5, 7,
    4, 4, 4, 7, 7, 8
), byrow=TRUE, ncol=6)

number_of_observations <- 13
