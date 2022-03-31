source('setup_example.R')

matrix_invariant_by_example_perm <- matrix(c(
    2, 1, 1, 3, 3, 4,
    1, 2, 1, 3, 3, 4,
    1, 1, 2, 3, 3, 4,
    3, 3, 3, 5, 6, 7,
    3, 3, 3, 6, 5, 7,
    4, 4, 4, 7, 7, 8
), byrow=TRUE, ncol=6)

test_that('orthogonal matrix fulfills condition from example 6',{
    U_Gamma <- prepare_orthogonal_matrix(example_perm, 6, example_basis)
    tested_matrix <- t(U_Gamma) %*% matrix_invariant_by_example_perm %*% U_Gamma
    tested_elements <- tested_matrix
    tested_elements[1:3, 1:3] <- 0
    tested_elements[4,4] <- 0
    tested_elements[5,5] <- 0
    tested_elements[6,6] <- 0
    expect_equal(tested_elements,
                 matrix(rep(0, 6*6), ncol=6))
    expect_true(isSymmetric(tested_matrix[1:3, 1:3]))
})

test_that('prepare_orthogonal_matrix works for example', {
    expect_equal(
        prepare_orthogonal_matrix(example_perm, 6, example_basis),
        example_orth_matrix
    )
})

test_that('get_v_matrix_for_subcycle works for 3-length cycle', {
    expect_equal(
        get_v_matrix_for_subcycle(
            example_perm,
            example_cycle_representatives[1],
            example_cycle_lengths[1],
            example_basis
        ),
        example_v_object[[1]]
    )})
test_that('get_v_matrix_for_subcycle works for 2-length cycle', {
    expect_equal(
        get_v_matrix_for_subcycle(
            example_perm,
            example_cycle_representatives[2],
            example_cycle_lengths[2],
            example_basis
        ),
        example_v_object[[2]]
    )})
test_that('get_v_matrix_for_subcycle works for identity', {
    expect_equal(
        get_v_matrix_for_subcycle(
            example_perm,
            example_cycle_representatives[3],
            example_cycle_lengths[3],
            example_basis
        ),
        example_v_object[[3]]
    )})

test_that('arrange_v_object works for example', {
    expect_equal(arrange_v_object(example_v_object),
                 example_orth_matrix)
})
