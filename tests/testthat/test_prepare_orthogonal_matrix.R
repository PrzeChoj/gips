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
    expect_equal(tested_matrix[4,4], tested_matrix[5,5])
})

test_that('prepare_orthogonal_matrix works for example', {
    expect_equal(
        prepare_orthogonal_matrix(example_perm, 6, example_basis),
        example_orth_matrix
    )
})

expected_permuted_representatives <- matrix(c(
    1,4,6,
    2,5,6,
    3,4,6
), nrow=3, byrow=TRUE)

test_that('get_permuted_representatives works for example_perm', {
    expect_equal(
        get_permuted_representatives(
            example_perm,
            example_cycle_representatives,
            example_cycle_lengths
        ),
        expected_permuted_representatives
    )})

perm_with_1_cycle_with_fixed <- to_perm(c(2,1,3))
perm_with_1_cycle <- to_perm(c(2,3,1))

test_that('get_permuted_representatives works for 1-cycle perm', {
    expect_equal(
        get_permuted_representatives(
            perm_with_1_cycle_with_fixed,
            c(1,3),
            c(2,1)
        ),
        matrix(c(
            1,3,
            2,3
        ),nrow=2, byrow = TRUE)
    )
    expect_equal(
        get_permuted_representatives(
            perm_with_1_cycle,
            1,
            3
        ),
        matrix(c(
            1,
            2,
            3
        ),nrow=3, byrow = TRUE)
    )
})

test_that('get_v_matrix_for_subcycle works for 3-length cycle', {
    expect_equal(
        get_v_matrix_for_subcycle(
            expected_permuted_representatives[1:3,1],
            example_basis
        ),
        example_v_object[[1]]
    )})
test_that('get_v_matrix_for_subcycle works for 2-length cycle', {
    expect_equal(
        get_v_matrix_for_subcycle(
            expected_permuted_representatives[1:2,2],
            example_basis
        ),
        example_v_object[[2]]
    )})
test_that('get_v_matrix_for_subcycle works for identity', {
    expect_equal(
        get_v_matrix_for_subcycle(
            expected_permuted_representatives[1,3],
            example_basis
        ),
        example_v_object[[3]]
    )})



test_that('arrange_v_object works for example', {
    expect_equal(arrange_v_object(example_v_object),
                 example_orth_matrix)
})
