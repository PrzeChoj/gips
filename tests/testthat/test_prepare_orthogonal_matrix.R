source('setup_example.R')

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
