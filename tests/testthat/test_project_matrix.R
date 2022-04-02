
values <- unique(as.integer(matrix_invariant_by_example_perm))
expected_equal_indices_by_example_perm <- lapply(values, function(v)
    which(as.integer(matrix_invariant_by_example_perm) == v))
full_double_indices <- matrix(c(rep(1:4, times=4),
                                rep(1:4, each=4)),
                              ncol=2)


test_that('get_equal_indices_by_perm works for example_perm', {
    expect_equal(get_equal_indices_by_perm(example_perm, 6),
                 expected_equal_indices_by_example_perm)
})


test_that('are_matrix_indices_symmetrical works for diagonal', {
    expect_true(are_matrix_indices_symmetrical(c(1, 5, 9), 3))
})

test_that('are_matrix_indices_symmetrical works for symmetrical', {
    expect_true(are_matrix_indices_symmetrical(c(1, 3, 7), 3))
})

test_that('are_matrix_indices_symmetrical works for nonsymmetrical', {
    expect_true(!are_matrix_indices_symmetrical(c(1, 3, 6), 3))
})

test_that('get_single_from_double_indices works',{
    expect_equal(get_single_from_double_indices(full_double_indices, 4),
                 1:16)
})

test_that('get_double_from_single_indices works',{
    expect_equal(get_double_from_single_indices(1:16, 4),
                 full_double_indices)
    expect_equal(get_double_from_single_indices(as.vector(matrix(1:16, ncol=4)), 4),
                 full_double_indices)
})
