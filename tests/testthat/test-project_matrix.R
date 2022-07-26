set.seed(2022)
values <- unique(as.integer(matrix_invariant_by_example_perm))
expected_equal_indices_by_example_perm <- lapply(values, function(v)
    which(as.integer(matrix_invariant_by_example_perm) == v))
full_double_indices <- matrix(c(rep(1:4, times=4),
                                rep(1:4, each=4)),
                              ncol=2)
S <- matrix(rnorm(36), ncol=6)

gips_example_perm <- gips_perm(example_perm, 6)


test_that('projected matrix is symvariant by example_perm', {
    projected <- project_matrix(S, gips_example_perm)
    expect_true(isSymmetric(projected))
    # now we need to check only lower diagonal
    expect_equal(rep(projected[1,1], 2),
                 c(projected[2,2], projected[3,3]))
    expect_equal(rep(projected[2,1],2),
                 projected[3, 1:2])
    expect_equal(rep(projected[4,1],6),
                 as.vector(projected[4:5,1:3]))
    expect_equal(rep(projected[6,1],2),
                 projected[6, 2:3])
    expect_equal(projected[4,4], projected[5,5])
    expect_equal(projected[6,4], projected[6,5])
})

test_that('get_equal_indices_by_perm works for example_perm', {
    actual_l <- lapply(get_equal_indices_by_perm(gips_example_perm),
                       sort)
    expected_l <- lapply(expected_equal_indices_by_example_perm,
                         sort)
    expect_setequal(actual_l, expected_l)
})

test_that('get_equal_indices_by_perm works for identity', {
    expect_setequal(get_equal_indices_by_perm(gips_perm(to_perm(1:3), 3)),
                 list(1, c(2,4), c(3,7), 5, c(6,8), 9))
})

test_that('get_single_from_double_indices works',{
    expect_equal(get_single_from_double_indices(full_double_indices, 4),
                 1:16)
})

test_that('get_single_from_double_indices works for 0 input',{
    expect_equal(get_single_from_double_indices(matrix(numeric(0), ncol=2), 4),
                 numeric(0))
})

test_that('get_double_from_single_indices works',{
    expect_equal(get_double_from_single_indices(1:16, 4),
                 full_double_indices)
    expect_equal(get_double_from_single_indices(as.vector(matrix(1:16, ncol=4)), 4),
                 full_double_indices)
})

test_that('get_double_from_single_indices works for 0 input',{
    expect_equal(get_double_from_single_indices(numeric(0), 4),
                 matrix(numeric(0), ncol=2))
})


