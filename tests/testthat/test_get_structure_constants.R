test_that('get_structure_constants for examples from paper',{
    expect_equal(get_structure_constants(example_perm, 6),
                 example_structure_constants)
    expect_equal(get_structure_constants(example_perm2, 5),
                 example_structure_constants2)
})

# Example 6 from the paper
test_that('calculate_r works for example from paper', {
    expect_equal(calculate_r(c(3,2,1), 6),
                 c(3,0,1,1))
})

test_that('calculate_r works for identity', {
    expect_equal(calculate_r(c(1,1,1,1), 1),
                 4)
})

test_that('calculate_d works for even perm_order', {
    expect_equal(calculate_d(6), c(1,2,2,1))
    expect_equal(calculate_d(4), c(1,2,1))
    expect_equal(calculate_d(2), c(1,1))
})

test_that('calculate_d works for odd perm_order', {
    expect_equal(calculate_d(5), c(1,2,2))
    expect_equal(calculate_d(3), c(1,2))
})

test_that('calculate_d works for identity', {
    expect_equal(calculate_d(1), 1)
})
