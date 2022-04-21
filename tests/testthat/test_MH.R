
block_diagonal_matrix <- matrix(c(
    2,1,0,0,0,
    1,2,0,0,0,
    0,0,2,0,0,
    0,0,0,5,0,
    0,0,0,0,0.2
), byrow=TRUE, ncol=5)


test_that('calculate_block_determinants works', {
    expect_equal(calculate_determinants_of_block_matrices(block_diagonal_matrix,
                                                       c(2,3,5)),
                 c(3,2,1))
})
