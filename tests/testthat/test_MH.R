block_diagonal_matrix <- matrix(c(
    2,1,0,0,0,
    1,2,0,0,0,
    0,0,2,0,0,
    0,0,0,5,0,
    0,0,0,0,0.2
), byrow=TRUE, ncol=5)

c <- permutations::as.cycle(permutations::as.word(c(2,1)))
cprim <- permutations::as.cycle(permutations::as.word(c(1,2)))

U1 <- matrix(c(1,0.5,0.5,2), nrow=2,byrow = TRUE)
U2 <- matrix(c(1.5,0.5,0.5,1.5), nrow=2,byrow = TRUE)


test_that('calculate_block_determinants works', {
    expect_equal(calculate_determinants_of_block_matrices(block_diagonal_matrix,
                                                          c(2,3,5)),
                 c(3,2,1))
})

test_that('goal_function returns proper values', {
  skip("Fix not yet implemented")  # PC: I think the variables Dc_exponent and DcUc_exponent in goal_function should be swapped
  
  expect_equal(goal_function(c, 2, 100, U1),
               6^(-103/2) * gamma(103/2) * gamma(103/2) / (pi / 4))
  expect_equal(goal_function(cprim, 2, 100, U1),
               (23/4)^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
})








