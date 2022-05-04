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


test_that('goal_function returns proper values', {
  # The value of goal_function on matrix should the same as the projection of the matrix
  # and U2 == pi_c(U1)
  expect_equal(goal_function(c, 2, 100, U1),
               goal_function(c, 2, 100, U2))

  # Those values were calculated by hand:
  expect_equal(goal_function(c, 2, 100, U1),
               6^(-103/2) * gamma(103/2) * gamma(103/2) / (pi / 4))
  expect_equal(goal_function(cprim, 2, 100, U1),
               (23/4)^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
  expect_equal(goal_function(cprim, 2, 100, U2),
               6^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
})


test_that('goal_function has desired property', {
  skip("Not yet work; see issue#4")

  # Example from the paper chapter 5

  p <- 10
  n <- 20

  mu <- numeric(p)
  sigma <- matrix(numeric(p*p), nrow=p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i,j] <- 1 - abs(i-j) / p
    }
    sigma[i,i] <- 1 + 1/p
  }

  Z <- MASS::mvrnorm(n, mu = mu, Sigma = sigma)
  U <- t(Z) %*% Z

  actual_permutation <- permutations::as.cycle(permutations::as.word(c(2:p, 1)))
  actual_permutation_function_value <- goal_function(real_permutation,
                                                           p, n, U/n)
  another_permutation_function_value <- goal_function(permutations::id,
                                                      p, n, U/n)

  # Example from the paper's Table 4:
  #another_permutation2 <- permutations::as.cycle(permutations::as.word(c(6,7,5,8,9,2,1,10,3,4)))
  #another_permutation2_function_value <- goal_function(another_permutation2,
  #                                                     p, n, U/n)

  # We want the goal function to have a bigger value for the real permutation than for the another
  expect_gt(real_permutation_function_value,
            another_permutation_function_value)
})

#TODO
test_that('calculate phi_part works', {
    expect_true(TRUE)
})

test_that('calculate_block_determinants works', {
    expect_equal(calculate_determinants_of_block_matrices(block_diagonal_matrix,
                                                          c(2,3,5)),
                 c(3,2,1))
})



