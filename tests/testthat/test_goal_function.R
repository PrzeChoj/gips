block_diagonal_matrix <- matrix(c(
    2,1,0,0,0,
    1,2,0,0,0,
    0,0,2,0,0,
    0,0,0,5,0,
    0,0,0,0,0.2
), byrow=TRUE, ncol=5)

c_perm <- permutations::as.cycle(permutations::as.word(c(2,1)))
id_perm <- permutations::id

U1 <- matrix(c(1,0.5,0.5,2), nrow=2, byrow = TRUE)
U2 <- matrix(c(1.5,0.5,0.5,1.5), nrow=2, byrow = TRUE)
D_matrix <- diag(nrow = 2)

test_that('goal_function returns proper values', {
  # The value of goal_function on matrix should the same as the projection of the matrix
  # and U2 == pi_c(U1)
  expect_equal(goal_function(c_perm, 100, U1, D_matrix=D_matrix),
               goal_function(c_perm, 100, U2, D_matrix=D_matrix))

  # Those values were calculated by hand:
  expect_equal(exp(goal_function(c_perm, 100, U1*2, D_matrix=D_matrix*2)),
               6^(-103/2) * gamma(103/2) * gamma(103/2) / (pi / 4))
  expect_equal(exp(goal_function(id_perm, 100, U1*2, D_matrix=D_matrix*2)),
               (23/4)^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
  expect_equal(exp(goal_function(id_perm, 100, U2*2, D_matrix=D_matrix*2)),
               6^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
})

test_that('goal_function has the desired property', {
  # Example from the paper chapter 5
  # This test is randomized.
    # It is mathematically possible the Z variables will be drawn such that
    # the test fails. See ISSUE#9 for discussion.

  p <- 10
  n <- 20

  mu <- numeric(p)
  sigma_matrix <- matrix(numeric(p*p), nrow=p)
  for(i in 1:p){
    for(j in 1:p){
      sigma_matrix[i,j] <- 1 - min(abs(i-j), p-abs(i-j)) / p
    }
    sigma_matrix[i,i] <- 1 + 1/p
  }

  Z <- MASS::mvrnorm(n, mu = mu, Sigma = sigma_matrix)
  U <- t(Z) %*% Z

  actual_permutation <- permutations::as.cycle(permutations::as.word(c(2:p, 1)))
  actual_permutation_function_value <- goal_function(actual_permutation, n, U)
  another_permutation_function_value <- goal_function(id_perm, n, U)

  # We want the goal function to have a bigger value for the real permutation than for the another
  expect_gt(actual_permutation_function_value,
            another_permutation_function_value)
})

gips_example_perm <- gips_perm(example_perm, 6)
D_matrix <- matrix(c(10, 1, 1, 2, 2, 3,
      1, 10, 1, 2, 2, 3,
      1, 1, 10, 2, 2, 3,
      2, 2, 2, 12, 4, 5,
      2, 2, 2, 4, 12, 5,
      3, 3, 3, 5, 5, 14), byrow = TRUE, ncol=6) * 2
delta <- 3
n <- 1
U_matrix <- diag(6) * 2
structure_constants <- get_structure_constants(gips_example_perm)
expected_phi_part <- log(2206^(-3) * 10^(-3/2) * 9^(-2)) - log(1680^(-5/2) * 8^(-3/2))

test_that('calculate phi_part returns proper values', {
    expect_equal(calculate_phi_part(gips_example_perm, n, U_matrix, delta, D_matrix,
                                    structure_constants),
                 expected_phi_part)
})

test_that('calculate_block_determinants returns proper values', {
    expect_equal(calculate_determinants_of_block_matrices(block_diagonal_matrix,
                                                          c(2,3,5)),
                 c(3,2,1))
})


