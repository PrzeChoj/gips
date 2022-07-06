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
S1 <- U1/100
U2 <- matrix(c(1.5,0.5,0.5,1.5), nrow=2, byrow = TRUE)
S2 <- U2/100
D_matrix <- diag(nrow = 2)


test_that('log_likelihood_of_perm returns proper values', {
  # The value of log_likelihood_of_perm on matrix should the same as the projection of the matrix
  # and U2 == pi_c(U1)
  expect_equal(log_likelihood_of_perm(c_perm, 100, S1, D_matrix=D_matrix),
               log_likelihood_of_perm(c_perm, 100, S2, D_matrix=D_matrix))

  # Those values were calculated by hand:
  expect_equal(exp(log_likelihood_of_perm(c_perm, 100, S1*2, D_matrix=D_matrix*2)),
               6^(-103/2) * gamma(103/2) * gamma(103/2) / (pi / 4))
  expect_equal(exp(log_likelihood_of_perm(id_perm, 100, S1*2, D_matrix=D_matrix*2)),
               (23/4)^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
  expect_equal(exp(log_likelihood_of_perm(id_perm, 100, S2*2, D_matrix=D_matrix*2)),
               6^(-52) * gamma(52) * gamma(51.5) * sqrt(2*pi) / (pi / sqrt(2)))
})

test_that('log_likelihood_of_perm has the desired property', {
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
  S <- U/n

  actual_permutation <- permutations::as.cycle(permutations::as.word(c(2:p, 1)))
  actual_permutation_function_value <- log_likelihood_of_perm(actual_permutation, n, S)
  another_permutation_function_value <- log_likelihood_of_perm(id_perm, n, S)

  # We want the goal function to have a bigger value for the real permutation than for the another
  expect_gt(actual_permutation_function_value,
            another_permutation_function_value)
})

test_that('calculate phi_part returns proper values', {
    skip("TODO")

    expect_true(TRUE)
})

test_that('calculate_block_determinants returns proper values', {
    expect_equal(calculate_determinants_of_block_matrices(block_diagonal_matrix,
                                                          c(2,3,5)),
                 c(3,2,1))
})


