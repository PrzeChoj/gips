# examples already tested

perm_size <- 6
mu <- numeric(perm_size)
# sigma is a matrix invariant under permutation (1,2,3,4,5,6)
sigma_matrix <- matrix(data = c(1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
                                0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
                                0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
                                0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
                                0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
                                0.8, 0.6, 0.4, 0.6, 0.8, 1.0),
                       nrow=perm_size, byrow = TRUE)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- (t(Z) %*% Z) / number_of_observations  # the theoretical mean is 0

test_that('Setting custom permutation in gips constructor works',{
  custom_perm1 <- gips_perm('(1,2)(3,4,5,6)', 6)
  g1 <- gips(S, number_of_observations, perm = custom_perm1)
  expect_identical(custom_perm1, g1[[1]])
  
  custom_perm2 <- permutations::permutation('(1,2)(3,4,5)(6)')  # notice, the `custom_perm2` does not remember the 6
  g2 <- gips(S, number_of_observations, perm = custom_perm2)
  expect_identical(gips_perm(custom_perm2, ncol(S)), g2[[1]])
})








