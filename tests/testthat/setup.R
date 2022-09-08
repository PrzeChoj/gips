require("MASS") # for mvrnorm()

example_perm <- to_perm(c(2, 3, 1, 5, 4, 6))
example_perm2 <- to_perm(c(2, 3, 4, 5, 1))

matrix_invariant_by_example_perm <- matrix(c(
  12, 1, 1, 3, 3, 4,
  1, 12, 1, 3, 3, 4,
  1, 1, 12, 3, 3, 4,
  3, 3, 3, 15, 6, 7,
  3, 3, 3, 6, 15, 7,
  4, 4, 4, 7, 7, 18
), byrow = TRUE, ncol = 6)
