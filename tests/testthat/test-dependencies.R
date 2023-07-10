test_that("permutations did not change the order in allperms", {
  all_perms_4 <- permutations::as.cycle(permutations::allperms(4))
  expect_equal(as.character(all_perms_4[5]), "(2,4,3)")
  expect_equal(as.character(all_perms_4[21]), "(1,4,3)")
})
