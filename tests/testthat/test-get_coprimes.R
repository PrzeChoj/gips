test_that("get_coprimes examples", {
  expect_equal(get_coprimes(10), c(1, 3, 7, 9))
  expect_equal(length(get_coprimes(100)), numbers::eulersPhi(100))
})
