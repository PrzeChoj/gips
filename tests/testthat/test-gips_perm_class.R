test_that("gips_perm works for no fixed elements", {
  perm_nofixed <- to_perm(c(2, 1, 4, 5, 3))
  gips_perm_nofixed <- structure(list(1:2, 3:5), size = 5, class = "gips_perm")

  expect_equal(
    gips_perm(perm_nofixed, 5),
    gips_perm_nofixed
  )
})

test_that("gips_perm works for some fixed elements", {
  perm_somefixed <- to_perm(c(2, 1, 3, 5, 4))
  gips_perm_somefixed <- structure(list(1:2, 3, 4:5), size = 5, class = "gips_perm")

  expect_equal(
    gips_perm(perm_somefixed, 5),
    gips_perm_somefixed
  )
})

test_that("gips_perm works for first fixed element", {
  perm_firstfixed <- to_perm(c(1, 3, 4, 5, 2))
  gips_perm_firstfixed <- structure(list(1, 2:5), size = 5, class = "gips_perm")

  expect_equal(
    gips_perm(perm_firstfixed, 5),
    gips_perm_firstfixed
  )
})

test_that("gips_perm works for last fixed element", {
  perm_lastfixed <- to_perm(c(2, 1, 4, 3, 5))
  gips_perm_lastfixed <- structure(list(1:2, 3:4, 5), size = 5, class = "gips_perm")

  expect_equal(
    gips_perm(perm_lastfixed, 5),
    gips_perm_lastfixed
  )
})

test_that("gips_perm works for identity", {
  perm_allfixed <- to_perm(1:5)
  gips_perm_allfixed <- structure(list(1, 2, 3, 4, 5), size = 5, class = "gips_perm")

  expect_equal(
    gips_perm(perm_allfixed, 5),
    gips_perm_allfixed
  )
})

test_that("gips_perm gives the same output for different type of input", {
  gperm1 <- gips_perm("(12)(45)", 5)
  gperm2 <- gips_perm("(1,2)(4,5)", 5)
  gperm3 <- gips_perm(as.matrix(c(2, 1, 3, 5, 4)), 5)
  gperm4 <- gips_perm(t(as.matrix(c(2, 1, 3, 5, 4))), 5) # both way works
  gperm5 <- gips_perm(list(list(c(2, 1), c(4, 5))), 5)
  gperm6 <- gips_perm(permutations::as.word(c(2, 1, 3, 5, 4)), 5)
  gperm7 <- gips_perm(permutations::as.cycle("(1,2)(4,5)"), 5)
  gperm8 <- gips_perm(gips(diag(5), 14, perm = "(12)(45)"), 5) # gips class

  # expect_equal and not expect_identical, because those are int / double:
  expect_equal(gperm1, gperm2)
  expect_equal(gperm1, gperm3)
  expect_equal(gperm1, gperm4)
  expect_equal(gperm1, gperm5)
  expect_equal(gperm1, gperm6)
  expect_equal(gperm1, gperm7)
  expect_equal(gperm1, gperm8)
})

test_that("gips_perm shows an error for wrong input", {
  g <- gips(diag(5), 14, perm = "(12)(45)")
  expect_error(gips_perm(g, 4), regexp = ", which in general is OK, but You also provided size = 4, which is different from attr")
})

test_that("constructor works for empty permutations", {
  expect_equal(
    new_gips_perm(list(), 0),
    structure(list(), size = 0, class = "gips_perm")
  )
})

test_that("gips_perm works for empty permutation", {
  expect_equal(
    gips_perm(permutations::nullword, 0),
    structure(list(), size = 0, class = "gips_perm")
  )
})

test_that("gips_perm warns when multiple permutations passed", {
  expect_warning(
    out <- gips_perm(c("(1,2,3)", "(1,3,2)"), 3),
    "multiple permutations"
  )
  expect_warning(
    out <- gips_perm(c("(1,2,3)", "(1,3,2)"), 3),
    "2 permutations"
  )
  expect_true(identical(
    gips_perm("(1,2,3)", 3),
    out
  ))
})

test_that("gips_perm coerces to permutations::cycle", {
  perm_nofixed <- to_perm(c(2, 1, 4, 5, 3))
  gips_perm_nofixed <- structure(list(1:2, 3:5), size = 5, class = "gips_perm")

  expect_equal(
    permutations::as.cycle(gips_perm_nofixed),
    perm_nofixed
  )


  perm_somefixed <- to_perm(c(2, 1, 3, 5, 4))
  gips_perm_somefixed <- structure(list(1:2, 3, 4:5), size = 5, class = "gips_perm")

  expect_equal(
    permutations::as.cycle(gips_perm_somefixed),
    perm_somefixed
  )


  perm_firstfixed <- to_perm(c(1, 3, 4, 5, 2))
  gips_perm_firstfixed <- structure(list(1, 2:5), size = 5, class = "gips_perm")

  expect_equal(
    permutations::as.cycle(gips_perm_firstfixed),
    perm_firstfixed
  )


  perm_lastfixed <- to_perm(c(2, 1, 4, 3, 5))
  gips_perm_lastfixed <- structure(list(1:2, 3:4, 5), size = 5, class = "gips_perm")

  expect_equal(
    permutations::as.cycle(gips_perm_lastfixed),
    perm_lastfixed
  )


  perm_allfixed <- to_perm(1:5)
  gips_perm_allfixed <- structure(list(1, 2, 3, 4, 5), size = 5, class = "gips_perm")

  expect_equal(
    permutations::as.cycle(gips_perm_allfixed),
    perm_allfixed
  )
})

test_that("gips_perm handles bad arguments", {
  perm_somefixed <- to_perm(c(2, 1, 3, 5, 4))

  expect_error(gips_perm(""), "You did not provide the `size` arg")
  expect_warning(
    gips_perm("", c(4, 3)),
    "Passing multiple sizes to \\`gips_perm\\(\\)\\` is not supported. Taking only the first one"
  )
  expect_error(
    gips_perm("", 3.2),
    "You provided `size == 3.2`"
  )
  expect_error(
    gips_perm(perm_somefixed, 3),
    "size"
  )
  expect_warning(
    gips_perm(c("", "(1,2)"), 3),
    "Passing multiple permutations to \\`gips_perm\\(\\)\\` is not supported. Taking only the first one"
  )
})


test_that("validate_gips works properly", {
  not_gips_perm <- "A"
  not_a_list <- structure("A", class = "gips_perm")
  without_size <- structure(list(), class = "gips_perm")
  wrong_size_type <- structure(list(), class = "gips_perm", size = "A")
  wrong_element <- structure(list("A"), class = "gips_perm", size = 3)
  duplicate_element <- structure(list(1:2, 1), class = "gips_perm", size = 3)
  not_sorted_internally <- structure(list(1:2, c(6, 5, 4)),
    class = "gips_perm", size = 6
  )
  cycles_not_sorted <- structure(list(3:5, 1:2),
    class = "gips_perm", size = 6
  )
  small_size <- structure(list(1:2, 3:5),
    class = "gips_perm",
    size = 3
  )


  gips_perm_somefixed <- structure(list(1:2, 3, 4:5), size = 5, class = "gips_perm")

  expect_silent(validate_gips_perm(gips_perm_somefixed))
  expect_error(validate_gips_perm(not_gips_perm), "class")
  expect_error(validate_gips_perm(not_a_list), "list")
  expect_error(validate_gips_perm(without_size), "size")
  expect_error(validate_gips_perm(wrong_size_type), "integer")
  expect_error(validate_gips_perm(wrong_element), "integer")
  expect_error(validate_gips_perm(duplicate_element), "repeat")
  expect_error(validate_gips_perm(not_sorted_internally), "First element .* minimum")
  expect_error(validate_gips_perm(cycles_not_sorted), "order determined by their first elements")
  expect_error(validate_gips_perm(small_size), "size")
})

test_that("object of gips_perm class can be printed", {
  expect_output(print(gips_perm("(1,2)(5,4)", 7)),
    regexp = "\\(12\\)\\(45\\)"
  )
})


test_that("identical() works with gips_perms", {
  gips_example_perm <- gips_perm(example_perm, 6)
  gips_example_perm_copy <- gips_perm(example_perm, 6)
  gips_example_perm_different_size <- gips_perm(example_perm, 5)
  gips_different_perm <- gips_perm(example_perm2, 6)

  expect_true(identical(gips_example_perm, gips_example_perm_copy))
  expect_false(identical(gips_example_perm, gips_example_perm_different_size))
  expect_false(identical(gips_example_perm, gips_different_perm))
})


# 36 tests, one for every (i, j) transposition
transpositions <- expand.grid(1:6, 1:6)
transpositions <- as.matrix(transpositions[transpositions[["Var1"]] < transpositions[["Var2"]], ])
for (i in 1:nrow(transpositions)) {
  transposition <- transpositions[i, ]
  test_that(paste0(
    "compose_with_transposition works for (",
    transposition[1], ", ", transposition[2], ")"
  ), {
    gips_example_perm <- gips_perm(example_perm, 6)
    tr_perm <- permutations::as.cycle(transposition)
    expected <- gips_perm(example_perm * tr_perm, 6)
    actual <- compose_with_transposition(
      gips_example_perm,
      transposition
    )
    expect_equal(actual, expected)
  })
}

test_that("rearrange_cycles() works properly", {
  cycles <- list(c(2, 4, 3), c(5, 1))
  rearranged <- rearrange_cycles(cycles)

  expect_true(identical(rearranged, list(c(1, 5), c(2, 4, 3))))
  expect_true(identical(rearrange_cycles(list()), list()))

  rand_size <- sample(1000, 1)
  x <- as.list(1:rand_size)
  expect_true(identical(rearrange_cycles(x), x))
})
