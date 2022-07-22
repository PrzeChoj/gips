perm_nofixed <- to_perm(c(2,1,4,5,3))
gips_perm_nofixed <- structure(list(1:2, 3:5), size=5, class='gips_perm')
perm_somefixed <- to_perm(c(2,1,3,5,4))
gips_perm_somefixed <- structure(list(1:2, 3, 4:5), size=5, class='gips_perm')
perm_firstfixed <- to_perm(c(1,3,4,5,2))
gips_perm_firstfixed <- structure(list(1, 2:5), size=5, class='gips_perm')
perm_lastfixed <- to_perm(c(2,1,4,3,5))
gips_perm_lastfixed <- structure(list(1:2, 3:4, 5), size=5, class='gips_perm')
perm_allfixed <- to_perm(1:5)
gips_perm_allfixed <- structure(list(1,2,3,4,5), size=5, class='gips_perm')

test_that('gips_perm works for no fixed elements',{
    expect_equal(gips_perm(perm_nofixed, 5),
                 gips_perm_nofixed)
})

test_that('gips_perm works for some fixed elements',{
    expect_equal(gips_perm(perm_somefixed, 5),
                 gips_perm_somefixed)
})

test_that('gips_perm works for first fixed element',{
    expect_equal(gips_perm(perm_firstfixed, 5),
                 gips_perm_firstfixed)
})

test_that('gips_perm works for last fixed element',{
    expect_equal(gips_perm(perm_lastfixed, 5),
                 gips_perm_lastfixed)
})

test_that('gips_perm works for identity',{
    expect_equal(gips_perm(perm_allfixed, 5),
                 gips_perm_allfixed)
})

test_that('constructor works for empty permutations',{
    expect_equal(new_gips_perm(list(), 0),
                 structure(list(), size=0, class='gips_perm'))
})

test_that('gips_perm works for empty permutation',{
    expect_equal(gips_perm(permutations::nullword, 0),
                 structure(list(), size=0, class='gips_perm'))
})

test_that('gips_perm warns when multiple permutations passed',{
    expect_warning(out <- gips_perm(c('(1,2,3)', '(1,3,2)'), 3),
                   'multiple permutations')
    expect_true(identical(gips_perm('(1,2,3)', 3),
                          out))
})

test_that('gips_perm coerces to permutations::cycle',{
    expect_equal(permutations::as.cycle(gips_perm_nofixed),
                 perm_nofixed)
    expect_equal(permutations::as.cycle(gips_perm_somefixed),
                 perm_somefixed)
    expect_equal(permutations::as.cycle(gips_perm_firstfixed),
                 perm_firstfixed)
    expect_equal(permutations::as.cycle(gips_perm_lastfixed),
                 perm_lastfixed)
    expect_equal(permutations::as.cycle(gips_perm_allfixed),
                 perm_allfixed)
})

gips_example_perm <- gips_perm(example_perm, 6)
gips_example_perm_copy <- gips_perm(example_perm, 6)
gips_example_perm_different_size <- gips_perm(example_perm, 5)
gips_different_perm <- gips_perm(example_perm2, 6)

test_that('identical() works with gips_perms', {
    expect_true(identical(gips_example_perm, gips_example_perm_copy))
    expect_false(identical(gips_example_perm, gips_example_perm_different_size))
    expect_false(identical(gips_example_perm, gips_different_perm))
})


# TODO more exact tests

transpositions <- expand.grid(1:6, 1:6)
transpositions <- as.matrix(transpositions[transpositions$Var1 < transpositions$Var2,])
for(i in 1:nrow(transpositions)){
    transposition <- transpositions[i,]
    test_that(paste0('compose_with_transposition works for (',
                     transposition[1], ', ', transposition[2], ')'),
              {
                  tr_perm <- permutations::as.cycle(transposition)
                  expected <- gips_perm(example_perm * tr_perm, 6)
                  actual <- compose_with_transposition(gips_example_perm,
                                                       transposition)
                  expect_equal(actual, expected)
              })
}

