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
