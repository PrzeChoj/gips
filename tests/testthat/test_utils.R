perm_nofixed <- to_perm(c(2,1,4,5,3))
perm_somefixed <- to_perm(c(2,1,3,5,4))
perm_firstfixed <- to_perm(c(1,3,4,5,2))
perm_lastfixed <- to_perm(c(2,1,4,3,5))
perm_allfixed <- to_perm(1:5)

test_that('get_subcycles works for no fixed elements',{
    expect_equal(get_subcycles(perm_nofixed, 5),
                 list(1:2, 3:5))
})

test_that('get_cycle_rep_lengths works for some fixed elements',{
    expect_equal(get_subcycles(perm_somefixed, 5),
                 list(1:2, 3, 4:5))
})

test_that('get_cycle_rep_lengths works for first fixed element',{
    expect_equal(get_subcycles(perm_firstfixed, 5),
                 list(1, 2:5))
})

test_that('get_cycle_rep_lengths works for last fixed element',{
    expect_equal(get_subcycles(perm_lastfixed, 5),
                 list(1:2, 3:4, 5))
})

test_that('get_cycle_rep_lengths works for identity',{
    expect_equal(get_subcycles(perm_allfixed, 5),
                 list(1,2,3,4,5))
})
