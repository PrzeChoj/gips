to_perm <- function(v) permutations::as.cycle(permutations::as.word(v))
perm_nofixed <- to_perm(c(2,1,5,3,4))
perm_somefixed <- to_perm(c(2,1,3,5,4))
perm_firstfixed <- to_perm(c(1,5,2,3,4))
perm_lastfixed <- to_perm(c(2,1,4,3,5))
perm_allfixed <- to_perm(1:5)
example_perm <- to_perm(c(3,1,2,5,4,6))

test_that('get_structure_constants for example from paper',{
    expect_equal(get_structure_constants(example_perm, 6),
                 list('r'=c(3,1,1),
                      'd'=c(1,2,1),
                      'k'=c(1,2,1),
                      'L'=3))
})

test_that('get_cycle_rep_lengths works for no fixed elements',{
    expect_equal(get_cycle_representatives_and_lengths(perm_nofixed, 5),
                     list('representatives' = c(1,3),
                          'cycle_lengths' = c(2,3)))
})

test_that('get_cycle_rep_lengths works for some fixed elements',{
    expect_equal(get_cycle_representatives_and_lengths(perm_somefixed, 5),
                     list('representatives' = c(1,3,4),
                          'cycle_lengths' = c(2,1,2)))
})

test_that('get_cycle_rep_lengths works for first fixed element',{
    expect_equal(get_cycle_representatives_and_lengths(perm_firstfixed, 5),
                     list('representatives' = c(1,2),
                          'cycle_lengths' = c(1,4)))
})

test_that('get_cycle_rep_lengths works for last fixed element',{
    expect_equal(get_cycle_representatives_and_lengths(perm_lastfixed, 5),
                 list('representatives' = c(1,3,5),
                      'cycle_lengths' = c(2,2,1)))
})

test_that('get_cycle_rep_lengths works for identity',{
    expect_equal(get_cycle_representatives_and_lengths(perm_allfixed, 5),
                 list('representatives' = 1:5,
                      'cycle_lengths' = rep(1,5)))
})

# Example 6 from paper
test_that('calculate_r works for example from paper', {
    expect_equal(calculate_r(c(3,2,1), 6),
                 c(3,0,1,1))
})

test_that('calculate_r works for identity', {
    expect_equal(calculate_r(c(1,1,1,1), 1),
                 4)
})

test_that('calculate_d works for even perm_order', {
    expect_equal(calculate_d(6),c(1,2,2,1))
    expect_equal(calculate_d(4), c(1,2,1))
    expect_equal(calculate_d(2), c(1,1))
})

test_that('calculate_d works for odd perm_order', {
    expect_equal(calculate_d(5), c(1,2,2))
    expect_equal(calculate_d(3), c(1,2))
})

test_that('calculate_d works for identity', {
    expect_equal(calculate_d(1), 1)
})
