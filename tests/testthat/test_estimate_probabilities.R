p1 <- permutations::as.word(c(2,1,3,6,4,5))
expected_p <- permutations::as.word(c(2,1,3,5,6,4))

points <- list(to_perm(c(2,3,1)),
               to_perm(c(3,1,2)),
               to_perm(c(1,3,2)),
               to_perm(c(1,3,2)),
               to_perm(c(1,2,3)))

expected_counts <- c(2/2, 2/1, 1)/5
expected_probabilities <- expected_counts / sum(expected_counts)
names(expected_probabilities) <- c('(1,2,3)', '(2,3)', '()')

test_that('estimate_probabilities works',{
    actual_probabilities <- estimate_probabilities(points)
    expect_setequal(names(actual_probabilities),
                    names(expected_probabilities))
    expect_equal(actual_probabilities[names(expected_probabilities)],
                 expected_probabilities)
})

test_that('get_representative examples',{
    expect_equal(permutations::cycle2word(get_group_representative(p1)),
                 expected_p)

    expect_equal(get_group_representative(permutations::as.word(numeric(0))),
                 permutations::as.cycle(permutations::as.word(numeric(0))))
})


