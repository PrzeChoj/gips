suppressMessages(library(permutations))

p1 <- as.word(c(2,1,3,6,4,5))
expected_p <- as.word(c(2,1,3,5,6,4))

test_that('get_representative example',{
    expect_equal(cycle2word(get_representative(p1)),
                 expected_p)
})


