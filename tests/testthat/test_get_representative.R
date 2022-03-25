p1 <- permutations::as.word(c(2,1,3,6,4,5))
expected_p <- permutations::as.word(c(2,1,3,5,6,4))

test_that('get_representative examples',{
    expect_equal(permutations::cycle2word(get_representative(p1)),
                 expected_p)
  
    expect_equal(get_representative(permutations::as.word(numeric(0))),
                 permutations::as.cycle(permutations::as.word(numeric(0))))
})


