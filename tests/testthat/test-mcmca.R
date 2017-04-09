context('mcmca')

## Tests for [.mcmca

test_that('indexing returns base type when iteration specified', {
    x <- mcmca(1 : 2)
    expect_equal(x[2], 2)
    expect_false(class(x[2]) == 'mcmca')

    x <- mcmca(matrix(1 : 4, nrow = 2))
    expect_equal(x[2, ], c(2, 4))
    expect_false(class(x[2, ]) == 'mcmca')

    x <- mcmca(array(1 : 8, dim = c(2, 2, 2)))
    expect_equal(x[2, , ], matrix(c(2, 4, 6, 8), nrow = 2))
    expect_false(class(x[2, , ]) == 'mcmca')
})

test_that('indexing returns mcmca when iteration not specified', {
    x <- mcmca(matrix(1 : 4, nrow = 2))
    expect_equal(x[, 2], mcmca(c(3, 4)))
    expect_is(x[, 2], 'mcmca')

    x <- mcmca(array(1 : 8, dim = c(2, 2, 2)))
    expect_equal(x[, 2, ], mcmca(matrix(c(3, 4, 7, 8), nrow = 2)))
    expect_is(x[, 2, ], 'mcmca')
})
