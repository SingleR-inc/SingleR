# library(testthat); library(SingleR); source("test-utils.R")

test_that("intersection works as expected for edge cases", {
    out <- SingleR:::.create_intersection(c("A", "B", "C", "a"), LETTERS)
    expect_identical(out$test, 1:3)
    expect_identical(out$reference, 1:3)

    out <- SingleR:::.create_intersection(c("a", "Z", "b", "G", "c", "B", "d"), LETTERS)
    expect_identical(out$test, c(2L, 4L, 6L))
    expect_identical(out$reference, c(26L, 7L, 2L))

    out <- SingleR:::.create_intersection(c("y", "y", "x", "z", NA, "z", "x"), letters)
    expect_identical(out$test, c(1L, 3L, 4L))
    expect_identical(out$reference, c(25L, 24L, 26L))
})
