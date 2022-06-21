# This tests the rebuildIndex function.
# library(testthat); library(SingleR); source("setup.R"); source("test-rebuild.R")

trained <- trainSingleR(training, training$label)

test_that("rebuildIndex works correctly", {
    ref <- classifySingleR(test, trained)

    tmp <- tempfile(fileext=".rds")
    saveRDS(trained, file=tmp)
    reloaded <- readRDS(tmp)
    expect_false(SingleR:::is_valid_built(reloaded$built))

    rebuilt <- rebuildIndex(reloaded)
    expect_true(SingleR:::is_valid_built(rebuilt$built))
    again <- classifySingleR(test, rebuilt)
    expect_equal(again, ref)
})

test_that("rebuildIndex is a no-op if required", {
    rebuilt <- rebuildIndex(trained)
    expect_identical(rebuilt, trained)

    tmp <- tempfile(fileext=".rds")
    saveRDS(trained, file=tmp)
    reloaded <- readRDS(tmp)
    expect_false(identical(reloaded, trained))
})

