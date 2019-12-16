# Test that all datasets can be loaded with their arguments.
# library(testthat); library(SingleR); source("test-datasets.R")

test_that("Ensembl coercion works for all datasets", {
    expect_error(X <- BlueprintEncodeData(ensembl=TRUE), NA)

    expect_error(X <- DatabaseImmuneCellExpressionData(ensembl=TRUE), NA)

    expect_error(X <- HumanPrimaryCellAtlasData(ensembl=TRUE), NA)

    expect_error(X <- ImmGenData(ensembl=TRUE), NA)

    expect_error(X <- MonacoImmuneData(ensembl=TRUE), NA)

    expect_error(X <- MouseRNAseqData(ensembl=TRUE), NA)

    expect_error(X <- NovershternHematopoieticData(ensembl=TRUE), NA)
})
