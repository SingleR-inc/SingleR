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

test_that("ontology mapping works for all datasets", {
    expect_error(X <- BlueprintEncodeData(cell.ont="nonna"), NA)

    expect_error(X <- DatabaseImmuneCellExpressionData(cell.ont="nonna"), NA)

    expect_error(X <- HumanPrimaryCellAtlasData(cell.ont="nonna"), NA)

    expect_error(X <- ImmGenData(cell.ont="nonna"), NA)

    expect_error(X <- MonacoImmuneData(cell.ont="nonna"), NA)

    expect_error(X <- MouseRNAseqData(cell.ont="nonna"), NA)

    expect_error(X <- NovershternHematopoieticData(cell.ont="nonna"), NA)
})
