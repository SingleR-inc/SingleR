# This tests the aggregator function.
# library(testthat); library(SingleR); source("test-aggregate.R")

library(scuttle)
sce <- mockSCE()
sce <- logNormCounts(sce)

test_that("aggregateReference works as expected for full aggregation", {
    labels <- sample(LETTERS, ncol(sce), replace=TRUE)

    aggr <- aggregateReference(sce, labels, power=0)
    out <- assay(aggr)
    colnames(out) <- sub("\\..*", "", colnames(out))

    ref <- colsum(logcounts(sce), labels)
    ref <- t(t(ref)/as.vector(table(labels)))
    m <- match(colnames(out), colnames(ref))

    expect_equivalent(out, ref[,m])
    expect_identical(aggr$label, colnames(ref)[m])
})

test_that("aggregateReference works as expected for partial aggregation", {
    labels <- sample(LETTERS[1:3], ncol(sce), replace=TRUE)

    aggr <- aggregateReference(sce, labels, power=0.5)
    expect_equal(ncol(aggr), sum(floor(sqrt(table(labels)))))

    expect_identical(aggr$label, sub("\\..*", "", colnames(aggr)))

    # Contrived example with three elements per label to check the averaging is done correctly.
    labels <- head(rep(seq_len(ceiling(ncol(sce)/3)), each=3), ncol(sce))

    aggr1 <- aggregateReference(sce, labels, power=0)
    aggr2 <- aggregateReference(sce, labels, power=0.5)
    expect_equal(aggr1, aggr2)
})

test_that("aggregateReference works as expected for no aggregation", {
    labels <- sample(LETTERS, ncol(sce), replace=TRUE)
    aggr <- aggregateReference(sce, labels, power=1)
    expect_equivalent(assay(aggr)[,order(aggr$label)], logcounts(sce)[,order(labels)])
})

test_that("aggregateReference skips PCA when the requested rank is too high", {
    labels <- sample(LETTERS, ncol(sce), replace=TRUE)
    ref <- aggregateReference(sce, labels)

    expect_error(out <- aggregateReference(sce, labels, rank=1000, BSPARAM="WHEEE"), NA) # should not even require the BSPARAM.
    expect_identical(ncol(ref), ncol(out))
})

test_that("aggregateReference works as expected for empty inputs", {
    labels <- sample(LETTERS, ncol(sce), replace=TRUE)
    aggr <- aggregateReference(sce[,0], labels[0])
    expect_identical(nrow(aggr), nrow(sce))
    expect_identical(ncol(aggr), 0L)
})
