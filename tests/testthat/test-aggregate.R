# This tests the aggregator function.
# library(testthat); library(SingleR); source("setup.R"); source("test-aggregate.R")

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

test_that("aggregateReference selects HVGs correctly", {
    mat <- assay(sce)
    keep <- sample(nrow(mat), 50)
    mat[keep,] <-  mat[keep,] * 100

    labels <- sample(LETTERS, ncol(sce), replace=TRUE)
    set.seed(10)
    aggr <- aggregateReference(mat, labels, ntop=50)

    set.seed(10)
    ref <- aggregateReference(mat[keep,], labels, ntop=Inf)

    expect_identical(aggr[keep,], ref)
})

test_that("aggregateReference seed setter behaves correctly", {
    # Seed is different for each set of labels. 
    set.seed(10)
    aggr <- aggregateReference(BiocGenerics::cbind(sce, sce), rep(1:2, each=ncol(sce)))

    N <- ncol(aggr)/2
    expect_false(identical(unname(assay(aggr)[,1:N]), unname(assay(aggr)[,N+1:N])))

    # Seed is different for different runs.
    labels <- sample(LETTERS, ncol(sce), replace=TRUE)

    set.seed(10)
    X1 <- aggregateReference(sce, labels) # double usage is intentional.
    X2 <- aggregateReference(sce, labels)
    expect_false(identical(X1, X2))

    # You get the same results with the same seed, and different results with different seeds.
    labels <- sample(LETTERS, ncol(sce), replace=TRUE)

    set.seed(20)
    ref <- aggregateReference(sce, labels)

    set.seed(20)
    out <- aggregateReference(sce, labels)
    expect_identical(ref, out)

    set.seed(30)
    out <- aggregateReference(sce, labels)
    expect_false(identical(ref, out))

    # You get same results regardless of parallelization. 
    # Note that SnowParam requires us to disable the failsafes,
    # otherwise setAutoBPPARAM() doesn't work properly.
    setAutoBPPARAM(SerialParam())

    set.seed(20)
    out <- aggregateReference(sce, labels, BPPARAM=BiocParallel::SnowParam(2))
    expect_identical(ref, out)

    # The seed is unset properly for downstream applications.
    set.seed(10)
    aggregateReference(sce, labels)
    test1 <- runif(10)

    set.seed(10)
    aggregateReference(sce, labels, BPPARAM=BiocParallel::SnowParam(2))
    test2 <- runif(10)
    expect_identical(test1, test2)

    setAutoBPPARAM(FAIL)
})
