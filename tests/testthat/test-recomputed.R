# This tests that combineRecomputedResults works as expected.
# library(testthat); library(SingleR); source("test-recomputed.R")

library(testthat); library(SingleR); 
set.seed(10000)

# Making up data (using an uneven distribution to avoid symmetry masking problems).
ref <- .mockRefData(nreps=8)
ref1 <- ref[,1:4%%4==0]
ref1 <- ref1[,sample(ncol(ref1))]

ref2 <- ref[,1:4%%4!=0]
ref2 <- ref2[,sample(ncol(ref2))]

ref2$label <- tolower(ref2$label)
test <- .mockTestData(ref)

# Computing scores.
test <- scuttle::logNormCounts(test)

ref1 <- scuttle::logNormCounts(ref1)
train1 <- trainSingleR(ref1, labels=ref1$label, test.genes=rownames(test))
pred1 <- classifySingleR(test, train1)

ref2 <- scuttle::logNormCounts(ref2)
train2 <- trainSingleR(ref2, labels=ref2$label, test.genes=rownames(test))
pred2 <- classifySingleR(test, train2)

test_that("combineRecomputedResults works as expected (light check)", {
    combined <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=test,
        trained=list(train1, train2))

    # Checking the sanity of the output.
    obs <- apply(combined$scores, 1, FUN=function(x) colnames(combined$scores)[!is.na(x)])
    ref <- rbind(pred1$labels, pred2$labels) 
    expect_identical(obs, ref)

    expect_true(all(combined$labels == pred1$labels | combined$labels==pred2$labels))
    expect_true(all(combined$first.labels == pred1$first.labels | combined$first.labels==pred2$first.labels))

    expect_true(all(
        combined$pruned.labels==pred1$pruned.labels | 
        combined$pruned.labels==pred2$pruned.labels |
        is.na(combined$pruned.labels)==is.na(pred1$pruned.labels) | 
        is.na(combined$pruned.labels)==is.na(pred2$pruned.labels)
    ))

    top <- apply(combined$scores, 1, FUN=function(x) colnames(combined$scores)[which.max(x)])
    expect_identical(top, combined$labels)
})

test_that("combineRecomputedResults matrix fragmentation works as expected", {
    combined1 <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=test,
        trained=list(train1, train2))

    # Testing that it works upon parallelization.
    combined1x <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=test,
        trained=list(train1, train2),
        BPPARAM=BiocParallel::MulticoreParam(2))
    expect_equal(combined1, combined1x)

    # Testing that it works for DA's, as well as when the DA
    # has memory limits that need to be respected.
    library(DelayedArray)
    DA <- DelayedArray(assay(test))
    combined2a <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=DA,
        trained=list(train1, train2))
    expect_equal(combined1, combined2a)

    old <- getAutoBlockSize()
    setAutoBlockSize(nrow(DA)*8L)
    combined2b <- combineRecomputedResults(
        results=list(pred1, pred2),
        test=DA,
        trained=list(train1, train2))
    expect_equal(combined1, combined2b)

    setAutoBlockSize(old)
})

test_that("combineRecomputedResults handles mismatches to rows and cells", {
    expect_error(combineRecomputedResults(
        results=list(pred1[1:10,], pred2), 
        test=test,
        trained=list(train1, train2)), "not identical")

    # Responds to differences in the cell names.
    colnames(test) <- seq_len(ncol(test))
    expect_error(combineRecomputedResults(
        results=list(pred1, pred2),
        test=test[,1],
        trained=list(train1, train2)), "not identical")
    colnames(test) <- NULL

    # Responds to mismatches in the genes.
    s <- sample(nrow(test))
    expect_error(combineRecomputedResults(
        results=list(pred1, pred2),
        test=test[s,],
        trained=list(train1, train2)), "test.genes")
})

test_that("combineRecomputedResults emits warnings when missing genes are present", {
    half <- nrow(ref) / 2

    # Spiking in some missing genes.
    train1b <- trainSingleR(ref1[seq_len(half),,drop=FALSE], labels=ref1$label, test.genes=rownames(test))
    pred1b <- classifySingleR(test, train1b)

    train2b <- trainSingleR(ref2[half + seq_len(half),,drop=FALSE], labels=ref2$label, test.genes=rownames(test))
    pred2b <- classifySingleR(test, train2b)

    expect_warning(out <- combineRecomputedResults(
        results=list(pred1b, pred2b), 
        test=test,
        trained=list(train1b, train2b)), "available in each reference")
})

test_that("combineRecomputedResults works with intersections", {
    tkeep <- sample(rownames(test), 500)
    rkeep <- sample(rownames(ref), 500)

    subtest <- test[tkeep,]
    train1b <- trainSingleR(ref1[rownames(ref1) %in% rkeep,,drop=FALSE], labels=ref1$label, test.genes=tkeep)
    pred1b <- classifySingleR(subtest, train1b)

    train2b <- trainSingleR(ref2[rownames(ref2) %in% rkeep,,drop=FALSE], labels=ref2$label, test.genes=tkeep)
    pred2b <- classifySingleR(subtest, train2b)

    out <- combineRecomputedResults(
        results=list(pred1b, pred2b), 
        test=subtest,
        trained=list(train1b, train2b)
    )

    # Comparing to some explicit subsets to the intersection.
    common <- intersect(tkeep, rkeep)
    subtest <- test[rownames(test) %in% common,]
    train1c <- trainSingleR(ref1[rownames(ref1) %in% common,], labels=ref1$label)
    pred1c <- classifySingleR(subtest, train1c)

    train2c <- trainSingleR(ref2[rownames(ref2) %in% common,], labels=ref2$label)
    pred2c <- classifySingleR(subtest, train2c)

    ref.out <- combineRecomputedResults(
        results=list(pred1c, pred2c), 
        test=subtest,
        trained=list(train1c, train2c)
    )

    expect_equal(out$scores, ref.out$scores)
    expect_identical(out$labels, ref.out$labels)
})

test_that("combineRecomputedResults is invariant to ordering", {
    ref3 <- .mockRefData(nreps=8)
    ref3 <- scuttle::logNormCounts(ref3)
    ref3$label <- paste0(ref3$label, "X")
    train3 <- trainSingleR(ref3, labels=ref3$label)
    pred3 <- classifySingleR(test, train3)

    combined <- combineRecomputedResults(
        results=list(pred1, pred2, pred3),
        test=test,
        trained=list(train1, train2, train3))

    flipped <- combineRecomputedResults(
        results=rev(list(pred1, pred2, pred3)),
        test=test,
        trained=rev(list(train1, train2, train3)))

    expect_identical(flipped$orig.results[,3], combined$orig.results[,1])
    expect_identical(flipped$orig.results[,2], combined$orig.results[,2])
    expect_identical(flipped$orig.results[,1], combined$orig.results[,3])
    expect_identical(flipped$labels, combined$labels)
})
