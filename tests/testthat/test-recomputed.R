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
train1 <- trainSingleR(ref1, labels=ref1$label)
pred1 <- classifySingleR(test, train1)

ref2 <- scuttle::logNormCounts(ref2)
train2 <- trainSingleR(ref2, labels=ref2$label)
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
        BPPARAM=BiocParallel::MulticoreParam(3))
    expect_equal(combined1, combined1x)

#    # Testing that it works for DA's, as well as when the DA
#    # has memory limits that need to be respected.
#    library(DelayedArray)
#    DA <- DelayedArray(assay(test))
#    combined2a <- combineRecomputedResults(
#        results=list(pred1, pred2), 
#        test=DA,
#        trained=list(train1, train2))
#    expect_equal(combined1, combined2a)
#
#    old <- getAutoBlockSize()
#    setAutoBlockSize(nrow(DA)*8L)
#    combined2b <- combineRecomputedResults(
#        results=list(pred1, pred2),
#        test=DA,
#        trained=list(train1, train2))
#    expect_equal(combined1, combined2b)
#
#    setAutoBlockSize(old)
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

    # Correctly reorders the gene universes.
    ref <- combineRecomputedResults(
        results=list(pred1, pred2),
        test=test,
        trained=list(train1, train2))

    s <- sample(nrow(test))
    out <- combineRecomputedResults(
        results=list(pred1, pred2),
        test=test[s,],
        trained=list(train1, train2))
    expect_equal(ref, out)
})

test_that("combineRecomputedResults emits warnings when missing genes are present", {
    # Spiking in some missing genes.
    ref1b <- ref1[c(1, seq_len(nrow(ref1))),]
    rownames(ref1b)[1] <- "BLAH"
    markers1 <- train1$markers$full
    markers1$A$B <- c(markers1$A$B, "BLAH")
    train1b <- trainSingleR(ref1b, labels=ref1$label, genes=markers1)

    ref2b <- ref2[c(1, seq_len(nrow(ref2))),]
    rownames(ref2b)[1] <- "WHEE"
    markers2 <- train2$markers$full
    markers2$A$B <- c(markers2$a$b, "WHEE")
    train2b <- trainSingleR(ref2b, labels=ref2$label, genes=markers2)

    test2 <- test[c(1,seq_len(nrow(test)),1),]
    rownames(test2)[1] <- "WHEE"
    rownames(test2)[length(rownames(test2))] <- "BLAH"

    expect_warning(out <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=test2,
        trained=list(train1b, train2b)), "differ in the universe")
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
