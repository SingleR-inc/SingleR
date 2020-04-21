# This tests that combineRecomputedResults works as expected.
# library(testthat); library(SingleR); source("test-recomputed.R")

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
test <- scater::logNormCounts(test)

ref1 <- scater::logNormCounts(ref1)
train1 <- trainSingleR(ref1, labels=ref1$label)
pred1 <- classifySingleR(test, train1)

ref2 <- scater::logNormCounts(ref2)
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

# Creating a reference function to test the recomputation.
REF <- function(scores, test, results, refs, subset, mode="de") {
    for (i in subset) {
        labs <- lapply(results, function(r) r[i,"labels"])

        if (mode=="de") {
            markers <- mapply(results, labs, FUN=function(r, l) metadata(r)$de.genes[[l]], SIMPLIFY=FALSE)
        } else {
            markers <- lapply(results, FUN=function(r) metadata(r)$common.genes)
        }
        all.markers <- unique(unlist(markers))

        keep <- mapply(refs, labs, FUN=function(r, l) r$label==l, SIMPLIFY=FALSE) 
        origins <- rep(seq_along(keep), vapply(keep, sum, 0L))
        all.refs <- mapply(refs, keep, FUN=function(R, k) R[all.markers,k,drop=FALSE], SIMPLIFY=FALSE)
        new.ref <- do.call(BiocGenerics::cbind, all.refs)

        out <- SingleR(test[all.markers,i], new.ref, genes="all", 
            label=origins, fine.tune=FALSE, prune=FALSE)

        for (j in seq_along(labs)) {
            expect_equivalent(scores[i,labs[[j]]], out$scores[,j])
        }
    }
}

test_that("combineRecomputedResults works as expected (thorough)", {
    combined <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=test,
        trained=list(train1, train2))

    REF(combined$scores, test, results=list(pred1, pred2), refs=list(ref1, ref2), subset=1:10)

    # Works for 3+ references.
    ref3 <- .mockRefData(nreps=8)
    ref3 <- scater::logNormCounts(ref3)
    ref3$label <- paste0(ref3$label, "X") # avoid problems with same column name in 'scores'.
    train3 <- trainSingleR(ref3, labels=ref3$label)
    pred3 <- classifySingleR(test, train3)

    combined <- combineRecomputedResults(
        results=list(pred1, pred2, pred3),
        test=test,
        trained=list(train1, train2, train3))
   
    REF(combined$scores, test, results=list(pred1, pred2, pred3), refs=list(ref1, ref2, ref3), subset=ncol(test) - 0:9)
})

test_that("combineRecomputedResults works as expected in SD mode", {
    train10 <- trainSingleR(ref1, labels=ref1$label, genes="sd")
    pred10 <- classifySingleR(test, train10)

    train20 <- trainSingleR(ref2, labels=ref2$label, genes="sd")
    pred20 <- classifySingleR(test, train20)

    combined <- combineRecomputedResults(
        results=list(pred10, pred20),
        test=test,
        trained=list(train10, train20))

    REF(combined$scores, test, results=list(pred10, pred20), refs=list(ref1, ref2), subset=11:20, mode="sd")
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

    # Correctly intersects the gene universes.
    ref <- combineRecomputedResults(
        results=list(pred1, pred2),
        test=test,
        trained=list(train1, train2))

    s <- sample(nrow(test))
    expect_warning(out <- combineRecomputedResults(
        results=list(pred1, pred2),
        test=test[s,],
        trained=list(train1, train2)), "differ in the universe")
    expect_identical(ref, out)
})

test_that("combineRecomputedResults is invariant to ordering", {
    ref3 <- .mockRefData(nreps=8)
    ref3 <- scater::logNormCounts(ref3)
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
