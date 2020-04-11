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

test_that("combineRecomputedResults works as expected (thorough)", {
    combined <- combineRecomputedResults(
        results=list(pred1, pred2), 
        test=test,
        trained=list(train1, train2))
    
    for (i in seq_len(10)) {
        lab1 <- pred1[i,"labels"]
        lab2 <- pred2[i,"labels"]
        markers1 <- metadata(pred1)$de.genes[[lab1]]
        markers2 <- metadata(pred2)$de.genes[[lab2]]
        all.markers <- union(unlist(markers1), unlist(markers2))

        keep1 <- ref1$label==lab1
        keep2 <- ref2$label==lab2
        origins <- rep(1:2, c(sum(keep1), sum(keep2)))
        new.ref <- BiocGenerics::cbind(ref1[all.markers, keep1, drop=FALSE], ref2[all.markers, keep2,drop=FALSE])
        out <- SingleR(test[all.markers,i], new.ref, genes="all", 
            label=origins, fine.tune=FALSE, prune=FALSE)

        expect_equivalent(combined$scores[i,lab1], out$scores[,1])
        expect_equivalent(combined$scores[i,lab2], out$scores[,2])
    }
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

    for (i in 10+seq_len(10)) {
        lab1 <- pred10[i,"labels"]
        lab2 <- pred20[i,"labels"]
        markers1 <- metadata(pred10)$common.genes
        markers2 <- metadata(pred20)$common.genes
        all.markers <- union(unlist(markers1), unlist(markers2))

        keep1 <- ref1$label==lab1
        keep2 <- ref2$label==lab2
        origins <- rep(1:2, c(sum(keep1), sum(keep2)))
        new.ref <- BiocGenerics::cbind(ref1[all.markers, keep1, drop=FALSE], ref2[all.markers, keep2, drop=FALSE])
        out <- SingleR(test[all.markers,i], new.ref, genes="all",
            label=origins, fine.tune=FALSE, prune=FALSE)

        expect_equivalent(combined$scores[i,lab1], out$scores[,1])
        expect_equivalent(combined$scores[i,lab2], out$scores[,2])
    }
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
