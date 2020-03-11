# This tests that combineCommonResults works as expected.
# library(testthat); library(SingleR); source("test-combine.R")

test_that("combineCommonResults works as expected", {
    scores <- rbind(
        c(0,0,0,0,1),
        c(0,0,0,2,0),
        c(0,0,1,0,0),
        c(0,2,0,0,0),
        c(1,0,0,0,0)
    )
    colnames(scores) <- LETTERS[1:5]
    
    labs <- colnames(scores)[max.col(scores)]
    results <- DataFrame(scores=I(scores), first.labels=paste0(labs, "1"),
        labels=labs, pruned.labels=tolower(labs))
    rownames(results) <- sprintf("CELL_%s", seq_len(nrow(results)))

    scores2 <- rbind(
        c(0,0,0,0,2),
        c(0,0,0,1,0),
        c(0,0,2,0,0),
        c(0,1,0,0,0),
        c(2,0,0,0,0)
    )
    colnames(scores2) <- LETTERS[6:10]
    
    labs2 <- colnames(scores2)[max.col(scores2)]
    results2 <- DataFrame(scores=I(scores2), first.labels=paste0(labs2, "1"),
        labels=labs2, pruned.labels=tolower(labs2))
    rownames(results2) <- sprintf("CELL_%s", seq_len(nrow(results2)))

    combined <- combineCommonResults(list(res1=results, res2=results2))
    expect_identical(combined$scores, cbind(scores, scores2))
    expect_identical(combined$first.labels, c("J1", "D1", "H1", "B1", "F1"))
    expect_identical(combined$labels, c("J", "D", "H", "B", "F"))
    expect_identical(combined$pruned.labels, c("j", "d", "h", "b", "f"))
    expect_equivalent(combined$orig.results$res1, results)
    expect_equivalent(combined$orig.results$res2, results2)

    combined <- combineCommonResults(list(res1=results, res2=results))
    standard <- c("first.labels", "labels", "pruned.labels")
    expect_identical(combined$scores, cbind(results$scores, results$scores))
    expect_identical(combined[,standard], results[,standard])
})

test_that("combineCommonResults handles non-matching cells", {
    scores <- diag(5)
    colnames(scores) <- LETTERS[1:5]
    results <- DataFrame(scores=I(scores), labels=colnames(scores)[max.col(scores)])

    scores2 <- diag(5)
    colnames(scores2) <- LETTERS[6:10]
    results2 <- DataFrame(scores=I(scores2), labels=colnames(scores2)[max.col(scores2)])

    expect_error(combineCommonResults(list(res1=results[1:2,], res2=results2)), "numbers .*are not identical")

    rownames(results) <- c("this", "is", "a", "neat", "test")
    rownames(results2) <- c("this", "is", "a", "neat", "package")
    expect_error(combineCommonResults(list(res1=results, res2=results2)), "cell/cluster names.*are not identical")
})

test_that("combineCommonResults handles non-common genes", {
    scores <- diag(5)
    colnames(scores) <- LETTERS[1:5]
    
    results <- DataFrame(scores=I(scores), labels=colnames(scores)[max.col(scores)])
    metadata(results)$common.genes <- c("i", "am", "a", "common", "gene")
    rownames(results) <- sprintf("GENE_%s", seq_len(nrow(results)))

    scores2 <- diag(5)
    colnames(scores2) <- LETTERS[6:10]
    
    results2 <- DataFrame(scores=I(scores2), labels=colnames(scores2)[max.col(scores2)])
    metadata(results2)$common.genes <- c("i", "am", "an", "uncommon", "gene")
    rownames(results2) <- sprintf("GENE_%s", seq_len(nrow(results2)))

    expect_warning(combined <- combineCommonResults(list(res1=results, res2=results2)), "common genes are not identical")
})

####################################################################

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

test_that("combineRecomputedResults works as expected", {
    # Combining results with recomputation of scores.
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

    # A more thorough check.
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
