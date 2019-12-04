# This tests that combineResults works as expected.
# library(testthat); library(SingleR); source("test-combine.R")

test_that("combineResults works as expected", {
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

    combined <- combineResults(list(res1=results, res2=results2))
    expect_identical(combined$scores, cbind(scores, scores2))
    expect_identical(combined$first.labels, c("J1", "D1", "H1", "B1", "F1"))
    expect_identical(combined$labels, c("J", "D", "H", "B", "F"))
    expect_identical(combined$pruned.labels, c("j", "d", "h", "b", "f"))
    expect_equivalent(combined$orig.results$res1, results)
    expect_equivalent(combined$orig.results$res2, results2)

    combined <- combineResults(list(res1=results, res2=results))
    expect_identical(combined[,-c(1, ncol(combined))], results[,-1])
})

test_that("combineResults handles non-matching cells", {
    scores <- diag(5)
    colnames(scores) <- LETTERS[1:5]
    
    results <- DataFrame(scores=I(scores), labels=colnames(scores)[max.col(scores)])
    rownames(results) <- c("this", "is", "a", "neat", "test")

    scores2 <- diag(5)
    colnames(scores2) <- LETTERS[6:10]
    
    results2 <- DataFrame(scores=I(scores2), labels=colnames(scores2)[max.col(scores2)])
    rownames(results2) <- c("this", "is", "a", "neat", "package")

    expect_error(combined <- combineResults(list(res1=results, res2=results2)), "cell/cluster names are not identical")
})

test_that("combineResults handles non-common genes", {
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

    expect_warning(combined <- combineResults(list(res1=results, res2=results2)), "common genes are not identical")
})
