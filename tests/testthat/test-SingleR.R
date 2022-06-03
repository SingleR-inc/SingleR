# End-to-end testing of SingleR on a dummy training/test dataset that should be near-perfectly assigned.
# library(testthat); library(SingleR); source("setup.R"); source("test-SingleR.R")

test_that("SingleR works in DE mode", {
    out <- SingleR(test=test, ref=training, labels=training$label, genes="de")
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.95)
})

test_that("SingleR works with custom gene selection", {
    all.labs <- sort(unique(training$label))
    collected <- rep(list(tail(rownames(training), 100)), length(all.labs))
    names(collected) <- all.labs

    more.collected <- rep(list(collected), length(all.labs))
    names(more.collected) <- all.labs

    out <- SingleR(test=test, ref=training, labels=training$label, genes=more.collected)
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.95)
    
    # We should get, in this case, the same result with a list of vectors.
    out2 <- SingleR(test=test, ref=training, labels=training$label, genes=collected)
    expect_identical(collected, lapply(metadata(out2)$de.genes, unlist, use.names=FALSE))
    metadata(out) <- metadata(out2) <- list()
    expect_identical(out, out2)
})

test_that("SingleR works when genes are not the same between test and training", {
    out <- SingleR(test=test[1:800,], ref=training[200:1000,], labels=training$label)
    ref <- SingleR(test=test[200:800,], ref=training[200:800,], labels=training$label)
    expect_identical(out, ref)

    expect_error(SingleR(test=test[1:200,], ref=training[800:1000,], labels=training$label), "no common genes")
})

library(Matrix)
test_that("SingleR works with non-ordinary matrices", {
    test.s <- test
    training.s <- training
    assay(test.s) <- as(assay(test.s), "dgCMatrix")
    assay(training.s) <- as(assay(training.s), "dgCMatrix")

    out <- SingleR(test=test.s, ref=training.s, labels=training.s$label)
    ref <- SingleR(test=test, ref=training, labels=training$label)
    expect_identical(out, ref)

    # Works with clustering mode.
    out <- SingleR(test=test.s, ref=training.s, labels=training.s$label, clusters=test.s$label)
    ref <- SingleR(test=test, ref=training, labels=training$label, clusters=test$label)
    expect_identical(out, ref)
})

test_that("SingleR works with multiple references", {
    # Handles mismatching row names.
    chosen0 <- sample(rownames(training), 900)
    chosen1 <- sample(rownames(training), 900)
    chosen2 <- sample(rownames(training), 900)

    # Works with recomputation.
    out <- SingleR(test[chosen0,], list(training[chosen1,], training[chosen2,]), 
        list(training$label, training$label), recompute=TRUE)
    inter <- Reduce(intersect, list(chosen0, chosen1, chosen2))
    ref <- SingleR(test[inter,], list(training[inter,], training[inter,]), 
        list(training$label, training$label), recompute=TRUE)

    out$reference <- ref$reference <- NULL # basically tied anyway.
    expect_identical(out, ref)
})

test_that("SingleR handles data.frame inputs", {
    set.seed(10)
    ref1 <- SingleR(test=test, ref=training, labels=training$label)
    set.seed(10)
    ref2 <- SingleR(test=data.frame(logcounts(test)), ref=data.frame(logcounts(training)), labels=training$label)
    
    rownames(ref2) <- NULL # as the data.frame coercion changes the cell's names.
    expect_identical(ref1, ref2)

    tmp <- data.frame(logcounts(test))
    tmp[,1] <- as.character(tmp[,1])
    expect_error(SingleR(test=tmp, ref=data.frame(logcounts(training)), labels=training$label, genes="de"), "numeric")
})

test_that("SingleR handles NAs in the labels", {
    populate <- rbinom(length(training$label), 1, 0.2)==1
    training$label[populate] <- NA

    set.seed(10)
    ref1 <- SingleR(test=test, ref=training, labels=training$label)
    set.seed(10)
    ref2 <- SingleR(test=test, ref=training[,!populate], labels=training$label[!populate])

    expect_identical(ref1, ref2)
})
