# End-to-end testing of SingleR on a dummy training/test dataset that should be near-perfectly assigned.
# library(testthat); library(SingleR); source("setup.R"); source("test-SingleR.R")

test_that("SingleR works in DE mode", {
    out <- SingleR(test=test, ref=training, labels=training$label, genes="de")
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.95)
})

test_that("SingleR works in SD mode", {
    out <- SingleR(test=test, ref=training, labels=training$label, genes="sd")
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.85) # not as good as 'de'.
})

test_that("SingleR works with custom gene selection", {
    all.labs <- unique(training$label)
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

    out <- SingleR(test=test.s, ref=training.s, labels=training.s$label, clusters=test.s$label, method="cluster")
    ref <- SingleR(test=test, ref=training, labels=training$label, clusters=test$label, method="cluster")
    expect_identical(out, ref)
})

test_that("SingleR works with multiple references", {
    out <- SingleR(test, list(training, training), list(training$label, training$label), recompute=FALSE)
    ref <- SingleR(test, training, training$label)
    expect_identical(out$scores, cbind(ref$scores, ref$scores))

    # Handles mismatching rownames.
    chosen0 <- sample(rownames(training), 900)
    chosen1 <- sample(rownames(training), 900)
    chosen2 <- sample(rownames(training), 900)

    out <- SingleR(test[chosen0,], list(training[chosen1,], training[chosen2,]), 
        list(training$label, training$label), recompute=FALSE)
    inter <- Reduce(intersect, list(chosen0, chosen1, chosen2))
    ref <- SingleR(test[inter,], list(training[inter,], training[inter,]), 
        list(training$label, training$label), recompute=FALSE)

    out$reference <- ref$reference <- NULL # basically tied anyway.
    expect_identical(out, ref)

    # Works with recomputation.
    out <- SingleR(test[chosen0,], list(training[chosen1,], training[chosen2,]), 
        list(training$label, training$label), recompute=TRUE)
    inter <- Reduce(intersect, list(chosen0, chosen1, chosen2))
    ref <- SingleR(test[inter,], list(training[inter,], training[inter,]), 
        list(training$label, training$label), recompute=TRUE)

    out$reference <- ref$reference <- NULL # basically tied anyway.
    expect_identical(out, ref)
})
