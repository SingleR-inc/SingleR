# This tests the classification *without* fine-tuning.
# library(testthat); library(SingleR); source("setup.R"); source("test-classify.R")

trained <- trainSingleR(training, training$label)

test_that("correlations are computed correctly by classifySingleR", {
    Q <- 0.8
    out <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)

    # Computing reference correlations between test and trained.
    y <- split(seq_along(training$label), training$label)
    collected <- matrix(0, ncol(test), length(y))
    colnames(collected) <- names(y)
    genes <- trained$common.genes

    for (x in seq_along(y)) {
        ref <- cor(assay(training)[genes,y[[x]]], assay(test)[genes,], method="spearman")
        k <- max(1, round((1-Q) * nrow(ref)))
       collected[,x] <- apply(ref, 2, FUN=function(x) sort(x, decreasing=TRUE)[k])
    }

    # Checking that they're the same.
    expect_equal(collected, out$scores[,colnames(collected)])

    # Checking that the correct label is chosen.
    expect_identical(colnames(collected)[max.col(collected)], out$labels)
})

test_that("classifySingleR behaves sensibly with very low 'quantile' settings", {
    Q <- 0
    out <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)

    # Computing reference correlations between test and trained.
    y <- split(seq_along(training$label), training$label)
    collected <- matrix(0, ncol(test), length(y))
    colnames(collected) <- names(y)
    genes <- trained$common.genes

    for (x in seq_along(y)) {
        ref <- cor(assay(training)[genes,y[[x]]], assay(test)[genes,], method="spearman")
        collected[,x] <- apply(ref, 2, FUN=min)
    }

    # Checking that they're the same.
    expect_equal(collected, out$scores[,colnames(collected)])
    expect_identical(colnames(collected)[max.col(collected)], out$labels)
})

test_that("classifySingleR behaves sensibly with very large 'quantile' settings", {
    Q <- 1
    out <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)

    # Computing reference correlations between test and trained.
    y <- split(seq_along(training$label), training$label)
    collected <- matrix(0, ncol(test), length(y))
    colnames(collected) <- names(y)
    genes <- trained$common.genes

    for (x in seq_along(y)) {
        ref <- cor(assay(training)[genes,y[[x]]], assay(test)[genes,], method="spearman")
        collected[,x] <- apply(ref, 2, FUN=max)
    }

    # Checking that they're the same.
    expect_equal(collected, out$scores[,colnames(collected)])
    expect_identical(colnames(collected)[max.col(collected)], out$labels)
})

test_that("classifySingleR behaves with no-variance cells", {
    sce <- test 
    logcounts(sce)[,1:10] <- 0

    Q <- 0.2
    out <- classifySingleR(sce, trained, fine.tune=FALSE, quantile=Q)
    expect_true(all(abs(out$scores[1:10,] - 0.5) < 1e-8)) # works out to 0.5, as a mathematical oddity.

    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_identical(out$scores[-(1:10),], ref$scores[-(1:10),])
    expect_identical(out$labels[-(1:10)], ref$labels[-(1:10)])
})

test_that("classifySingleR behaves with missing values", {
    # Can't just set the first entry to NA, as we need to ensure 
    # that the test set contains a superset of genes in the training set.
    sce <- BiocGenerics::rbind(test[1,], test)
    logcounts(sce)[1,1] <- NA

    Q <- 0.8
    expect_warning(out <- classifySingleR(sce, trained, fine.tune=FALSE, quantile=Q), 'missing values')
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_identical(out, ref)
})

test_that("classifySingleR behaves with silly inputs", {
    out <- classifySingleR(test[,0], trained, fine.tune=FALSE)
    expect_identical(nrow(out$scores), 0L)
    expect_identical(length(out$labels), 0L)
    expect_error(classifySingleR(test[0,], trained, fine.tune=FALSE), "does not contain")

    trained <- trainSingleR(training[,0], training$label[0])
    pred <- classifySingleR(test, trained, fine.tune=FALSE)
    expect_identical(ncol(pred$scores), 0L)
    expect_true(all(is.na(pred$labels)))
})
