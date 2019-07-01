# This tests the classification *without* fine-tuning.
# library(testthat); library(SingleR); source("setup.R"); source("test-classify.R")

trained <- trainSingleR(training, training$label)

test_that("correlations are computed correctly by classifySingleR", {
    Q <- 0.2
    test <- classifySingleR(test, trained, fine.tune=FALSE, quant.thresh=Q)
    out <- colData(test)$SingleR

    # Computing reference correlations between test and trained.
    y <- split(seq_along(training$label), training$label)
    collected <- matrix(0, ncol(test), length(y))
    colnames(collected) <- names(y)
    genes <- trained$common.genes

    for (x in seq_along(y)) {
        ref <- cor(assay(training)[genes,y[[x]]], assay(test)[genes,], method="spearman")
        k <- max(1, round(Q * nrow(ref)))
       collected[,x] <- apply(ref, 2, FUN=function(x) sort(x, decreasing=TRUE)[k])
    }

    # Checking that they're the same.
    expect_equal(collected, out$scores[,colnames(collected)])

    # Checking that the correct label is chosen.
    expect_identical(colnames(collected)[max.col(collected)], out$labels)
})

test_that("classifySingleR behaves sensibly with very low 'quant.thresh' settings", {
    Q <- 0
    test <- classifySingleR(test, trained, fine.tune=FALSE, quant.thresh=Q)
    out <- colData(test)$SingleR

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

test_that("classifySingleR behaves sensibly with very large 'quant.thresh' settings", {
    Q <- 1
    test <- classifySingleR(test, trained, fine.tune=FALSE, quant.thresh=Q)
    out <- colData(test)$SingleR

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

test_that("classifySingleR behaves with no-variance cells", {
    sce <- test 
    assay(sce)[,1:10] <- 0

    Q <- 0.2
    sce <- classifySingleR(sce, trained, fine.tune=FALSE, quant.thresh=Q)
    out <- colData(sce)$SingleR
    expect_true(all(is.na(out2$scores[1:10,])))
    expect_true(all(is.na(out2$labels[1:10])))

    test <- classifySingleR(test, trained, fine.tune=FALSE, quant.thresh=Q)
    ref <- colData(test)$SingleR
    expect_identical(out$scores[-(1:10),], ref$scores[-(1:10),])
    expect_identical(out$labels[-(1:10)], ref$labels[-(1:10)])
})

test_that("classifySingleR behaves with silly inputs", {
    test <- classifySingleR(test[,0], trained, fine.tune=FALSE)
    out <- colData(test)$SingleR
    expect_identical(nrow(out$scores), 0L)
    expect_identical(length(out$labels), 0L)

    expect_error(classifySingleR(test[0,], trained, fine.tune=FALSE), "does not contain")
})


