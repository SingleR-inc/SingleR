# This tests the fine-tuning functionality.
# library(testthat); library(SingleR); source("setup.R"); source("test-fine-tune.R")

test_that("fine-tuning correlation calculator is consistent with classifySingleR", {
    trained <- trainSingleR(training, training$label, genes='all')

    Q <- 0.8
    out <- SingleR:::.compute_label_scores_manual(assay(test)[,1], trained$original.exprs, 
        top.labels=names(trained$original.exprs), common=rownames(test), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[1,])

    Q <- 0.2
    out <- SingleR:::.compute_label_scores_manual(assay(test)[,2], trained$original.exprs,
        top.labels=names(trained$original.exprs), common=rownames(test), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[2,])

    Q <- 0.5
    out <- SingleR:::.compute_label_scores_manual(assay(test)[,3], trained$original.exprs, 
        top.labels=names(trained$original.exprs), common=rownames(test), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[3,])
})

test_that("fine-tuning correlation calculator responds to gene subsets", {
    trained <- trainSingleR(training, training$label, genes='de')

    Q <- 0.8
    out <- SingleR:::.compute_label_scores_manual(assay(test)[,1], trained$original.exprs, 
        common=trained$common.genes, top.labels=names(trained$original.exprs), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[1,])

    out <- SingleR:::.compute_label_scores_manual(assay(test)[,2], trained$original.exprs, 
        common=trained$common.genes, top.labels=names(trained$original.exprs)[1:2], quantile=Q)
    expect_equal(out, ref$scores[2,1:2])
})

###################################
###################################

# Splitting each label into two for some more variety.
subset <- sample(1:2, ncol(training), replace=TRUE)
training$label <- paste0(training$label, subset)

test_that("fine-tuning by DE runs without errors", {
    trained <- trainSingleR(training, training$label, genes='de')
    ref <- classifySingleR(test, trained, fine.tune=FALSE)

    Q <- 0.8
    thresh <- 0.05
    tuned <- SingleR:::.fine_tune_de(assay(test), ref$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, de.info=trained$search$extra)

    expect_identical(length(tuned), nrow(ref))
    is.diff <- tuned!=ref$labels
    expect_true(any(is.diff))

    maxed <- ref$scores[cbind(seq_len(nrow(ref)), max.col(ref$scores))]
    expect_false(any(is.diff & rowSums(ref$scores >= maxed - thresh)==1L))
})

test_that("fine-tuning DE marker selection works", {
    trained <- trainSingleR(training, training$label, genes='de')

    stuff <- as.list(trained$search$extra)
    de.info <- do.call(cbind, stuff)
    out <- SingleR:::.fine_tune_de_genes(c("A1", "B2"), de.info)
    expect_identical(sort(out), sort(union(stuff$A1$B2, stuff$B2$A1)))

    # Throwing if the constructed matrix doesn't have dimnames.
    Q <- 0.8
    thresh <- 0.05
    expect_error(SingleR:::.fine_tune_de(assay(test), ref$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, de.info=unname(trained$search$extra)), "named")
    expect_error(SingleR:::.fine_tune_de(assay(test), ref$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, de.info=lapply(trained$search$extra, unname)), "named")
})

test_that("fine-tuning by SD runs without errors", {
    trained <- trainSingleR(training, training$label, genes='sd')
    ref <- classifySingleR(test, trained, fine.tune=FALSE)

    Q <- 0.8
    thresh <- 0.05
    tuned <- SingleR:::.fine_tune_sd(assay(test), ref$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, median.mat=trained$search$extra, sd.thresh=1)

    expect_identical(length(tuned), nrow(ref))
    is.diff <- tuned!=ref$labels
    expect_true(any(is.diff))

    maxed <- ref$scores[cbind(seq_len(nrow(ref)), max.col(ref$scores))]
    expect_false(any(is.diff & rowSums(ref$scores >= maxed - thresh)==1L))
})

test_that("fine-tuning SD selection works", {
    trained <- trainSingleR(training, training$label, genes='sd')

    mat <- trained$search$extra
    lab <- c("A1", "A2", "B1", "B2")
    ref <- DelayedMatrixStats::rowSds(mat[,lab])

    # Responds to various settings correctly
    # (difficult to test exactly without regurgitating the code).
    out <- SingleR:::.fine_tune_sd_genes(lab, mat, 1, sd.n=100)
    expect_true(min(ref[match(out, rownames(mat))]) > max(ref[-match(out, rownames(mat))]))

    out2 <- SingleR:::.fine_tune_sd_genes(lab, mat, 1, sd.n=200)
    expect_true(length(out2) > length(out))
    expect_true(min(ref[match(out2, rownames(mat))]) > max(ref[-match(out2, rownames(mat))]))

    out3 <- SingleR:::.fine_tune_sd_genes(lab, mat, 0.5, sd.n=100)
    expect_true(length(out3) > length(out))
    expect_true(min(ref[match(out3, rownames(mat))]) > max(ref[-match(out3, rownames(mat))]))
})

###################################
###################################

test_that("fine-tuning works in the body of the function", {
    trained <- trainSingleR(training, training$label, genes='de')
    ref <- classifySingleR(test, trained)
    expect_true("first.labels" %in% colnames(ref))

    trained <- trainSingleR(training, training$label, genes='sd')
    ref <- classifySingleR(test, trained)
    expect_true("first.labels" %in% colnames(ref))
})
