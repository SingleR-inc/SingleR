# This tests the fine-tuning functionality.
# library(testthat); library(SingleR); source("setup.R"); source("test-fine-tune.R")

# A reference implementation in pure R.
.fine_tune_cell <- function(i, exprs, scores, references, quantile, tune.thresh, commonFUN, ...) {
    cur.exprs <- exprs[,i]
    cur.scores <- scores[i,]
    top.labels <- names(cur.scores)[cur.scores >= max(cur.scores) - tune.thresh]
    old.labels <- character(0)

    # Need to compare to old.labels, to avoid an infinite loop 
    # if the correlations are still close after fine tuning.
    while (length(top.labels) > 1L && !identical(top.labels, old.labels)) {
        common <- commonFUN(top.labels, ...)
        cur.scores <- .compute_label_scores_manual(common, top.labels, cur.exprs, references, quantile=quantile)
        old.labels <- top.labels
        cur.scores <- cur.scores[!is.na(cur.scores)]
        top.labels <- names(cur.scores)[cur.scores >= max(cur.scores) - tune.thresh] 
    }

    if (length(top.labels)==1L) {
        label <- top.labels
    } else if (length(top.labels)==0L) {
        label <- NA_character_
    } else {
        label <- names(cur.scores)[which.max(cur.scores)]
    }
    list(label=label, best=max(cur.scores), second=-sort.int(-cur.scores)[2])
}

.compute_label_scores_manual <- function(common, top.labels, cur.exprs, references, quantile) {
    cur.exprs <- cur.exprs[common]
    cur.scores <- numeric(length(top.labels))
    names(cur.scores) <- top.labels

    for (u in top.labels) {
        ref <- references[[u]]
        ref <- as.matrix(ref[common,,drop=FALSE])

        # We use a 'k'-based method for selecting the quantile, for consistency with classifySingleR.
        k <- max(1, round((1-quantile) * ncol(ref)))
        cur.cor <- cor(cur.exprs, ref, method="spearman")
        cur.scores[u] <- -sort(-cur.cor, partial=k)[k]
    }

    cur.scores
}

.fine_tune_de_genes <- function(top.labels, de.info) {
    # Finding the subset of genes (assuming 'extras' is a matrix of lists).
    all.combos <- expand.grid(top.labels, top.labels)
    Reduce(union, de.info[as.matrix(all.combos)])
}

.fine_tune_sd_genes <- function(top.labels, extras, sd.thresh) {
    sds <- DelayedMatrixStats::rowSds(extras, col=match(top.labels, colnames(extras)))
    rownames(extras)[sds > sd.thresh]
}

.fine_tune_de_ref <- function(exprs, scores, references, quantile, tune.thresh, de.info) { 
    de.info <- do.call(cbind, de.info)
    if (is.null(colnames(de.info)) || !identical(colnames(de.info), rownames(de.info))) {
        stop("marker list should be named during training")
    }
    out <- lapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell, exprs=exprs, scores=scores, 
        references=references, quantile=quantile, tune.thresh=tune.thresh, 
        commonFUN=.fine_tune_de_genes, de.info=de.info)
    do.call(mapply, c(list(FUN=c, SIMPLIFY=FALSE, USE.NAMES=FALSE), out))
}

.fine_tune_sd_ref <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh) {
    out <- lapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell, exprs=exprs, scores=scores, 
        references=references, quantile=quantile, tune.thresh=tune.thresh, 
        commonFUN=.fine_tune_sd_genes, median.mat, sd.thresh=sd.thresh)
    do.call(mapply, c(list(FUN=c, SIMPLIFY=FALSE, USE.NAMES=FALSE), out))
}

############################
############################
# Checking that our test functions are actually correct. 

test_that("fine-tuning correlation calculator is consistent with classifySingleR", {
    trained <- trainSingleR(training, training$label, genes='all')

    Q <- 0.8
    out <- .compute_label_scores_manual(assay(test)[,1], trained$original.exprs, 
        top.labels=names(trained$original.exprs), common=rownames(test), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[1,])

    Q <- 0.2
    out <- .compute_label_scores_manual(assay(test)[,2], trained$original.exprs,
        top.labels=names(trained$original.exprs), common=rownames(test), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[2,])

    Q <- 0.5
    out <- .compute_label_scores_manual(assay(test)[,3], trained$original.exprs, 
        top.labels=names(trained$original.exprs), common=rownames(test), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[3,])
})

test_that("fine-tuning correlation calculator responds to gene subsets", {
    trained <- trainSingleR(training, training$label, genes='de')

    Q <- 0.8
    out <- .compute_label_scores_manual(assay(test)[,1], trained$original.exprs, 
        common=trained$common.genes, top.labels=names(trained$original.exprs), quantile=Q)
    ref <- classifySingleR(test, trained, fine.tune=FALSE, quantile=Q)
    expect_equal(out, ref$scores[1,])

    out <- .compute_label_scores_manual(assay(test)[,2], trained$original.exprs, 
        common=trained$common.genes, top.labels=names(trained$original.exprs)[1:2], quantile=Q)
    expect_equal(out, ref$scores[2,1:2])
})

test_that("fine-tuning DE marker selection works", {
    trained <- trainSingleR(training, training$label, genes='de')

    stuff <- as.list(trained$search$extra)
    de.info <- do.call(cbind, stuff)
    out <- .fine_tune_de_genes(c("A", "B"), de.info)
    expect_identical(sort(out), sort(union(stuff$A$B, stuff$B$A)))

    # Throwing if the constructed matrix doesn't have dimnames.
    Q <- 0.8
    thresh <- 0.05
    expect_error(.fine_tune_de_ref(assay(test), ref$scores, trained$original.exprs,
         quantile=Q, tune.thresh=thresh, de.info=unname(trained$search$extra)), "named")
    expect_error(.fine_tune_de_ref(assay(test), ref$scores, trained$original.exprs,
         quantile=Q, tune.thresh=thresh, de.info=lapply(trained$search$extra, unname)), "named")
})

test_that("fine-tuning SD selection works", {
    trained <- trainSingleR(training, training$label, genes='sd', sd.thresh=0.5)

    mat <- trained$search$extra
    lab <- c("A", "B", "C")
    ref <- DelayedMatrixStats::rowSds(mat[,lab])

    # Responds to various settings correctly
    # (difficult to test exactly without regurgitating the code).
    out <- .fine_tune_sd_genes(lab, mat, 1)
    expect_true(min(ref[match(out, rownames(mat))]) > max(ref[-match(out, rownames(mat))]))

    out2 <- .fine_tune_sd_genes(lab, mat, 0.5)
    expect_true(length(out2) > length(out))
    expect_true(min(ref[match(out2, rownames(mat))]) > max(ref[-match(out2, rownames(mat))]))
})

###################################
###################################

# Splitting each label into two for some more variety.
subset <- sample(1:2, ncol(training), replace=TRUE)
training$label <- paste0(training$label, subset)

test_that("fine-tuning by DE runs without errors", {
    trained <- trainSingleR(training, training$label, genes='de')
    pred <- classifySingleR(test, trained, fine.tune=FALSE)

    # Testing minor offsets to avoid problems with numerical precision.
    for (Q in c(0, 0.21, 0.51, 0.81, 1)) { 
        for (thresh in c(0, 0.05, 0.1)) {
            tuned <- SingleR:::.fine_tune_de(assay(test), pred$scores, trained$original.exprs, 
                 quantile=Q, tune.thresh=thresh, de.info=trained$search$extra)
            ref <- .fine_tune_de_ref(assay(test), pred$scores, trained$original.exprs, 
                 quantile=Q, tune.thresh=thresh, de.info=trained$search$extra)

            expect_equal(colnames(pred$scores)[tuned[[1]]+1], ref[[1]])
            expect_equal(tuned[[2]], ref[[2]])
            expect_equal(tuned[[3]], unname(ref[[3]]))
        }
    }

    # Sanity checking of the dimensions and output.
    Q <- 0.8
    thresh <- 0.05
    tuned <- SingleR:::.fine_tune_de(assay(test), pred$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, de.info=trained$search$extra)

    expect_identical(lengths(tuned), rep(nrow(pred), length(tuned)))
    is.diff <- colnames(pred$scores)[tuned[[1]]+1]!=pred$labels
    expect_true(any(is.diff))

    maxed <- pred$scores[cbind(seq_len(nrow(pred)), max.col(pred$scores))]
    expect_false(any(is.diff & rowSums(pred$scores >= maxed - thresh)==1L))

    # Works with parallelization.
    multi <- SingleR:::.fine_tune_de(assay(test), pred$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, de.info=trained$search$extra, 
         BPPARAM=BiocParallel::SnowParam(3))
    expect_identical(multi, tuned)
})

test_that("fine-tuning by SD runs without errors", {
    trained <- trainSingleR(training, training$label, genes='sd', sd.thresh=0.5) # turning down the threshold.
    pred <- classifySingleR(test, trained, fine.tune=FALSE)

    # Minor offsets to avoid problems with numerical precision.
    for (Q in c(0, 0.21, 0.51, 0.81, 1)) { 
        for (thresh in c(0, 0.05, 0.1)) {
            tuned <- SingleR:::.fine_tune_sd(assay(test), pred$scores, trained$original.exprs, 
                 quantile=Q, tune.thresh=thresh, median.mat=trained$search$extra, sd.thresh=0.5)
            ref <- .fine_tune_sd_ref(assay(test), pred$scores, trained$original.exprs, 
                 quantile=Q, tune.thresh=thresh, median.mat=trained$search$extra, sd.thresh=0.5)

            expect_equal(colnames(pred$scores)[tuned[[1]]+1], ref[[1]])
            expect_equal(tuned[[2]], ref[[2]])
            expect_equal(tuned[[3]], unname(ref[[3]]))
        }
    }

    # Sanity checking of the dimensions and output.
    Q <- 0.8
    thresh <- 0.05
    tuned <- SingleR:::.fine_tune_sd(assay(test), pred$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, median.mat=trained$search$extra, sd.thresh=1)

    expect_identical(lengths(tuned), rep(nrow(pred), length(tuned)))
    is.diff <- colnames(pred$scores)[tuned[[1]]+1]!=pred$labels
    expect_true(any(is.diff))

    maxed <- pred$scores[cbind(seq_len(nrow(pred)), max.col(pred$scores))]
    expect_false(any(is.diff & rowSums(pred$scores >= maxed - thresh)==1L))

    # Works with parallelization.
    multi <- SingleR:::.fine_tune_sd(assay(test), pred$scores, trained$original.exprs, 
         quantile=Q, tune.thresh=thresh, median.mat=trained$search$extra, sd.thresh=1,
         BPPARAM=BiocParallel::SnowParam(3))
    expect_identical(multi, tuned)
})

test_that("fine-tuning handles the edge cases sensibly", {
    # Only one label available:
    trained <- trainSingleR(training[,1], training$label[1])
    pred <- classifySingleR(test, trained, fine.tune=FALSE)

    tuned <- SingleR:::.fine_tune_de(assay(test), pred$scores, trained$original.exprs, 
         quantile=0.5, tune.thresh=0.1, de.info=trained$search$extra)
    expect_true(all(is.na(tuned[[3]])))
    expect_false(any(is.na(tuned[[2]])))

    # No labels available:
    trained <- trainSingleR(training[,0], training$label[0])
    pred <- classifySingleR(test, trained, fine.tune=FALSE)

    tuned <- SingleR:::.fine_tune_de(assay(test), pred$scores, trained$original.exprs, 
         quantile=0.5, tune.thresh=0.1, de.info=trained$search$extra)
    expect_true(all(is.na(tuned[[3]])))
    expect_true(all(is.na(tuned[[2]])))
})

###################################
###################################

test_that("fine-tuning works in the body of the function", {
    trained <- trainSingleR(training, training$label, genes='de')
    ref <- classifySingleR(test, trained)
    expect_true("first.labels" %in% colnames(ref))

    trained <- trainSingleR(training, training$label, genes='sd', sd.thresh=0.5)
    ref <- classifySingleR(test, trained)
    expect_true("first.labels" %in% colnames(ref))

    # Stress-checking edge cases.
    trained <- trainSingleR(training[,0], training$label[0])
    pred <- classifySingleR(test, trained)
    expect_identical(sum(is.na(pred$labels)), ncol(test))

    trained <- trainSingleR(training[0,], training$label)
    expect_error(pred <- classifySingleR(test, trained), NA)
})
