# This tests the classic marker detection scheme.
# library(testthat); library(SingleR); source("setup.R"); source("test-markers.R")

library(DelayedArray)
median_by_label <- function(mat, labels) {
    ulabels <- SingleR:::.get_levels(labels)
    output <- matrix(0, nrow(mat), length(ulabels))
    rownames(output) <- rownames(mat)
    colnames(output) <- ulabels

    for (u in ulabels) {
        output[,u] <- apply(as.matrix(mat[,u==labels]), 1, median)
    }
    output
}

REF <- function(ref, labels, de.n=NULL) {
    mat <- median_by_label(ref, labels)
    if (is.null(de.n)) {
        de.n <- round(500*(2/3)^log2(ncol(mat)))
    }

    ulabels <- SingleR:::.get_levels(labels)
    collected <- list()
    for (i in ulabels) {
        subcollected <- list()

        for (j in ulabels) {
            s <- sort(mat[,i] - mat[,j], decreasing=TRUE)
            s <- s[s>0]
            subcollected[[j]] <- as.character(head(names(s), de.n))
        }
        collected[[i]] <- subcollected
    }

    collected
}

test_that("getClassicMarkers works as expected", {
    ref <- REF(logcounts(training), training$label)
    out <- getClassicMarkers(logcounts(training), training$label)
    expect_identical(ref, out)

    # Unaffected by monotonic transforms.
    out <- getClassicMarkers(training, training$label, assay.type="logcounts")
    expect_identical(ref, out)

    # Same results if the labels are reversed.
    shuffle <- rev(seq_len(ncol(training)))
    out <- getClassicMarkers(training[,shuffle], training$label[shuffle], assay.type="logcounts")
    expect_identical(ref, out)

    # Responds to a custom de.n=.
    ref <- REF(logcounts(training), training$label, de.n=20)
    out <- getClassicMarkers(training, training$label, assay.type="logcounts", de.n=20)
    expect_identical(ref, out)

    # Caps out correctly.
    ref <- REF(logcounts(training), training$label, de.n=1e6)
    out <- getClassicMarkers(training, training$label, assay.type="logcounts", de.n=1e6)
    expect_identical(ref, out)
})

test_that("getClassicMarkers works with blocking", {
    # Blocking emits the same results in the simple case.
    out <- getClassicMarkers(logcounts(training), training$label)
    out2 <- getClassicMarkers(list(logcounts(training), logcounts(training)), 
        list(training$label, training$label))
    expect_identical(out, out2)

    # Blocking is robust to training sets that don't have the labels. 
    out3 <- getClassicMarkers(list(logcounts(training), logcounts(training)[,0]), 
        list(training$label, training$label[0]))
    expect_identical(out, out3)

    # Blocking actually makes a difference.
    set.seed(100)
    shuffled <- sample(training$label)
    out2 <- getClassicMarkers(list(logcounts(training), logcounts(training)),
        list(training$label, shuffled))
    expect_false(identical(out, out2))
})
