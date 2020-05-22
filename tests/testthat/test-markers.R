# This tests the classic marker detection scheme.
# library(testthat); library(SingleR); source("setup.R"); source("test-markers.R")

REF <- function(ref, labels, de.n=NULL) {
    mat <- SingleR:::.median_by_label(ref, labels)
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

    out <- getClassicMarkers(training, training$label, assay.type="logcounts")
    expect_identical(ref, out)

    ref <- REF(logcounts(training), training$label, de.n=20)
    out <- getClassicMarkers(training, training$label, assay.type="logcounts", de.n=20)
    expect_identical(ref, out)
})

test_that("getClassicMarkers works with blocking", {
    # Blocking emits the same results in the simple case.
    out <- getClassicMarkers(logcounts(training), training$label)
    out2 <- getClassicMarkers(list(logcounts(training), logcounts(training)), 
        list(training$label, training$label))
    expect_identical(out, out2)

    # Blocking actually makes a difference.
    set.seed(100)
    shuffled <- sample(training$label)
    out2 <- getClassicMarkers(list(logcounts(training), logcounts(training)),
        list(training$label, shuffled))
    expect_false(identical(out, out2))
})
