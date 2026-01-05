# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-score-dist.R")

# pred1 = no pruned.labels
pred1 <- SingleR(test=test, ref=training, labels=training$label, genes="de")

# pred2 = 16 pruned labels
pred2 <- pred1
pred2$pruned.labels[seq(5,80,5)]<-NA

#######################################
### Single-reference checks, scores ###
#######################################

test_that("plotScoreDistribution works for no pruned, one label", {
    expect_s3_class(
        plotScoreDistribution(results = pred1, labels.use = "A"),
        "ggplot")

    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A", dots.on.top = TRUE),
        "ggplot")
})

test_that("plotScoreDistribution works for no pruned, many labels", {
    expect_s3_class(
        plotScoreDistribution(results = pred1),
        "ggplot")
})

test_that("plotScoreDistribution works with pruned, single label", {
    expect_s3_class(
        plotScoreDistribution(results = pred2, labels.use = "A"),
        "ggplot")

    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A", dots.on.top = TRUE),
        "ggplot")
})

test_that("plotScoreDistribution works with pruned, multiple labels", {
    expect_s3_class(
        plotScoreDistribution(results = pred2),
        "ggplot")
})

test_that("plotScoreDistribution warns if no 'labels.use' are present", {
    # Should give message but still output plot
    expect_warning(plotScoreDistribution(results = pred1,
        labels.use = c("a")),
        "ignoring 'labels.use'")
})

#######################################
### Single-reference checks, delta  ###
#######################################

test_that("plotDeltaDistribution works for no pruned, one label", {
    expect_s3_class(
        plotDeltaDistribution(results = pred1, labels.use = "A"),
        "ggplot")

    expect_s3_class(
        plotDeltaDistribution(
            results = pred1, show = "delta.next", labels.use = "A"),
        "ggplot")

    expect_s3_class(
        plotDeltaDistribution(
            results = pred1, labels.use = "A", dots.on.top = TRUE),
        "ggplot")
})

test_that("plotDeltaDistribution works for no pruned, many labels", {
    expect_s3_class(
        plotDeltaDistribution(results = pred1),
        "ggplot")

    expect_s3_class(
        plotDeltaDistribution(results = pred1, show = "delta.next"),
        "ggplot")
})

test_that("plotDeltaDistribution works with pruned, single label", {
    expect_s3_class(
        plotDeltaDistribution(results = pred2, labels.use = "A"),
        "ggplot")

    expect_s3_class(
        plotDeltaDistribution(
            results = pred2, labels.use = "A", show = "delta.next"),
        "ggplot")

    expect_s3_class(
        plotDeltaDistribution(
            results = pred2, labels.use = "A", dots.on.top = TRUE),
        "ggplot")
})

test_that("plotDeltaDistribution works with pruned, multiple labels", {
    expect_s3_class(
        plotDeltaDistribution(results = pred2),
        "ggplot")

    expect_s3_class(
        plotDeltaDistribution(results = pred2, show="delta.next"),
        "ggplot")
})

test_that("plotDeltaDistribution warns if no 'labels.use' are present", {
    # Should give message but still output plot
    expect_warning(plotDeltaDistribution(results = pred1,
        labels.use = c("a")),
        "ignoring 'labels.use'")
})

#######################################
### Prep for multi-reference checks ###
#######################################

ref <- .mockRefData(nreps=8)
ref1 <- ref[,1:4%%4==0]
ref1 <- ref1[,sample(ncol(ref1))]
ref1 <- scrapper::normalizeRnaCounts.se(ref1)

ref2 <- ref[,1:4%%4!=0]
ref2 <- ref2[,sample(ncol(ref2))]
ref2 <- scrapper::normalizeRnaCounts.se(ref2)

ref2$label <- tolower(ref2$label)

combined <- SingleR(
    test, ref = list(smallRef = ref1, largeRef = ref2),
    labels = list(ref1$label, ref2$label))

combined_prunedRef1 <- combined
combined_prunedRef1$orig.results$smallRef$pruned.labels[1:3%%3==0] <- NA_character_

ref1.pruned <- is.na(combined_prunedRef1$orig.results$smallRef$pruned.labels)
ref1.title <- "smallRef"
ref2.title <- "largeRef"

######################################
### Multi-reference checks, scores ###
######################################

test_that("plotScoreDistribution works for multi-ref runs - combined", {
    expect_s3_class(plotScoreDistribution(results = combined, references = 0),
        "ggplot")
})

test_that("plotScoreDistribution works for multi-ref runs - individual", {
    expect_s3_class(plotScoreDistribution(results = combined, references = 1),
        "ggplot")
})

test_that("plotScoreDistribution works for multi-ref runs - multiple", {
    expect_s3_class(plotScoreDistribution(results = combined, references = 0:1),
        "gtable")

    expect_s3_class(plotScoreDistribution(results = combined, references = NULL),
        "gtable")

    expect_equal(
        length(
            plotScoreDistribution(results = combined, grid.vars = NULL,
                references = NULL)),
        length(
            plotScoreDistribution(results = combined, grid.vars = NULL,
                references = 0:2)
        )
    )
})

test_that("plotScoreDistribution works with grid.vars control", {
    expect_s3_class(plotScoreDistribution(results = combined, references = NULL,
        grid.vars = NULL)[[1]],
        "ggplot")

    expect_s3_class(plotScoreDistribution(results = combined, references = NULL,
        grid.vars = list(ncol = 2)),
        "gtable")
})

test_that("plotScoreDistribution works with other typical adjustments", { 
    expect_s3_class(plotScoreDistribution(results = combined, labels.use = c("A", "a")),
        "gtable")

    expect_s3_class(plotScoreDistribution(results = combined, size = 3),
        "gtable")

    expect_s3_class(plotScoreDistribution(results = combined, ncol = 10),
        "gtable")

    expect_s3_class(plotScoreDistribution(results = combined, dots.on.top = FALSE),
        "gtable")
})

######################################
### Multi-reference checks, deltas ###
######################################

test_that("plotDeltaDistribution works for multi-ref runs - combined", {
    expect_error(plotDeltaDistribution(results = combined,
        references = 0, show = "delta.med"),
        "deltas cannot be shown")

    expect_error(plotDeltaDistribution(results = combined,
        references = 0, show = "delta.next"),
        "deltas cannot be shown")
})

test_that("plotDeltaDistribution works for multi-ref runs - individual", {
    expect_s3_class(plotDeltaDistribution(results = combined,
        references = 1, show = "delta.med"),
        "ggplot")

    expect_s3_class(plotDeltaDistribution(results = combined,
        references = 1, show = "delta.next"),
        "ggplot")
})

test_that("plotDeltaDistribution can be made for multi-ref runs - multiple", {
    expect_error(plotDeltaDistribution(results = combined,
        references = 0:1),
        "deltas cannot be shown")

    expect_s3_class(plotDeltaDistribution(results = combined,
        references = NULL),
        "gtable")

    # when references = NULL combined plot should not be requested.
    expect_equal(
        length(
            plotDeltaDistribution(results = combined, grid.vars = NULL,
                show = "delta.med", references = NULL)),
        length(
            plotDeltaDistribution(results = combined, grid.vars = NULL,
                show = "delta.med", references = 1:2))
        )
})

test_that("plotDeltaDistribution responds to chosen.only", {
    # Individual references
    expect_s3_class(plotDeltaDistribution(results = combined_prunedRef1, references = 1,
        chosen.only = FALSE, show = "delta.med"),
        "ggplot")

    # All references
    expect_s3_class(plotDeltaDistribution(results = combined_prunedRef1,
        chosen.only = FALSE, references = NULL, show = "delta.med"),
        "gtable")
})

test_that("plotDeltaDistribution responds to grid.vars control", {
    expect_s3_class(plotDeltaDistribution(results = combined, references = NULL,
        grid.vars = NULL)[[1]],
        "ggplot")

    expect_s3_class(plotDeltaDistribution(results = combined, references = NULL,
        grid.vars = list(ncol = 2)),
        "gtable")
})

test_that("plotDeltaDistribution responds to other typical adjustments", {
    expect_s3_class(plotDeltaDistribution(results = combined,
        labels.use = c("A", "a")),
        "gtable")

    expect_s3_class(plotDeltaDistribution(results = combined,
        size = 3),
        "gtable")

    expect_s3_class(plotDeltaDistribution(results = combined,
        ncol = 10),
        "gtable")

    expect_s3_class(plotDeltaDistribution(results = combined,
        dots.on.top = FALSE),
        "gtable")

    expect_s3_class(plotDeltaDistribution(results = combined_prunedRef1, 
        this.color = "blue", pruned.color = "red"),
        "gtable")
})
