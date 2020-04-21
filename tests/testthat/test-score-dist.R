# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-score-dist.R")

# pred1 = no pruned.labels
pred1 <- SingleR(test=test, ref=training, labels=training$label, genes="de")

# pred2 = 16 pruned labels
pred2 <- pred1
pred2$pruned.labels[seq(5,80,5)]<-NA

# pred3 = no tuning run
pred3 <- SingleR(test=test, ref=training, labels=training$label, genes="de", fine.tune = FALSE)

test_that("plot-dist - no pruned, one label", {
    expect_s3_class(
        plotScoreDistribution(results = pred1, labels.use = "A"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A", size = 10),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A", dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A",
            this.color = "blue",
            pruned.color = "yellow",
            other.color = "orange"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A",
            size = 5),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A",
            show = "scores"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels.use = "A",
            show = "delta.next"),
        "ggplot")
})

test_that("plot-dist - no pruned, many labels", {
    expect_s3_class(
        plotScoreDistribution(results = pred1),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1,
            this.color = "blue",
            pruned.color = "yellow",
            other.color = "orange"),
        "ggplot")
})

test_that("plot-dist - with pruned, single label", {
    expect_s3_class(
        plotScoreDistribution(results = pred2, labels.use = "A"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A", size = 10),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A", dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A",
            this.color = "blue",
            pruned.color = "yellow",
            other.color = "orange"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A",
            size = 5),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A",
            show = "scores"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels.use = "A",
            show = "delta.next"),
        "ggplot")
})

test_that("score-dist, with pruned, many labels", {
    expect_s3_class(
        plotScoreDistribution(results = pred2),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2,
            this.color = "blue",
            pruned.color = "yellow",
            other.color = "orange"),
        "ggplot")
})

test_that("plot-scores, no tuning run", {
    expect_s3_class(
        plotScoreDistribution(results = pred3),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred3, show = "delta.med"),
        "ggplot")
    expect_error(
        plotScoreDistribution(results = pred3, show = "delta.next"),
        "lacks fine-tuning diagnostics")
})

test_that("plot-scores - cutoff lines can be added/adjusted", {
    # nmads removed or moved
    expect_s3_class(
        plotScoreDistribution(
            results = pred3, show = "delta.med", show.nmads = NULL),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred3, show = "delta.med", show.nmads = 1),
        "ggplot")

    # min.diff(.med)
    expect_s3_class(
        plotScoreDistribution(
            results = pred3, show = "delta.med", show.min.diff = 0.05),
        "ggplot")

    # min.diff(.next)
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, show = "delta.next", show.min.diff = 0.05),
        "ggplot")
})

test_that("plot-dist ignores 'labels.use' if it yields 0 labels", {
    # Should give message but still output plot
    expect_message(plotScoreDistribution(results = pred1,
        labels.use = c("a")),
        paste0("No 'labels.use' in scores. Ignoring input."))
})

#######################################
### Prep for multi-reference checks ###
#######################################

ref <- .mockRefData(nreps=8)
ref1 <- ref[,1:4%%4==0]
ref1 <- ref1[,sample(ncol(ref1))]
ref1 <- scater::logNormCounts(ref1)

ref2 <- ref[,1:4%%4!=0]
ref2 <- ref2[,sample(ncol(ref2))]
ref2 <- scater::logNormCounts(ref2)

ref2$label <- tolower(ref2$label)

combined <- SingleR(
    test, ref = list(smallRef = ref1, largeRef = ref2),
    labels = list(ref1$label, ref2$label))

combined_prunedRef1 <- combined
combined_prunedRef1$orig.results$smallRef$pruned.labels[1:3%%3==0] <- NA_character_

ref1.pruned <- is.na(combined_prunedRef1$orig.results$smallRef$pruned.labels)
ref1.title <- "smallRef"
ref2.title <- "largeRef"

test_that("dist-plot can be made for multi-ref runs - combined", {
    expect_s3_class(plotScoreDistribution(results = combined,
        scores.use = 0, show = "scores"),
        "ggplot")
    expect_error(plotScoreDistribution(results = combined,
        scores.use = 0, show = "delta.med"),
        "'show = delta.med' cannot be used for combined results.")
    expect_error(plotScoreDistribution(results = combined,
        scores.use = 0, show = "delta.next"),
        "'show = delta.next' cannot be used for combined results.")
})

test_that("dist-plot can be made for multi-ref runs - individual", {
    expect_s3_class(plotScoreDistribution(results = combined,
        scores.use = 1, show = "scores"),
        "ggplot")
    expect_s3_class(plotScoreDistribution(results = combined,
        scores.use = 1, show = "delta.med"),
        "ggplot")
    expect_s3_class(plotScoreDistribution(results = combined,
        scores.use = 1, show = "delta.next"),
        "ggplot")
})

test_that("dist-plot can be made for multi-ref runs - multiple", {
    expect_s3_class(plotScoreDistribution(results = combined,
        scores.use = 0:1),
        "gtable")
    expect_s3_class(plotScoreDistribution(results = combined,
        scores.use = NULL),
        "gtable")

    # when scores.use = NULL combined plot should not be requested unless 'show = scores'
    expect_equal(
        length(
            plotScoreDistribution(results = combined, grid.vars = NULL,
                scores.use = NULL)),
        length(
            plotScoreDistribution(results = combined, grid.vars = NULL,
                scores.use = 0:2))
        )
    expect_equal(
        length(
            plotScoreDistribution(results = combined, grid.vars = NULL,
                show = "delta.med", scores.use = NULL)),
        length(
            plotScoreDistribution(results = combined, grid.vars = NULL,
                show = "delta.med", scores.use = 1:2))
        )
})

test_that("plot-dist multi-ref - calls & pruned calls can be selected with calls.use & pruned.use", {
    # Individual scores.use
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1, scores.use = 1,
        calls.use = 1, pruned.use = 1),
        "ggplot")
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1, scores.use = 1,
        calls.use = 1, pruned.use = 1,
        show = "delta.med"),
        "ggplot")
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1, scores.use = 1,
        calls.use = 1, pruned.use = 1,
        show = "delta.next"),
        "ggplot")

    # All scores.use
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1,
        calls.use = 1, pruned.use = 1, scores.use = NULL),
        "gtable")
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1,
        calls.use = 1, pruned.use = 1, scores.use = NULL,
        show = "delta.med"),
        "gtable")
    # Outputs but with message
    expect_s3_class(suppressMessages(plotScoreDistribution(results = combined_prunedRef1,
        calls.use = 1, pruned.use = 1, scores.use = NULL,
        show = "delta.next")),
        "gtable")
    expect_message(plotScoreDistribution(results = combined_prunedRef1,
        calls.use = 1, pruned.use = 1, scores.use = NULL,
        show = "delta.next"),
        "'calls.use', and 'scores.use' must be the same for 'show=\"delta.next\"'. 'calls.use' updated.")

    # Multiple calls.use
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1,
        calls.use = 0:2, pruned.use = 1, scores.use = NULL),
        "gtable")

    # Multiple pruned.use
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1,
        calls.use = 1, pruned.use = 0:2, scores.use = NULL),
        "gtable")
})

test_that("plot-dist multi-ref - grid.vars control", {
    expect_s3_class(plotScoreDistribution(results = combined, scores.use = NULL,
        grid.vars = NULL)[[1]],
        "ggplot")
    expect_s3_class(plotScoreDistribution(results = combined, scores.use = NULL,
        grid.vars = list(ncol = 2)),
        "gtable")
})

test_that("plot-dist multi-ref - Other typical adjustments throw no unexpected errors", {
    expect_s3_class(plotScoreDistribution(results = combined,
        labels.use = c("A", "a")),
        "gtable")
    expect_s3_class(plotScoreDistribution(results = combined,
        size = 3),
        "gtable")
    expect_s3_class(plotScoreDistribution(results = combined,
        ncol = 10),
        "gtable")
    expect_s3_class(plotScoreDistribution(results = combined,
        dots.on.top = FALSE),
        "gtable")
    expect_s3_class(plotScoreDistribution(results = combined_prunedRef1, pruned.use = 1,
        this.color = "blue", pruned.color = "red", other.color = "white"),
        "gtable")
})
