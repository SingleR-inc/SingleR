# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-scoreDistributions.R")

colnames(test) <- sprintf("cell_%i", seq_len(ncol(test)))
pred1 <- SingleR(test=test, ref=training, labels=training$label, genes="de")
pred2 <- pred1
pred1$pruned.labels <- pred1$labels
pred2$pruned.labels[seq(5,80,5)]<-NA
pred3 <- SingleR(test=test, ref=training, labels=training$label, genes="de", fine.tune = FALSE)

# pred1 = no pruned.labels
# pred2 = 16 pruned labels
# pred3 = no tuning run

test_that("we can produce single label scoreDistributions when no labels were pruned", {
    expect_s3_class(
        plotScoreDistribution(results = pred1, labels = "A"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels = "A", size = 10),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels = "A", dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels = "A",
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels = "A",
            size = 5),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels = "A",
            show = "delta.med"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, labels = "A",
            show = "delta.next"),
        "ggplot")
    expect_error(
        plotScoreDistribution(
            results = pred1, labels = "A",
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce multi label scoreDistributions when no labels were pruned", {
    expect_s3_class(
        plotScoreDistribution(results = pred1),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred1, size = 10),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred1, dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred1, labels = c("A","B","D")),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred1, ncol = 3),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1,
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, size = 5),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, show = "delta.med"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred1, show = "delta.next"),
        "ggplot")
    expect_error(
        plotScoreDistribution(
            results = pred1,
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce single label scoreDistributions", {
    expect_s3_class(
        plotScoreDistribution(results = pred2, labels = "A"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels = "A", size = 10),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels = "A", dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels = "A",
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels = "A",
            size = 5),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels = "A",
            show = "delta.med"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, labels = "A",
            show = "delta.next"),
        "ggplot")
    expect_error(
        plotScoreDistribution(
            results = pred2, labels = "A",
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce multi label scoreDistributions", {
    expect_s3_class(
        plotScoreDistribution(results = pred2),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred2, size = 10),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred2, dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred2, labels = c("A","B","D")),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred2, ncol = 3),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2,
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, size = 5),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, show = "delta.med"),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred2, show = "delta.next"),
        "ggplot")
    expect_error(
        plotScoreDistribution(
            results = pred2,
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce multi label scoreDistributions when tuning was not run.", {
    expect_s3_class(
        plotScoreDistribution(results = pred3),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(results = pred3, show = "delta.med"),
        "ggplot")
    expect_error(
        plotScoreDistribution(results = pred3, show = "delta.next"),
        NA)
})

test_that("we can add cutoffs to single and multi-label plots.", {
    expect_s3_class(
        plotScoreDistribution(
            results = pred3, show = "delta.med", show.nmads = 3),
        "ggplot")
    expect_s3_class(
        plotScoreDistribution(
            results = pred3, show = "delta.med", show.min.diff = 0.05),
        "ggplot")
    expect_error(
        plotScoreDistribution(
            results = pred3, show = "delta.next", show.min.diff = 0.05),
        NA)
})

