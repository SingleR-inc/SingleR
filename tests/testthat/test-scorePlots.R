# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-scorePlots.R")

colnames(test) <- sprintf("cell_%i", seq_len(ncol(test)))
pred <- SingleR(test=test, ref=training, labels=training$label, genes="de", prune = FALSE)
pred2 <- SingleR(test=test, ref=training, labels=training$label, genes="de")

test_that("we can produce single cell scorePlots without pruneScores run", {
    expect_s3_class(
        plotScoresSingleCell(results = pred, cell.id = 1),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred, cell.id = 1, labels.use = c("A")),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred, cell.id = 1, labels.use = c("A","B")),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred, cell.id = 1, size = 10),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred, cell.id = 1,
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_error(
        plotScoresSingleCell(
            results = pred, cell.id = 1,
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce single cell scorePlots with pruned scores", {
    expect_s3_class(
        plotScoresSingleCell(results = pred2, cell.id = 1),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred2, cell.id = 1, labels.use = c("A")),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred2, cell.id = 1, labels.use = c("A","B")),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred2, cell.id = 1, size = 10),
        "ggplot")
    expect_s3_class(
        plotScoresSingleCell(
            results = pred2, cell.id = 1,
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_error(
        plotScoresSingleCell(
            results = pred2, cell.id = 1,
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce single label scorePlots without pruneScores run", {
    expect_s3_class(
        plotScoresSingleLabel(results = pred, label = "A"),
        "ggplot")
    expect_s3_class(
        plotScoresSingleLabel(
            results = pred, label = "A", size = 10),
        "ggplot")
    expect_s3_class(
        plotScoresSingleLabel(
            results = pred, label = "A", dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoresSingleLabel(
            results = pred, label = "A",
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_error(
        plotScoresSingleLabel(
            results = pred, label = "A",
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce single label scorePlots with pruned scores", {
    expect_s3_class(
        plotScoresSingleLabel(results = pred2, label = "A"),
        "ggplot")
    expect_s3_class(
        plotScoresSingleLabel(
            results = pred2, label = "A", size = 10),
        "ggplot")
    expect_s3_class(
        plotScoresSingleLabel(
            results = pred2, label = "A", dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoresSingleLabel(
            results = pred2, label = "A",
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_error(
        plotScoresSingleLabel(
            results = pred2, label = "A",
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce multi label scorePlots without pruneScores run", {
    expect_s3_class(
        plotScoresMultiLabels(results = pred),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred, size = 10),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred, dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred, labels.use = c("A","B","D")),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred, ncol = 3),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(
            results = pred,
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_error(
        plotScoresMultiLabels(
            results = pred,
            colors = c("blue", "yellow")),
        NULL)
})

test_that("we can produce multi label scorePlots with pruned scores", {
    expect_s3_class(
        plotScoresMultiLabels(results = pred2),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred2, size = 10),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred2, dots.on.top = TRUE),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred2, labels.use = c("A","B","D")),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(results = pred2, ncol = 3),
        "ggplot")
    expect_s3_class(
        plotScoresMultiLabels(
            results = pred2,
            colors = c("blue", "yellow", "orange")),
        "ggplot")
    expect_error(
        plotScoresMultiLabels(
            results = pred2,
            colors = c("blue", "yellow")),
        NULL)
})
