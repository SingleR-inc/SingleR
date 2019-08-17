# Tests for visualization functions
# library(SingleR); library(testthat); source("../GitHub/SingleR/tests/testthat/setup.R"); source("../GitHub/SingleR/tests/testthat/test-visualizations.R")

colnames(test) <- sprintf("cell_%i", seq_len(ncol(test)))
pred <- SingleR(test=test, ref=training, labels=training$label, genes="de")
test.exp <- assay(test,1)
training.exp <- assay(training,1)

test_that("we can produce expression comparison scatterplots with plotCellVsReference", {
    expect_error(plotCellVsReference(
        test = test, test.id = 1,
        ref = training, ref.id = 1,
        assay.type.test = 5, assay.type.ref = 5), NULL)
    expect_s3_class(plotCellVsReference(
        test = test, test.id = 1,
        ref = training, ref.id = 1,
        assay.type.test = 1, assay.type.ref = 1), "ggplot")
    expect_s3_class(P <- plotCellVsReference(
        test = test, test.id = 1,
        ref = training, ref.id = 1), "ggplot")
    expect_s3_class(plotCellVsReference(
        test = test.exp, test.id = 1,
        ref = training.exp, ref.id = 1), "ggplot")
    expect_s3_class(plotCellVsReference(
        test = test, test.id = 1,
        ref = training.exp, ref.id = 1,
        assay.type.test = 1, assay.type.ref = 5), "ggplot")
    expect_s3_class(plotCellVsReference(
        test = test.exp, test.id = 1,
        ref = training, ref.id = 1,
        assay.type.test = 5, assay.type.ref = 1), "ggplot")
    expect_false(isTRUE(all.equal(
        P,
        P2 <- plotCellVsReference(
            test = test, test.id = 2,
            ref = training, ref.id = 1))))
    expect_false(isTRUE(all.equal(
        P, 
        P3 <- plotCellVsReference(
            test = test, test.id = 1,
            ref = training, ref.id = 2))))
    expect_false(isTRUE(all.equal(P2, P3)))
})

test_that("We can produce heatmaps of scores with plotScoreHeatmap", {
    expect_s3_class(plotScoreHeatmap(results = pred), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, normalize = FALSE), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, cube.normalize = FALSE), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, cells.use = 1:50), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, cells.use = rownames(pred)[1:50]), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, cells.order=seq_len(nrow(pred))), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, labels.use = levels(as.factor(pred$labels))[1:3]), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, max.labels = length(levels(as.factor(pred$labels)))-1), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, clusters = pred$labels), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, clusters = pred$labels, order.by.clusters=TRUE), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, silent=TRUE), "pheatmap")
})

test_that("We can pass excess pheatmap::pheatmap parameters through plotScoreHeatmap.", {
    expect_s3_class(plotScoreHeatmap(results = pred, cutree_col = 3), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, fontsize.row = 5), "pheatmap")
})
