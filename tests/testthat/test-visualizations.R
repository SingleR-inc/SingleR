# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-visualization.R")

colnames(test) <- paste0("name",seq_len(ncol(counts(test))))
test.exp <- assay(test,1)
training.exp <- assay(training,1)
pred <- SingleR(test=test.exp, training=training.exp, labels=training$label, genes="de")

test_that("we can produce expression comparison scatterplots with plotCellVsReference", {
  expect_error(plotCellVsReference(test,1,training,1), NULL) # looks for logcounts, but the dataset does not have logcounts
  expect_s3_class(P <- plotCellVsReference(test, 1, training, 1, 1, 1), "ggplot")
  expect_s3_class(plotCellVsReference(test.exp, 1, training.exp, 1), "ggplot")
  expect_s3_class(plotCellVsReference(test, 1, training.exp, 1,1,NULL), "ggplot")
  expect_s3_class(plotCellVsReference(test.exp, 1, training, 1,NULL,1), "ggplot")
  expect_false(isTRUE(all.equal(P,
                                P2 <- plotCellVsReference(test, 2, training, 1, 1, 1))))
  expect_false(isTRUE(all.equal(P, 
                                P3 <- plotCellVsReference(test, 1, training, 2, 1, 1))))
  expect_false(isTRUE(all.equal(P2, P3)))
})

test_that("We can produce heatmaps of scores with plotScoreHeatmap", {
  expect_s3_class(plotScoreHeatmap(results = pred),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                          cells.use = 1:50),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                          cells.use = rownames(pred)[1:50]),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   labels.use = levels(as.factor(pred$labels))[1:3]),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   clusters = pred$labels),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   max.labels = length(levels(as.factor(pred$labels)))-1),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   normalize = FALSE),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   cells.order=seq_len(nrow(pred))),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   clusters = pred$labels,
                                   order.by.clusters=TRUE),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   silent=TRUE),
                  "pheatmap")
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   fontsize.row = 5),
                  "pheatmap")
})

test_that("We can pass excess pheatmap::pheatmap parameters through plotScoreHeatmap.", {
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   cutree_col = 3),
                  "pheatmap")
})
