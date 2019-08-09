# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-visualization.R")

test.exp <- assay(test,1)
training.exp <- assay(training,1)

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

example("SingleR")
colnames(test) <- paste0("name",seq_len(ncol(counts(test))))
pred <- SingleR(test, sce, labels=sce$label)

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
})

test_that("Explicit visualization adjustments don't break plotScoreHeatmap plotting.",{
  expect_s3_class(plotScoreHeatmap(results = pred,
                                   fontsize.row = 5),
                  "pheatmap")
})

test_that("We can pass excess pheatmap::pheatmap parameters through plotScoreHeatmap.", {
  plotScoreHeatmap(results = pred,
                   cutree_col = 3)
})

# I'm not sure how best to compare two plots (ggplot nor pheatmap).
# test_that("Providing indices versus cell names in plotScoreHeatmap 'cells.use' does the same thing",{
#   hm1 <- plotScoreHeatmap(results = pred,
#                           cells.use = 1:50)
#   hm2 <- plotScoreHeatmap(results = pred,
#                           cells.use = rownames(pred)[1:50])
#   # These outputs of these two plots, made above, should be the same, but some names underneath the hood are different.
#   # expect_equal(hm1, hm2) #This fails because it checks unplotted info.
# })
# test_that("When order.by.clusters = TRUE, cells.order is overrided in plotScoreHeatmap 'cells.use'",{
#   hm1 <- plotScoreHeatmap(results = pred,
#                           clusters = pred$labels,
#                           order.by.clusters=TRUE)
#   hm2 <- plotScoreHeatmap(results = pred,
#                           clusters = pred$labels,
#                           cells.order=seq_len(nrow(pred)), #not utilized because order.by.clusters takes precedence
#                           order.by.clusters=TRUE)
#   hm3 <- plotScoreHeatmap(results = pred,
#                           clusters = pred$labels,
#                           cells.order=seq_len(nrow(pred)),
#                           order.by.clusters=FALSE)
#   # I'm not sure how to do this! The info inside a pheatmap is hard to parse and probably not great to access anyway.
#   # hm1 and hm2 shoul dbe the same plot
#   # hm3 should be a different plot
# })
