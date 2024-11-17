# library(SingleR); library(testthat); source("setup.R"); source("test-plotMarkerHeatmap.R")

colnames(test) <- sprintf("cell_%i", seq_len(ncol(test)))
pred <- SingleR(test=test, ref=training, labels=training$label, genes="de")

test_that("We can produce heatmaps of markers", {
    expect_s3_class(plotMarkerHeatmap(pred, test, pred$labels[1]), "pheatmap")
})

test_that("marker heatmaps handle pruned NAs correctly", {
    pred$pruned.labels[1] <- NA
    lab <- na.omit(unique(pred$pruned.labels))[1]
    out <- plotMarkerHeatmap(pred, test, lab, use.pruned=TRUE)
    expect_s3_class(out, "pheatmap")
})

test_that("marker heatmaps work with a subset of other labels", {
    expect_s3_class(plotMarkerHeatmap(pred, test, pred$labels[1], other.labels=unique(pred$labels)[1:2]), "pheatmap")
})

test_that("marker heatmaps work with other names", {
    alt.names <- paste0("X_", rownames(test))
    expect_s3_class(plotMarkerHeatmap(pred, test, pred$labels[1], display.row.names=alt.names), "pheatmap")
})

test_that("marker heatmap falls back to average abundances", {
    lab <- pred$labels[1]
    keep <- pred$labels == lab
    pred <- pred[keep,]
    test <- test[,keep,drop=FALSE]
    expect_s3_class(plotMarkerHeatmap(pred, test, lab), "pheatmap")
})
