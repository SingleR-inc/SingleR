# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-visualizations.R")

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
  
    expect_s3_class(plotScoreHeatmap(results = pred, cells.use = 1:50), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, cells.use = rownames(pred)[1:50]), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, cells.order=seq_len(nrow(pred))), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, labels.use = levels(as.factor(pred$labels))[1:3]), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, max.labels = length(levels(as.factor(pred$labels)))-1), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, clusters = pred$labels), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, clusters = pred$labels, order.by.clusters=TRUE), "pheatmap")
    
    expect_s3_class(plotScoreHeatmap(results = pred, show.pruned = TRUE), "pheatmap")
  
    expect_s3_class(plotScoreHeatmap(results = pred, silent=TRUE), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred))), "pheatmap")
})

test_that("cells.use can be combined with annotations & annotations can be combined with eachother", {
    expect_s3_class(plotScoreHeatmap(results = pred, cells.use = 1:50, clusters = pred$labels), "pheatmap")

    expect_s3_class(plotScoreHeatmap(
        results = pred, cells.use = 1:50, clusters = pred$labels,
        show.pruned = TRUE), "pheatmap")

    expect_s3_class(plotScoreHeatmap(
        results = pred, cells.use = 1:50, clusters = pred$labels,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred)),
        show.pruned = TRUE), "pheatmap")
})

test_that("cells.use can be combined with ordering (by cells or by cluster)", {
    expect_s3_class(plotScoreHeatmap(
        results = pred, cells.use = 1:50, clusters = pred$labels,
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred)),
        cells.order = 1:100), "pheatmap")

    expect_s3_class(plotScoreHeatmap(
        results = pred, cells.use = 1:50, clusters = pred$labels,
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred)),
        order.by.clusters = TRUE), "pheatmap")
})

test_that("We can pass excess pheatmap::pheatmap parameters through plotScoreHeatmap.", {
    expect_s3_class(plotScoreHeatmap(results = pred, cutree_col = 3), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, fontsize.row = 5), "pheatmap")
})

####################################
#### Manual Visualization Check ####
####################################

test_that("Annotations stay linked, even with cells.use, cells.order, or order.by.clusters = TRUE", {
    # Make prune.call TRUE for every 10th value.  (We need known order for testing annotation placement.)
    pred$pruned.labels <- rep(c(rep(FALSE,9),NA),nrow(pred)/10)
    
    #Reference plot: Every tenth cell, pruned = TRUE. Clusters from 100:1. annot from 1:100.
    expect_s3_class(plotScoreHeatmap(
        results = pred,
        cells.order = seq_len(nrow(pred)),
        # order.by.clusters = TRUE,
        # cells.use = 1:50,
        clusters = seq(nrow(pred),1),
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred))),
        "pheatmap")

    #Reversed order: First, 11th, 21st... cell, pruned = TRUE. Clusters from 1:100. annot from 100:1.
    expect_s3_class(plotScoreHeatmap(
        results = pred,
        # cells.order = seq_len(nrow(pred)),
        order.by.clusters = TRUE,
        # cells.use = 1:50,
        clusters = seq(nrow(pred),1),
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred))),
        "pheatmap")

    #Reference plot, but only half: Every tenth cell, pruned = TRUE. Clusters from 50:1. annot from 100:51.
    expect_s3_class(plotScoreHeatmap(
        results = pred,
        cells.order = seq_len(nrow(pred)),
        # order.by.clusters = TRUE,
        cells.use = 1:50,
        clusters = seq(nrow(pred),1),
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred))),
        "pheatmap")

    #Reference plot, but with annot flipped 100:1 because it's rownames were flipped.
    expect_s3_class(plotScoreHeatmap(
        results = pred,
        cells.order = seq_len(nrow(pred)),
        # order.by.clusters = TRUE,
        # cells.use = 1:50,
        clusters = seq(nrow(pred),1),
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred)[seq(nrow(pred),1)])),
        "pheatmap")
})

test_that("Row and Column annotation coloring works", {
    # Make prune.call TRUE for every 10th value.
    pred$pruned.labels <- rep(c(rep(FALSE,9),NA),nrow(pred)/10)
    
    #When works:
        # Clusters and Continuous are shades of the same color
        # Pruned and Discrete are many discrete colors
    expect_s3_class(plotScoreHeatmap(
        results = pred,
        cells.order = seq_len(nrow(pred)),
        clusters = seq(nrow(pred),1),
        show.pruned = TRUE,
        annotation_row = data.frame(
            Discrete = as.character(seq_len(ncol(pred$scores))),
            Continuous = as.numeric(seq_len(ncol(pred$scores))),
            row.names = colnames(pred$scores))),
        "pheatmap")
})
