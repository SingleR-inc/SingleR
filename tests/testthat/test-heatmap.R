# Tests for visualization functions
# library(SingleR); library(testthat); source("setup.R"); source("test-heatmap.R")

colnames(test) <- sprintf("cell_%i", seq_len(ncol(test)))
pred <- SingleR(test=test, ref=training, labels=training$label, genes="de")
test.exp <- assay(test,1)
training.exp <- assay(training,1)

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
    expect_s3_class(plotScoreHeatmap(results = pred, show.labels = TRUE), "pheatmap")

    expect_s3_class(plotScoreHeatmap(results = pred, silent=TRUE), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred))), "pheatmap")
})

test_that("heatmap - 'cells.use' can be combined with annotations & annotations can be combined with eachother", {
    # default labels annot still there and displayed
    expect_s3_class(plotScoreHeatmap(results = pred,
        cells.use = 1:50), "pheatmap")

    # clusters and labels
    expect_s3_class(plotScoreHeatmap(
        results = pred, cells.use = 1:50, clusters = pred$labels), "pheatmap")

    # annot & clusters and labels
    expect_s3_class(plotScoreHeatmap(
        results = pred, cells.use = 1:50, clusters = pred$labels,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred))), "pheatmap")
})

test_that("heatmap - Error is thrown when order.by = `clusters` but no clusters are given.", {
    expect_error(plotScoreHeatmap(
        results = pred, cells.use = 1:50,
        order.by = "clusters"),
        "'clusters' input is required when 'order.by=\"clusters\"'")
})


test_that("heatmap - can pass excess pheatmap::pheatmap parameters through plotScoreHeatmap.", {
    expect_s3_class(plotScoreHeatmap(results = pred, cutree_col = 3), "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = pred, fontsize_row = 5), "pheatmap")
    expect_equal(plotScoreHeatmap(results = pred, silent = TRUE,
        fontsize_row = 5, return.data = TRUE)$fontsize_row,
        5)
})

test_that("heatmap scores color can be adjusted, regardless of 'normalize' value", {
    expect_equal(
        plotScoreHeatmap(results = pred, silent = TRUE, return.data = TRUE,
            normalize = FALSE,
            color = colorRampPalette(c("red", "blue"))(33))$color,
        colorRampPalette(c("red", "blue"))(33))
    expect_equal(
        plotScoreHeatmap(results = pred, silent = TRUE, return.data = TRUE,
            normalize = TRUE,
            color = colorRampPalette(c("red", "blue"))(33))$color,
        colorRampPalette(c("red", "blue"))(33))
})

test_that("heatmap allows users to adjust breaks, legend_breaks, legend_labels", {
    expect_s3_class(
        plotScoreHeatmap(results = pred, silent = TRUE,
            normalize = FALSE,
            color = colorRampPalette(c("red", "blue"))(33),
            breaks = seq(-5, 5, length.out = 34),
            legend_breaks = c(-5, 0, 5),
            legend_labels = c("manually", "set", "labels")),
        "pheatmap")
    non_norm_args <- plotScoreHeatmap(results = pred, silent = TRUE, return.data = TRUE,
        normalize = FALSE,
        color = colorRampPalette(c("red", "blue"))(33),
        breaks = seq(-5, 5, length.out = 34),
        legend_breaks = c(-5, 0, 5),
        legend_labels = c("manually", "set", "labels"))
    expect_equal(non_norm_args$breaks, seq(-5, 5, length.out = 34))
    expect_equal(non_norm_args$legend_breaks, c(-5, 0, 5))
    expect_equal(non_norm_args$legend_labels, c("manually", "set", "labels"))
    norm_args <- plotScoreHeatmap(results = pred, silent = TRUE, return.data = TRUE,
        normalize = TRUE,
        color = colorRampPalette(c("red", "blue"))(33),
        breaks = seq(-5, 5, length.out = 34),
        legend_breaks = c(-5, 0, 5),
        legend_labels = c("manually", "set", "labels"))
    expect_equal(norm_args$breaks, seq(-5, 5, length.out = 34))
    expect_equal(norm_args$legend_breaks, c(-5, 0, 5))
    expect_equal(norm_args$legend_labels, c("manually", "set", "labels"))
})

test_that("heatmap is adjusted properly when 'labels.use' yields 1 or 0 labels", {
    # Should give message but still output plot
    expect_warning(plotScoreHeatmap(results = pred,
        labels.use = c("A")),
        paste0("disabling normalization"))

    expect_equal(
        suppressMessages(plotScoreHeatmap(results = pred, silent = TRUE,
            labels.use = c("A"),
            color = colorRampPalette(c("red", "blue"))(33), # proximal to normalization being turned off
            return.data = TRUE)$color),
        colorRampPalette(c("red", "blue"))(33))

    # Should give message but still output plot
    expect_warning(plotScoreHeatmap(results = pred,
        labels.use = c("a")),
        paste0("ignoring 'labels.use'"))

    expect_equal(
        nrow(suppressMessages(plotScoreHeatmap(results = pred, silent = TRUE,
            labels.use = c("a"),
            return.data = TRUE)$mat)),
        5)
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
        # order.by = "clusters",
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
        order.by = "clusters",
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
        # order.by = "clusters",
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
        # order.by = "clusters",
        # cells.use = 1:50,
        clusters = seq(nrow(pred),1),
        show.pruned = TRUE,
        annotation_col = data.frame(
            annot = seq_len(nrow(pred)),
            row.names = row.names(pred)[seq(nrow(pred),1)])),
        "pheatmap")
})

test_that("Row and Column annotation coloring works", {
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

#######################################
### Prep for multi-reference checks ###
#######################################

ref <- .mockRefData(nreps=8)
ref1 <- ref[,1:4%%4==0]
ref1 <- ref1[,sample(ncol(ref1))]
ref1 <- scuttle::logNormCounts(ref1)

ref2 <- ref[,1:4%%4!=0]
ref2 <- ref2[,sample(ncol(ref2))]
ref2 <- scuttle::logNormCounts(ref2)

ref2$label <- tolower(ref2$label)

combined <- SingleR(
    test, ref = list(smallRef = ref1, largeRef = ref2),
    labels = list(ref1$label, ref2$label))

combined_prunedRef1 <- combined
combined_prunedRef1$orig.results$smallRef$pruned.labels[1:3%%3==0] <- NA_character_

ref1.pruned <- is.na(combined_prunedRef1$orig.results$smallRef$pruned.labels)
ref1.title <- "smallRef"
ref2.title <- "largeRef"

test_that("heatmap can be made for multi-ref runs - combined", {
    expect_s3_class(plotScoreHeatmap(results = combined, silent = TRUE,
        scores.use = 0),
        "pheatmap")
    # title correct
    expect_equal(plotScoreHeatmap(results = combined, silent = TRUE,
        scores.use = 0, return.data = TRUE)$main,
        "Combined Scores")
})

test_that("heatmap can be made for multi-ref runs - individual", {
    expect_s3_class(plotScoreHeatmap(results = combined, silent = TRUE,
        scores.use = 1),
        "pheatmap")
    # title correct
    expect_equal(plotScoreHeatmap(results = combined, silent = TRUE,
        scores.use = 1, return.data = TRUE)$main,
        paste(ref1.title,"Scores"))
})

test_that("heatmap can be made for multi-ref runs - multiple", {
    expect_s3_class(plotScoreHeatmap(results = combined,
        scores.use = 0:1),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        scores.use = NULL),
        "gtable")
    expect_equal(
        length(
            plotScoreHeatmap(results = combined, silent = TRUE, grid.vars = NULL,
                scores.use = NULL)),
        length(
            plotScoreHeatmap(results = combined, silent = TRUE, grid.vars = NULL,
                scores.use = 0:2))
        )
})

test_that("heatmap multi-ref - calls & pruned calls can be selected with calls.use", {
    # Individual
    expect_s3_class(plotScoreHeatmap(results = combined_prunedRef1, scores.use = 1,
        calls.use = 1, show.pruned = TRUE),
        "pheatmap")
    # Correct annotation title
    expect_true("smallRef Labels" %in%
        names(plotScoreHeatmap(results = combined_prunedRef1, scores.use = 1, silent = TRUE,
            calls.use = 1, show.pruned = TRUE, return.data = TRUE)$annotation_col))
    # Correct prune calls added
    expect_equal(sum(ref1.pruned),
        sum(plotScoreHeatmap(results = combined_prunedRef1, scores.use = 1, silent = TRUE,
            calls.use = 1, show.pruned = TRUE, return.data = TRUE)$annotation_col$Pruned==TRUE))

    # All
    expect_s3_class(plotScoreHeatmap(results = combined,
        calls.use = 1, show.pruned = TRUE,
        scores.use = NULL),
        "gtable")

    # Multiple calls.use
    expect_s3_class(plotScoreHeatmap(results = combined,
        calls.use = 0:2, show.pruned = TRUE,
        scores.use = NULL),
        "gtable")
})

test_that("heatmap multi-ref - grid.vars control", {
    expect_s3_class(plotScoreHeatmap(results = combined, scores.use = NULL,
        grid.vars = NULL)[[1]],
        "pheatmap")
    expect_s3_class(plotScoreHeatmap(results = combined, scores.use = NULL,
        grid.vars = list(ncol = 2)),
        "gtable")
})

test_that("heatmap multi-ref - 'na.color'", {
    expect_equal(
        tail(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE,
            scores.use = 0,
            na.color = "#000000")$color, 1),
        "#000000")
})

test_that("heatmap multi-ref - labels with no scores are removed", {
    combined$scores <- cbind(combined$scores, "f" = NA)
    expect_true("f" %in% colnames(combined$scores))
    expect_false("f" %in% rownames(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE,
            scores.use = 0)$mat))
})

test_that("heatmap multi-ref - labels with least calls/calcs are removed by 'max.labels'", {
    combined$scores <- cbind(combined$scores, "neverCalled" = 1) # actual score is immaterial
    combined$scores <- cbind(combined$scores, "rarelyCalc" = NA)
    combined$scores[1,"rarelyCalc"] <- 1 # Needs at least one score to not be removed anyway.
    expect_true(all(c("neverCalled", "rarelyCalc") %in% colnames(combined$scores)))

    # Both there with no trimming
    expect_true(all(c("neverCalled", "rarelyCalc") %in% rownames(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE, scores.use = 0,
        max.labels = 40)$mat)))

    # The rarely picked for calculation "rarelyCalc" label should be removed first
    expect_true("neverCalled" %in% rownames(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE, scores.use = 0,
        max.labels = 11)$mat))
    expect_false("rarelyCalc" %in% rownames(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE, scores.use = 0,
        max.labels = 11)$mat))

    # The never picked as final label "neverCalled" label should be removed next
    expect_false("neverCalled" %in% rownames(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE, scores.use = 0,
        max.labels = 10)$mat))
    expect_false("rarelyCalc" %in% rownames(plotScoreHeatmap(results = combined, silent = TRUE, return.data = TRUE, scores.use = 0,
        max.labels = 10)$mat))
})

test_that("heatmap multi-ref - Other typical adjustments throw no unexpected errors", {
    # Our vars
    expect_s3_class(plotScoreHeatmap(results = combined,
        normalize = FALSE),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        labels.use = c("A", "a")),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        max.labels = 3),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        clusters = g),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        order.by = "clusters", clusters = g),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        cluster_col = TRUE),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        cells.order = seq_len(nrow(combined))),
        "gtable")
    expect_s3_class(plotScoreHeatmap(results = combined,
        cells.use = 1:20),
        "gtable")

    # pheatmap var
    expect_s3_class(plotScoreHeatmap(results = combined,
        treeheight_row = 5),
        "gtable")
})

test_that("heatmap - max.labels trim when duplicate labels", {
    combined_dup1 <- SingleR(
        test, ref = list(smallRef = ref1, smallRef2 = ref1, largeRef = ref2),
        labels = list(ref1$label, ref1$label, ref2$label))
    expect_s3_class(plotScoreHeatmap(results = pred, max.labels = 10), "pheatmap")
})

test_that("heatmap - rows.order sets row order of heatmap and warns on missing values", {
    # Combined
    expect_s3_class(plotScoreHeatmap(results = combined,
        rows.order = c("A","a","B","b","C","c","D","d","E","e")),
        "gtable")
    # Single (extra labels okay)
    expect_s3_class(plotScoreHeatmap(results = pred,
        rows.order = c("A","a","B","b","C","c","D","d","E","e")),
        "pheatmap")
    # Warn on missing labels
    expect_warning(plotScoreHeatmap(results = pred,
        rows.order = c("A","a","B","b","C","c")),
        "Label(s) of Scores missing from 'rows.order' will not be plotted: D, E", fixed = TRUE)
})
