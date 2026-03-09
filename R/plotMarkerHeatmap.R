#' Plot a heatmap of the markers for a label
#'
#' Create a heatmap of the log-normalized expression for the most interesting markers of a particular label.
#'
#' @param results A \link[S4Vectors]{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}} on a single reference.
#' @param test A numeric matrix of log-normalized expression values where rows are genes and columns are cells.
#' Each row should be named with the same gene name that was used to compute \code{results}.
#'
#' Alternatively, a \link[SummarizedExperiment]{SummarizedExperiment} object containing such a matrix.
#' @param label String specifying the label of interest.
#' @param other.labels Character vector specifying the other labels to be compared to \code{label} when finding interesting markers.
#' Defaults to all available labels.
#' @param assay.type Integer scalar or string specifying the matrix of expression values to use if \code{test} is a \link[SummarizedExperiment]{SummarizedExperiment}.
#' @param use.pruned Logical scalar indicating whether the pruned labels should be used instead.
#' @param order.by.effect String specifying the effect size from \code{\link[scrapper]{scoreMarkers}} with which to sort for interesting markers.
#' @param order.by.summary String specifying the summary statistic from \code{\link[scrapper]{scoreMarkers}} with which to sort for interesting markers.
#' @param score.args Named list of additional arguments to pass to \code{\link[scrapper]{scoreMarkers}}, e.g., \code{block}, \code{threshold}.
#' @param display.row.names Character vector of length equal to the number of rows of \code{test},
#' containing the names of the features to show on the heatmap (e.g., to replace IDs with symbols).
#' If \code{NULL}, the existing row names of \code{test} are used.
#' @param top Integer scalar indicating the top most interesting markers to show in the heatmap.
#' @param average Boolean indicating whether to visualize the average expression profile for each label.
#' If \code{FALSE}, the expression profile for each cell in \code{test} is shown.
#' @param center Boolean indicating whether to center the expression profiles for each gene at zero.
#' @param num.threads Integer scalar specifying the number to threads to use.
#' @param BPPARAM Deprecated, use \code{num.threads} instead.
#' @param ... Additional parameters for heatmap control passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return 
#' For \code{plotMarkerHeatmap}, the output of \code{\link[pheatmap]{pheatmap}} is returned showing the heatmap on the current graphics device.
#'
#' For \code{configureMarkerHeatmap}, a list is returned containing:
#' \itemize{
#' \item \code{rows}, an integer vector of row indices for the markers of \code{label}, ordered from most to least interesting.
#' \item \code{columns}, an integer vector of column indices to show in the heatmap.
#' This is ordered by the predicted labels so that cells assigned to the same label are contiguous.
#' \item \code{predictions}, a character vector of predicted labels for cells to be shown in the heatmap.
#' Each entry corresponds to an entry of \code{columns}.
#' The labels in this vector are guaranteed to be sorted.
#' }
#'
#' @details
#' The \code{plotMarkerHeatmap} function creates a heatmap where each row is a marker gene for \code{label} and each column is a cell in the test dataset.
#' The aim is to check the effectiveness of the reference-derived markers for distinguishing between labels in the test dataset.
#' \dQuote{Interesting} markers should show strong upregulation in cells assigned to \code{label} compared to cells assigned to all \code{other.labels}. 
#' We identify such markers by scoring all reference-derived markers with \code{\link[scrapper]{scoreMarkers}} on the \code{test} expression.
#' The \code{top} markers based on the specified \code{order.by.*} fields are shown in the heatmap.
#' If only one label is present, markers are ranked by average abundance intead. 
#'
#' The \code{configureMarkerHeatmap} function performs all the calculations underlying \code{plotMarkerHeatmap}.
#' This can be used to construct a similar heatmap in other plotting frameworks, e.g., \pkg{scater}, \pkg{dittoSeq}.
#'
#' @author Aaron Lun
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' plotMarkerHeatmap(pred, test, pred$labels[1])
#' plotMarkerHeatmap(pred, test, pred$labels[1], use.pruned=TRUE)
#' plotMarkerHeatmap(pred, test, pred$labels[1], order.by.effect="auc")
#' plotMarkerHeatmap(pred, test, pred$labels[1], center=TRUE)
#'
#' # Manually configuring a simpler heatmap by label:
#' config <- configureMarkerHeatmap(pred, test, pred$labels[1])
#' mat <- assay(test, "logcounts")[head(config$rows, 20), config$columns]
#' aggregated <- scrapper::aggregateAcrossCells(mat, list(config$predictions))
#' averages <- t(t(aggregated$sums) / aggregated$counts)
#' colnames(averages) <- aggregated$combinations[,1]
#' pheatmap::pheatmap(averages, cluster_col=FALSE)
#'
#' # ... which is basically the same as this:
#' plotMarkerHeatmap(pred, test, pred$labels[1], average=TRUE)
#' 
#' @export
#' @importFrom Matrix rowMeans
#' @importFrom utils head
plotMarkerHeatmap <- function(
    results,
    test,
    label,
    other.labels = NULL,
    assay.type = "logcounts",
    display.row.names = NULL,
    use.pruned = FALSE,
    order.by.effect = "cohens.d",
    order.by.summary = "min.rank",
    score.args = list(),
    average = FALSE,
    center = FALSE,
    top = 20,
    num.threads = 1,
    BPPARAM = NULL,
    ... 
) {
    num.threads <- .get_num_threads(num.threads, BPPARAM)
    test <- .to_clean_matrix(test, assay.type, check.missing=FALSE, num.threads=num.threads)
    config <- configureMarkerHeatmap(
        results, 
        test,
        label=label,
        other.labels=other.labels,
        assay.type=assay.type,
        use.pruned=use.pruned,
        order.by.effect=order.by.effect,
        order.by.summary=order.by.summary,
        score.args=score.args,
        num.threads=num.threads
    )

    to.show <- head(config$rows, top)
    predictions <- config$predictions
    test <- test[to.show, config$columns, drop=FALSE]
    if (!is.null(display.row.names)) {
        rownames(test) <- display.row.names[to.show]
    }

    if (average) {
        aggregated <- scrapper::aggregateAcrossCells(test, list(predictions), num.threads=num.threads)
        test <- t(t(aggregated$sums) / aggregated$counts)
        predictions <- aggregated$combinations[,1]
    }

    N <- 25
    if (center) {
        test <- test - rowMeans(test)
        color <- grDevices::colorRampPalette(c("blue", "white", "red"))(N)
        maxed <- max(abs(test)) 
        breaks <- seq(-maxed, maxed, length.out=N+1)
    } else {
        color <- viridis::viridis(N)
        limits <- range(test, na.rm=TRUE)
        breaks <- seq(limits[1], limits[2], length.out=N+1)
    }

    colnames(test) <- seq_len(ncol(test))
    pheatmap::pheatmap(
        test,
        breaks=breaks,
        color=color,
        annotation_col=data.frame(labels=predictions, row.names=colnames(test)),
        cluster_col=FALSE,
        show_colnames=FALSE,
        ...
    )
}

#' @export
#' @rdname plotMarkerHeatmap
#' @importFrom S4Vectors metadata
#' @importFrom Matrix rowMeans
configureMarkerHeatmap <- function(
    results,
    test,
    label,
    other.labels = NULL,
    assay.type = "logcounts",
    use.pruned = FALSE,
    order.by.effect = "cohens.d",
    order.by.summary = "min.rank",
    score.args = list(),
    num.threads=1
) {
    test <- .to_clean_matrix(test, assay.type, check.missing=FALSE, num.threads=num.threads)
    all.markers <- metadata(results)$de.genes[[label]]

    if (use.pruned) {
        labfield <- "pruned.labels"
    } else {
        labfield <- "labels"
    }
    predictions <- results[[labfield]]

    ckeep <- seq_len(ncol(test))
    if (!is.null(other.labels)) {
        ckeep <- which(predictions %in% other.labels)
        predictions <- predictions[ckeep]
        for (n in names(all.markers)) {
            if (!(n %in% other.labels) && n != label) {
                all.markers[[n]] <- NULL
            }
        }
    } else if (anyNA(predictions)) {
        ckeep <- which(!is.na(predictions))
        predictions <- predictions[ckeep]
    }

    rkeep <- which(rownames(test) %in% unlist(all.markers, use.names=FALSE))
    test <- test[rkeep,ckeep,drop=FALSE]

    # Prioritize the markers with interesting variation in the test data for
    # visualization. If we only have one label, we use the most abundant markers.
    if (length(unique(predictions)) > 1L) {
        score.args$num.threads <- num.threads
        score.args$compute.auc <- (order.by.effect == "auc")
        score.args$compute.cohens.d <- (order.by.effect == "cohens.d")
        score.args$compute.delta.mea <- (order.by.effect == "delta.mean")
        score.args$compute.delta.detected <- (order.by.effect == "delta.detected")

        interesting <- do.call(scrapper::scoreMarkers, c(list(test, predictions), score.args))
        stats <- interesting[[order.by.effect]][[label]][[order.by.summary]]
        decreasing <- (order.by.summary!="min.rank")

    } else {
        stats <- rowMeans(test)
        decreasing <- TRUE
    }

    ro <- order(stats, decreasing=decreasing)
    co <- order(predictions)
    list(rows=rkeep[ro], columns=ckeep[co], predictions=predictions[co])
}
