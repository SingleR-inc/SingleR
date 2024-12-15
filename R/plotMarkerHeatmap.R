#' Plot a heatmap of the markers for a label
#'
#' Create a heatmap of the log-normalized expression for the most interesting markers of a particular label.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}},
#' \code{\link{classifySingleR}}, \code{\link{combineCommonResults}}, or \code{\link{combineRecomputedResults}}.
#' @param test A numeric matrix of log-normalized expression values where rows are genes and columns are cells.
#' Each row should be named with the same gene name that was used to compute \code{results}.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param label String specifying the label of interest.
#' @param other.labels Character vector specifying the other labels to be compared to \code{label} when finding interesting markers.
#' Defaults to all available labels.
#' @param assay.type Integer scalar or string specifying the matrix of expression values to use if \code{test} is a \linkS4class{SummarizedExperiment}.
#' @param use.pruned Logical scalar indicating whether the pruned labels should be used instead.
#' @param order.by.effect String specifying the effect size from \code{\link[scrapper]{scoreMarkers}} with which to sort for interesting markers.
#' @param order.by.summary String specifying the summary statistic from \code{\link[scrapper]{scoreMarkers}} with which to sort for interesting markers.
#' @param display.row.names Character vector of length equal to the number of rows of \code{test},
#' containing the names of the features to show on the heatmap (e.g., to replace IDs with symbols).
#' If \code{NULL}, the existing row names of \code{test} are used.
#' @param top Integer scalar indicating the top most interesting markers to show in the heatmap.
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
#' This can be used to apply the same general approach with other plots, e.g., using functions from \pkg{scuttle} or \pkg{dittoSeq}.
#'
#' @author Aaron Lun
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' plotMarkerHeatmap(pred, test, pred$labels[1])
#' plotMarkerHeatmap(pred, test, pred$labels[1], use.pruned=TRUE)
#' plotMarkerHeatmap(pred, test, pred$labels[1], order.by.effect="auc")
#'
#' # Manually configuring a simpler heatmap by label:
#' config <- configureMarkerHeatmap(pred, test, pred$labels[1])
#' mat <- assay(test, "logcounts")[head(config$rows, 20), config$columns]
#' aggregated <- scuttle::summarizeAssayByGroup(mat, config$predictions)
#' pheatmap::pheatmap(assay(aggregated), cluster_col=FALSE)
#' 
#' @export
#' @importFrom BiocParallel SerialParam bpnworkers
#' @importFrom Matrix rowMeans
#' @importFrom utils head
plotMarkerHeatmap <- function(
    results,
    test,
    label,
    other.labels=NULL,
    assay.type="logcounts",
    display.row.names=NULL,
    use.pruned=FALSE,
    order.by.effect="cohens.d",
    order.by.summary="min.rank",
    top=20,
    num.threads=bpnworkers(BPPARAM),
    BPPARAM = SerialParam(),
    ...) 
{
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
        num.threads=num.threads
    )

    to.show <- head(config$rows, top)
    predictions <- config$predictions
    test <- test[to.show, config$columns, drop=FALSE]
    if (!is.null(display.row.names)) {
        rownames(test) <- display.row.names[to.show]
    }

    limits <- range(test, na.rm=TRUE)
    colnames(test) <- seq_len(ncol(test))
    pheatmap::pheatmap(
        test,
        breaks=seq(limits[1], limits[2], length.out=26),
        color=viridis::viridis(25),
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
    other.labels=NULL,
    assay.type="logcounts",
    use.pruned=FALSE,
    order.by.effect="cohens.d",
    order.by.summary="min.rank",
    num.threads=1)
{
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
        interesting <- scrapper::scoreMarkers(
            test,
            predictions,
            num.threads=num.threads,
            compute.auc=(order.by.effect=="auc"),
            compute.cohens.d=(order.by.effect=="cohens.d"),
            compute.delta.mean=(order.by.effect=="delta.mean"),
            compute.delta.detected=(order.by.effect=="delta.detected")
        )

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
