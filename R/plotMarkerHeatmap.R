#' Plot a heatmap of the markers for a label
#'
#' Create a heatmap of the log-normalized expression for the most interesting markers of a particular label.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}},
#' \code{\link{classifySingleR}}, \code{\link{combineCommonResults}}, or \code{\link{combineRecomputedResults}}.
#' @param test A numeric matrix of log-normalized expression values where rows are genes and columns are cells.
#' Each row should be named with the gene name.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param label String specifying the label of interest.
#' @param other.labels Character vector specifying the other labels to be compared to \code{label} when finding interesting markers.
#' Defaults to all available labels.
#' @param assay.type Integer scalar or string specifying the matrix of expression values to use if \code{test} is a \linkS4class{SummarizedExperiment}.
#' @param use.pruned Logical scalar indicating whether the pruned labels should be used instead.
#' @param order.by String specifying the column of the output of \code{\link[scran]{scoreMarkers}} with which to sort for interesting markers.
#' @param top Integer scalar indicating the top most interesting markers to show in the heatmap.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization scheme to use for marker scoring.
#'
#' @return 
#' The output of \code{\link[pheatmap]{pheatmap}} is returned showing the heatmap on the current graphics device.
#'
#' @details
#' This function creates a heatmap where each row is a marker gene for \code{label} and each column is a cell in the test dataset.
#' The aim is to check the effectiveness of the reference-derived markers for distinguishing between labels in the test dataset.
#' Useful markers from the reference should show strong upregulation in \code{label} compared to all \code{other.labels}. 
#' We identify such markers by scoring all reference-derived markers with \code{\link[scran]{scoreMarkers}} on the \code{test} expression.
#' The \code{top} markers based on the specified \code{order.by} field are shown in the heatmap.
#' If only one label is present, markers are ranked by average abundance intead. 
#'
#' @author Aaron Lun
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' plotMarkerHeatmap(pred, test, pred$labels[1])
#' plotMarkerHeatmap(pred, test, pred$labels[1], use.pruned=TRUE)
#' plotMarkerHeatmap(pred, test, pred$labels[1], order.by="rank.AUC")
#' 
#' @export
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix rowMeans
#' @importFrom utils head
plotMarkerHeatmap <- function(
    results,
    test,
    label,
    other.labels=NULL,
    assay.type="logcounts",
    use.pruned=FALSE,
    order.by="rank.logFC.cohen",
    top=20,
    BPPARAM = SerialParam()) 
{
    test <- .to_clean_matrix(test, assay.type, check.missing=FALSE)
    all.markers <- metadata(results)$de.genes[[label]]

    if (use.pruned) {
        labfield <- "pruned.labels"
    } else {
        labfield <- "labels"
    }
    predictions <- results[[labfield]]

    ckeep <- seq_len(ncol(test))
    if (!is.null(other.labels)) {
        ckeep <- predictions %in% other.labels
        predictions <- predictions[ckeep]
        for (n in names(all.markers)) {
            if (!(n %in% other.labels) && n != label) {
                all.markers[[n]] <- NULL
            }
        }
    } else if (anyNA(predictions)) {
        ckeep <- !is.na(predictions)
        predictions <- predictions[ckeep]
    }

    rkeep <- rownames(test) %in% unlist(all.markers)
    test <- test[rkeep,ckeep,drop=FALSE]

    # Prioritize the markers with interesting variation in the test data for
    # visualization. If we only have one label, we use the most abundant markers.
    if (length(unique(predictions)) > 1L) {
        interesting <- scran::scoreMarkers(test, predictions, BPPARAM=BPPARAM)
        stats <- interesting[[label]]
        o <- order(stats[[order.by]], decreasing=!startsWith(order.by, "rank."))
        to.show <- rownames(test)[o]
    } else {
        abundance <- rowMeans(test)
        to.show <- names(abundance)[order(abundance, decreasing=TRUE)]
    }

    to.show <- head(to.show, top)
    limits <- range(test, na.rm=TRUE)
    colnames(test) <- seq_len(ncol(test))
    col.order <- order(predictions)

    pheatmap::pheatmap(
        test[to.show,col.order,drop=FALSE],
        breaks=seq(limits[1], limits[2], length.out=26),
        color=viridis::viridis(25),
        annotation_col=data.frame(labels=predictions, row.names=colnames(test)),
        cluster_col=FALSE,
        show_colnames=FALSE
    )
}
