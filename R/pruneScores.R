#' Prune out low-quality assignments
#'
#' Remove low-quality assignments based on the cell-label score matrix returned by \code{\link{classifySingleR}}.
#'
#' @param scores SingleR output, or a Matrix of scores per cell (row) and label (column) generated during assignment by \code{\link{classifySingleR}}.
#' @param min.diff.vs.median Numeric scalar specifying the minimum difference of each cell's maximum score from the median score.
#' @param min.percent.diff.vs.next Numeric scalar specifying the minimum percent difference between the best score and the next best score.
#' @param nmads Numeric scalar specifying the number of MADs to use to define cells with low outlier scores per label.
#' 
#' @return A logical vector specifying which assignments should be ignored.
#'
#' @details
#' By itself, the SingleR algorithm will always assign a label to every cell.
#' This occurs even if the cell's true label is not represented in the reference set of labels,
#' resulting in assignment of an incorrect label to that cell.
#' The \code{pruneScores} function aims to mitigate this effect by removing poor-quality assignments with \dQuote{low} scores.
#'
#' We define a low score using two orthogonal measures, per-cell and per-label.
#' For the per-cell check, we see whether the maximum score is less than \code{min.diff} above the median score for each cell.
#' If so, this indicates that the cell matches all labels with the same confidence 
#' such that the one reported label is not particularly meaningful.
#' This can also be justified in high-dimensional analyses where,
#' in the absence of any strong similarity to a single label, all distances converge to the same value.
#'  
#' For the per-label check, we identify cells that are small outliers in the distribution of scores for each label.
#' (This only includes cells that pass the per-cell check.)
#' Specifically, cells that are more than \code{nmads} below the median score for each label are ignored.
#' This assumes that most cells are correctly assigned to their true label
#' and that cells of the same label have a smooth distribution of scores.
#' Thus, small outliers represent poor-quality assignments from a distinct subpopulation that should be removed.
#'
#' The defaults for both of these parameters are largely arbitrary and chosen based on experience.
#' Smaller values for \code{min.diff} and larger values for \code{nmads} will reduce the stringency of the pruning.
#' In particular:
#' \itemize{
#' \item The scores do not consider the effects of fine-tuning (as scores are not comparable across different fine-tuning steps).
#' In situations involving a majority of labels with only subtle distinctions,
#' it is possible for the scores to be relatively similar but for the labels to be correctly assigned after fine-tuning.
#' In such cases, the default setting of \code{min.diff} may be too stringent.
#' \item It is possible for the per-label score distribution to be multimodal yet still correct,
#' e.g., due to cells belong to subtypes when the (correct) label corresponds to the main type.
#' In such cases, the default \code{nmads} may be too stringent as it will remove minor subpopulations with low scores.
#' }
#' Note that decreasing \code{min.diff} may actually \emph{increase} the stringency of the per-label check, 
#' depending on whether the additional retained cells decrease the MAD.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{classifySingleR}}, to generate \code{scores}.
#' 
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' summary(pruneScores(pred$scores))
#'
#' # Less stringent:
#' summary(pruneScores(pred$scores, min.diff=0))
#' summary(pruneScores(pred$scores, nmads=5))
#' 
#' @export
#' @importFrom DelayedMatrixStats rowMedians 
#' @importFrom DelayedArray DelayedArray rowMaxs
#' @importFrom stats median mad
pruneScores <- function(scores, min.diff.vs.median = 0.05, nmads=3, min.percent.diff.vs.next = 0.05) {
    if (is(scores, "DataFrame")){
      scores <- scores$scores
    }

    maxed <- rowMaxs(DelayedArray(scores))
    delta <- maxed - rowMedians(DelayedArray(scores))
    keep <- delta >= min.diff.vs.median
    keep <- keep & (rowSums(scores/maxed > (1-min.percent.diff.vs.next)) < 2)

    best <- max.col(scores)
    by.label <- split(which(keep), best[keep])

    for (l in by.label) {
        current <- maxed[l]
        med <- median(current)
        MAD <- mad(current, center=med)
        keep[l] <- (current >= med - nmads * MAD)
    }

    !keep
}
