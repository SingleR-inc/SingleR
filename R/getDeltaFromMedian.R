#' Compute the difference from median
#'
#' Compute the delta value for each cell,
#' defined as the difference between the score for the assigned label and the and median score across all labels.
#' 
#' @inheritParams pruneScores
#'
#' @details
#' This funciton computes the same delta value that is used in \code{\link{pruneScores}},
#' for users who want to apply more custom filters or visualizations.
#'
#' @author Aaron Lun
#'
#' @return A numeric vector containing delta values for each cell in \code{results}.
#'
#' @seealso
#' \code{\link{pruneScores}}, where the delta values are used.
#'
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' summary(getDeltaFromMedian(pred))
#' 
#' @export
#' @importFrom beachmat initializeCpp tatami.row.medians
getDeltaFromMedian <- function(results) {
    scores <- results$scores
    labels <- results$labels
    assigned <- scores[cbind(seq_along(labels), match(labels, colnames(scores)))]
    assigned - tatami.row.medians(initializeCpp(scores), num.threads=1)
}
