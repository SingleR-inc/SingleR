#' Plot score distributions of labels.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param show A string selecting the scores representation to show.
#' @param labels String vector indicating one or more labels to show.
#' If \code{labels} is left \code{NULL}, all labels available in \code{results} are presented.
#' @param dots.on.top Logical which sets whether cell dots are plotted on top of, versus behind, the violin plots.
#' @param colors String vector that sets the colors.
#' Order of colors should be: `this label`, `this label - pruned`, `any label`.
#' @param size Scalar which sets the size of the dots
#' @param ncol Integer number of labels to display per row
#' 
#' @return a \link[ggplot2]{ggplot} object showing SingleR scores in a dot and/or violin plot representation.
#' 
#' @details
#' 
#' This function creates jitter and violin plots showing the representations of the scores of all cells across a single label or multiple labels,
#' and may be useful for visualizing and tuning the \code{nmads}, \code{min.diff.med}, and \code{min.diff.next} cutoffs of the \code{\link{pruneScores}} function.
#'  
#' The \code{show} input determines what representation to show.  Options are:
#' \itemize{
#' \item "scores" = shows the raw scores
#' \item "delta.med" = shows the difference between max score of a cell and the median of all scores for that cell
#' \item "delta.next" = shows the difference between best and second-best tuning scores of each cell
#' }
#' 
#' For a given label X, cells distributions in several categories are shown:
#' \itemize{
#' \item Was assigned to label X, and the label was not pruned away.
#' \item Was assigned to label X, and the label was pruned away.
#' \item For \code{show="scores"} only, Was assigned as any label, including label X.
#' }
#' Each category is grouped and colored separately.
#' When \code{show = "delta.med"} or \code{show = "delta.next"}, representations are only shown for a cell's final label.
#' 
#' Note that when \code{show="scores"} or \code{show="delta.med"}, the function focuses on initial scores only, i.e., prior to fine tuning.
#' However, the labels may be defined after fine-tuning in \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' Thus, the best score for an individual cell may not be its final label.
#'
#' @seealso
#' \code{\link{SingleR}}, to generate scores.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores.
#' 
#' @author Daniel Bunis
#' @examples
#' example(SingleR, echo=FALSE)
#' 
#' # To show the distribution of scores grouped by label:
#' plotScoreDistribution(results = pred)
#' # We can display a particular label using the label
#' plotScoreDistribution(results = pred, labels = "B")
#'
#' # To show the distribution of deltas between cells' maximum and median scores,
#' #   grouped by label, change `show` to "delta.med":
#' #   This is useful for checking/adjusting nmads and min.diff.med
#' plotScoreDistribution(results = pred, show = "delta.med")
#' # The nmads cutoff can be displayed using show.nmads.cutoff.
#' plotScoreDistribution(results = pred, show = "delta.med",
#'     show.nmads.cutoff = 3)
#'
#' # To show the distribution of deltas between cells' top 2 fine-tuning scores,
#' #   grouped by label, change `show` to "delta.next":
#' #   This is useful for checking/adjusting min.diff.next
#' plotScoreDistribution(results = pred, show = "delta.next")
#' 
#' @export
plotScoreDistribution <- function(results, show = "scores",
    labels = colnames(results$scores), size = 0.5, ncol = 5,
    dots.on.top = TRUE, colors = c("#F0E442", "#E69F00", "gray60")) {
    
    show <- match.arg(show, c("scores", "delta.med", "delta.next"))

    if (length(colors)<3) {
        stop("3 colors are expected.")
    }
    if (is.null(names(colors))) {
        names(colors) <- 
            c('this label', 'this label - pruned',
            'any label')
    }
    if (show == "scores") {
        labels <- labels[labels %in% colnames(results$scores)]
    } else {
        labels <- labels[labels %in% levels(as.factor(results$labels))]
    }

    # Gathere the scores data in a dataframe and trim to 'labels'
    df <- .scores_data_gather(results, show)
    df <- df[df$label %in% labels,]
    
    # Make the plot
    p <- ggplot2::ggplot(
            data = df,
            ggplot2::aes_string(
                x = "cell.calls", y = show, fill = "cell.calls")) + 
        ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(name = "Cell Calls", values = colors) +
        ggplot2::facet_wrap(facets = ~label, ncol = ncol)
    if (length(labels) == 1) {
        p <- p + ggplot2::scale_x_discrete(name = NULL, labels = NULL)   
    } else {
        p <- p + ggplot2::scale_x_discrete(name = "Labels", labels = NULL)
    }
    if (dots.on.top) {
        p <- p+ ggplot2::geom_violin()
    }
    p <- p + ggplot2::geom_jitter(
        height = 0, width = 0.3, color = "black", shape = 16,size = size)
    if (!dots.on.top) {
        p <- p + ggplot2::geom_violin()
    }
    
    p
}

#' @importFrom DelayedMatrixStats rowMedians 
#' @importFrom DelayedArray DelayedArray rowMaxs
.scores_data_gather <- function(results, show) {
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    
    scores <- results$scores
    maxed <- rowMaxs(DelayedArray(scores))
    delta <- maxed - DelayedMatrixStats::rowMedians(DelayedArray(scores))
    
    # Create a dataframe with separate rows for each score in scores.
    df <- data.frame(
        #cell id of the cell
        id = rep(rownames(results), each=ncol(scores)),
        #final call of the cell
        called = rep(results$labels, each=ncol(scores)),
        #label of the current score
        label = rep(
            colnames(results$scores),
            nrow(results)),
        scores = as.numeric(t(scores)),
        delta.med = as.numeric(t(delta)),
        stringsAsFactors = FALSE)
    
    if (!is.null(results$tuning.scores)){
        df$delta.next <- 
            results$tuning.scores$first - results$tuning.scores$second
    }
    
    # Add whether this label is the final label given to each cell.
    df$cell.calls <- "any label"
    df$cell.calls[df$label == df$called] <- "this label"
    
    if (!is.null(results$pruned.labels)){
        # Retrieve if cells' calls were scored as to be prunes versus not,
        #  then add this to df$cell.calls, but only when =="this label"
        prune.calls <- is.na(results$pruned.labels)
        prune.string <- as.character(prune.calls)
        prune.string[prune.calls] <- " - pruned"
        prune.string[!prune.calls] <- ""
        df$cell.calls[df$cell.calls=="this label"] <- paste0(
            df$cell.calls[df$cell.calls=="this label"],
            rep(prune.string, each=ncol(scores))[df$cell.calls=="this label"])
    }
    
    #Duplicate the "this label" data, and change df$cell.calls to "any label"
    df.thisLab <- 
        df[df$cell.calls %in% c("this label", "this label - pruned"),
           ,drop=FALSE]
    df$cell.calls <- 'any label'
    
    if (show %in% c("delta.med", "delta.next")) {
        df <- df.thisLab
    } else {
        df <- rbind(df.thisLab, df)
    }
    
    df
}
