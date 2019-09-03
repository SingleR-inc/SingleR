#' Plot score distributions of labels.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param show String specifying whether to show the scores, the difference from the median or the difference from the next-best score.
#' @param labels String vector indicating one or more labels to show.
#' If \code{NULL}, all labels available in \code{results} are presented.
#' @param dots.on.top Logical specifying whether cell dots should be plotted on top of the violin plots.
#' @param this.color String specifying the color for cells that were assigned to the label.
#' @param pruned.color String specifying the color for cells that were assigned to the label but pruned.
#' @param other.color String specifying the color for other cells not assigned to the label.
#' Only used when \code{show="scores"}.
#' @param size Numeric scalar to set the size of the dots.
#' @param ncol Integer scalar to set the number of labels to display per row.
#' @param show.nmads Numeric scalar that shows the threshold that would be used for pruning with \code{\link{pruneScores}}.
#' Only used when \code{show="delta.med"}.
#' 
#' @return A \link[ggplot2]{ggplot} object showing assignment scores in violin plots.
#' 
#' @details
#' This function creates jitter and violin plots showing assignment scores or related values for all cells across one or more labels.
#' It is intended for visualizing and adjusting the \code{nmads}, \code{min.diff.med}, and \code{min.diff.next} cutoffs of the \code{\link{pruneScores}} function.
#'  
#' The \code{show} argument determines what values to show on the y-axis.
#' Options are:
#' \itemize{
#' \item \code{"delta.med"}, the difference between the score of the assigned label and the median of all scores for each cell.
#' \item \code{"delta.next"}, the difference between best and second-best tuning scores of each cell.
#' \item \code{"scores"}, the raw assignment scores prior to fine-tuning.
#' }
#' 
#' For a given label X, cells distributions in several categories are shown:
#' \itemize{
#' \item Was assigned to label X, and the label was not pruned away.
#' \item Was assigned to label X, and the label was pruned away.
#' \item Was assigned as any label, including label X.
#' }
#' Each category is grouped and colored separately based on \code{this.color} and related parameters.
#'
#' Values are stratified according to the assigned labels in \code{results$labels}.
#' If any fine-tuning was performed, the highest scoring label for an individual cell may not be its final label.
#' This may manifest as negative values when \code{show="delta.med"}.
#' 
#' Also note that \code{\link{pruneScores}} trims based on the \code{min.diff.med} and \code{min.diff.next} cutoffs first,
#' before calculating the first-labels' delta medians.
#' Thus, the actual \code{nmads} cut-off used in \code{\link{pruneScores}} may vary from the one portrayed in the plot.
#'
#' @seealso
#' \code{\link{SingleR}}, to generate scores.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores, and to see more about the quailty cutoffs.
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
#' # The nmads cutoff can be displayed using show.nmads.
#' plotScoreDistribution(results = pred, show = "delta.med",
#'     show.nmads = 3)
#' # A min.diff.med cutoff can be shown using show.min.diff
#' plotScoreDistribution(results = pred, show = "delta.med",
#'     min.diff = 0.03)
#'
#' # To show the distribution of deltas between cells' top 2 fine-tuning scores,
#' #   grouped by label, change `show` to "delta.next":
#' #   This is useful for checking/adjusting min.diff.next
#' plotScoreDistribution(results = pred, show = "delta.next")
#' # A min.diff.med cutoff can be shown using show.min.diff
#' plotScoreDistribution(results = pred, show = "delta.next",
#'     min.diff = 0.03)
#' 
#' @export
plotScoreDistribution <- function(results,
    show = c("delta.med", "delta.next", "scores"),
    labels = colnames(results$scores), 
    size = 0.5, ncol = 5, dots.on.top = TRUE, 
    this.color = "#F0E442", 
    pruned.color = "#E69F00", 
    other.color = "gray60",
    show.nmads = NULL, show.min.diff = NULL) 
{
    show <- match.arg(show)
    if (show!="delta.next") {
        df <- .scores_data_gather(results, show)
    } else {
        df <- .next_data_gather(results)
    }

    labels <- labels[labels %in% colnames(results$scores)]
    df <- df[df$label %in% labels,]

    # Creating the plot object:
    p <- ggplot2::ggplot(data = df, 
            ggplot2::aes_string(x = "cell.calls", y = "values", fill = "cell.calls")) + 
        ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(name = "Cell Calls", 
            values = c(assigned=this.color, pruned=pruned.color, other=other.color)) +
        ggplot2::facet_wrap(facets = ~label, ncol = ncol) +
        ggplot2::ylab(show)

    if (length(labels) == 1) {
        p <- p + ggplot2::scale_x_discrete(name = NULL, labels = NULL)   
    } else {
        p <- p + ggplot2::scale_x_discrete(name = "Labels", labels = NULL)
    }
    
    # Adding the frills:
    if (!dots.on.top) {
        p <- p + ggplot2::geom_jitter(
            height = 0, width = 0.3, color = "black", shape = 16,size = size)
    }
    p <- p + ggplot2::geom_violin(na.rm=TRUE)

    if (grepl("delta",show) && !(is.null(show.nmads)) || !(is.null(show.min.diff))) {
        p <- .add_cutoff_lines(p, results, df, show, show.nmads, show.min.diff)
    }

    if (dots.on.top) {
        p <- p + ggplot2::geom_jitter(
            height = 0, width = 0.3, color = "black", shape = 16,size = size)
    }
    
    p
}

#' @importFrom DelayedMatrixStats rowMedians 
.scores_data_gather <- function(results, show) {
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }

    values <- results$scores
    if (show=="delta.med") {
        values <- values - DelayedMatrixStats::rowMedians(DelayedArray(values))
    }
    
    # Create a dataframe with separate rows for each score in values.
    df <- data.frame(
        id = rep(rownames(results), each=ncol(values)),
        called = rep(results$labels, each=ncol(values)),
        label = rep(colnames(results$scores), nrow(results)),
        values = as.numeric(t(values)),
        stringsAsFactors = FALSE)
    
    # Add whether this label is the final label given to each cell.
    df$cell.calls <- rep("other", nrow(df)) # rep() protects when nrow(df)=0.
    df$cell.calls[df$label == df$called] <- "assigned"
    
    if (!is.null(results$pruned.labels)) {
        is.pruned <- rep(is.na(results$pruned.labels), each=ncol(values))
        df$cell.calls[is.pruned & df$cell.calls=="assigned"] <- "pruned"
    }
        
    df
}

.next_data_gather <- function(results) {
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    if (is.null(results$tuning.scores)) {
        stop("'results' lacks fine-tuning diagnostics for 'show=\"delta.next\"'")
    }

    df <- data.frame(
        id = rownames(results), 
        called = results$labels,
        value = results$tuning.scores$first - results$tuning.scores$second,
        cell.calls = rep("assigned", nrow(results)), # don't unrep, protects when nrow(results)=0.
        stringsAsFactors = FALSE)

    if (!is.null(results$pruned.labels)) {
        df$cell.calls[is.na(results$pruned.labels)] <- "pruned"
    }

    df
}

.add_cutoff_lines <- function(
    p, results, df, show, show.nmads, show.min.diff) {

    if (show == "delta.med" && !(is.null(show.nmads))) {
        # Add cutoff line for nmads cutoff
            # Uses geom_error (error bars, but with ymin = ymax)
        p <- .add_nmads_lines(p, results, df, show.nmads)
    }

    if (grepl("delta",show) && !(is.null(show.min.diff))) {
        # Add cutoff line for min.diff cutoff
        df_min.diff <- data.frame(color = "2.min.diff")
        p <- p + 
            ggplot2::geom_hline(
                data = df_min.diff, na.rm = TRUE, size = 1.1,
                ggplot2::aes(yintercept = show.min.diff, color = color))
    }
    # Set the colors and a=labels for combined legend.
    p <- p + ggplot2::scale_color_manual(
        name = "Lines",
        values = c(
            '1.nmad' = "#0072B2",
            '2.min.diff' = "red"),
        labels = c(
            '1.nmad' = "nmad cutoff",
            '2.min.diff' = "min.diff cutoff")
        ) + 
        ggplot2::guides(
            color = ggplot2::guide_legend(override.aes = list(size=2)))

    p
}

#' @importFrom stats median mad
.add_nmads_lines <- function(p, results, df, show.nmads) {
    labels <- levels(as.factor(df$label))
    
    if (is.null(results$first.labels)) {
        df$first.labels <- results$labels
    } else {
        df$first.labels <- results$first.labels
    }
    
    # Calculate the nmad cutoff for each label, (ignoring min.diff cutoffs).
    df_bars <- data.frame(
        labels = labels,
        medians = vapply(
            labels,
            function (X) median(df$delta.med[df$first.labels == X]),
            FUN.VALUE = numeric(1))
        )
    df_bars$nmads.cutoff <- vapply(
        labels,
        function (X) df_bars$medians[df_bars$labels == X] - show.nmads *
            mad(df$delta.med[df$first.labels == X]),
        FUN.VALUE = numeric(1))
    # Add to main df
    df$nmads.cutoff <- df_bars$nmads.cutoff[match(df$label, df_bars$label)]
    
    # Add dummy var for setting color later utilizing ggplot scale_color.
    df$mad <- "1.nmad"
    
    p <- p +
        ggplot2::geom_errorbar(
            data = df, na.rm=TRUE,
            ggplot2::aes_string(
                x = "cell.calls", ymin= "nmads.cutoff", ymax= "nmads.cutoff",
                color = "mad"),
            width = 1, size = 1.1,
            show.legend = c(color = TRUE, fill = FALSE))
    p
}
