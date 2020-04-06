#' Plot score distributions of labels.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param show String specifying whether to show the scores, the difference from the median or the difference from the next-best score.
#' @param labels String vector indicating one or more labels to show.
#' If \code{NULL}, all labels available in \code{results} are presented.
#' @param scores.use,calls.use,pruned.use Integer which sets which scores, final label calls, and pruning calls to use when \code{results} is the output of a multiple reference \code{\link{SingleR}} or \code{\link{classifySingleR}} run.\
#' A value of 0 points to the overall results, while any other integer indicates the index of the original \code{ref} that should be targetted.
#' @param dots.on.top Logical specifying whether cell dots should be plotted on top of the violin plots.
#' @param this.color String specifying the color for cells that were assigned to the label.
#' @param pruned.color String specifying the color for cells that were assigned to the label but pruned.
#' @param other.color String specifying the color for other cells not assigned to the label.
#' @param size Numeric scalar to set the size of the dots.
#' @param ncol Integer scalar to set the number of labels to display per row.
#' @param show.nmads Numeric scalar that shows the threshold that would be used for pruning with \code{\link{pruneScores}}.
#' Only used when \code{show="delta.med"}.
#' @param show.min.diff Numeric scalar that shows the threshold that would be used for pruning with \code{\link{pruneScores}}.
#' Only used when \code{show="delta.med"} or \code{"delta.next"}.
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
#' @author Daniel Bunis and Aaron Lun
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
#'     show.min.diff = 0.03)
#'
#' # To show the distribution of deltas between cells' top 2 fine-tuning scores,
#' #   grouped by label, change `show` to "delta.next":
#' #   This is useful for checking/adjusting min.diff.next
#' plotScoreDistribution(results = pred, show = "delta.next")
#' # A min.diff.med cutoff can be shown using show.min.diff
#' plotScoreDistribution(results = pred, show = "delta.next",
#'     show.min.diff = 0.03)
#'
#' # When SingleR is run with multiple references, investigation of how scores
#' # and calls from particular references contributed to the final outcome can
#' # be achieved through utilization of calls.use, scores.use, and pruned.use
#' example(combineRecomputedResults, echo = FALSE)
#' plotScoreDistribution(results = combined, show = "scores", scores.use = 1)
#' plotScoreDistribution(results = combined, show = "delta.med", scores.use = 1,
#'     show.nmads = 3,
#'     show.min.diff = 0.03)
#' plotScoreDistribution(results = combined, show = "delta.next", scores.use = 1,
#'     show.min.diff = 0.03)
#'
#' @export
plotScoreDistribution <- function(
    results,
    show = c("delta.med", "delta.next", "scores"),
    labels = colnames(.grab_results(results, scores.use)$scores),
    scores.use = 0,
    calls.use = scores.use,
    pruned.use = calls.use,
    size = 0.5,
    ncol = 5,
    dots.on.top = TRUE,
    this.color = "#F0E442",
    pruned.color = "#E69F00",
    other.color = "gray60",
    show.nmads = NULL,
    show.min.diff = NULL)
{
    show <- match.arg(show)
    if (show!="delta.next") {
        df <- .scores_data_gather(results, show, scores.use, calls.use, pruned.use)
    } else {
        df <- .next_data_gather(results, scores.use, calls.use, pruned.use)
    }

    labels <- labels[labels %in% colnames(.grab_results(results, scores.use)$scores)]
    df <- df[df$label %in% labels,]

    # Creating the plot object:
    p <- ggplot2::ggplot(data = df,
            ggplot2::aes_string(x = "cell.calls", y = "values", fill = "cell.calls")) +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(name = .legend_title(calls.use),
            values = c(assigned=this.color, pruned=pruned.color, other=other.color)) +
        ggplot2::scale_y_continuous(name = .y_title(show, scores.use)) +
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
            height = 0, width = 0.3, color = "black", shape = 16,size = size,
            na.rm = TRUE)
    }
    p <- p + ggplot2::geom_violin(na.rm=TRUE)

    if (grepl("delta",show) && !(is.null(show.nmads)) || !(is.null(show.min.diff))) {
        p <- .add_cutoff_lines(
            p, .grab_results(results, scores.use),
            df, show, show.nmads, show.min.diff)
    }

    if (dots.on.top) {
        p <- p + ggplot2::geom_jitter(
            height = 0, width = 0.3, color = "black", shape = 16,size = size,
            na.rm = TRUE)
    }

    p
}

.grab_results <- function(results, index) {
    if (index == 0) {
        return(results)
    } else {
        return(results$orig.results[[index]])
    }
}

.ensure_named <- function(results) {
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    results
}

#' @importFrom DelayedMatrixStats rowMedians
.scores_data_gather <- function(results, show, scores.use, calls.use, pruned.use) {

    score.res <- .ensure_named(.grab_results(results, scores.use))
    call.res <- .ensure_named(.grab_results(results, calls.use))
    pruned.res <- .ensure_named(.grab_results(results, pruned.use))

    values <- score.res$scores
    if (show=="delta.med") {
        values <- values - DelayedMatrixStats::rowMedians(DelayedArray(values), na.rm = TRUE)
    }

    # Create a dataframe with separate rows for each score in values.
    df <- data.frame(
        label = rep(colnames(score.res$scores), nrow(score.res)),
        values = as.numeric(t(values)),
        stringsAsFactors = FALSE)

    # Add whether this label is the final label given to each cell.
    df$cell.calls <- rep("other", nrow(df)) # rep() protects when nrow(df)=0.
    is.called <- df$label == rep(call.res$labels, each=ncol(values))
    df$cell.calls[is.called] <- "assigned"

    if (!is.null(pruned.res$pruned.labels)) {
        is.pruned <- rep(is.na(pruned.res$pruned.labels), each=ncol(values))
        df$cell.calls[is.pruned & is.called] <- .prune_label(pruned.use, calls.use)
    }

    df
}

.next_data_gather <- function(results, scores.use, calls.use, pruned.use) {

    score.res <- .ensure_named(.grab_results(results, scores.use))
    call.res <- .ensure_named(.grab_results(results, calls.use))
    pruned.res <- .ensure_named(.grab_results(results, pruned.use))

    if (is.null(score.res$tuning.scores)) {
        stop("Target 'results' lacks fine-tuning diagnostics for 'show=\"delta.next\"'")
    }
    if (calls.use != scores.use || pruned.use != scores.use) {
        stop("'calls.use', 'pruned.use', and 'scores.use' should be the same for 'show=\"delta.next\"'")
    }

    df <- data.frame(
        values = score.res$tuning.scores$first - score.res$tuning.scores$second,
        label = call.res$labels,
        cell.calls = rep("assigned", nrow(results)), # don't unrep, protects when nrow(results)=0.
        stringsAsFactors = FALSE)

    if (!is.null(pruned.res$pruned.labels)) {
        df$cell.calls[is.na(pruned.res$pruned.labels)] <- .prune_label(pruned.use, calls.use)
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
        df_min.diff <- data.frame(color = "2.min.diff", show.min.diff = show.min.diff)
        p <- p +
            ggplot2::geom_hline(
                data = df_min.diff, na.rm = TRUE, size = 1.1,
                ggplot2::aes_string(yintercept = "show.min.diff", color = "color"))
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
    to.show <- pruneScores(results, get.thresholds=TRUE)

    df <- data.frame(label=names(to.show), nmads.cutoff=to.show,
        values=0, # to keep the inherited aes_string() happy.
        cell.calls=rep("assigned", length(to.show)),
        mad=rep("1.nmad", length(to.show))) # Add dummy var for setting color later.

    p + ggplot2::geom_errorbar(data = df, na.rm=TRUE,
        ggplot2::aes_string(x = "cell.calls", ymin= "nmads.cutoff", ymax= "nmads.cutoff", color = "mad"),
        width = 1, size = 1.1, show.legend = c(color = TRUE, fill = FALSE))
}

# Sets the Title for the color legend based on which ref's calls are shown
.legend_title <- function(calls.use){
    switch(as.character(calls.use==0),
        "TRUE" = "Final Calls",
        "FALSE" = paste0("Ref #", calls.use, " Calls"))
}

# Sets the Title for the scores axis based which ref's scores are shown
.y_title <- function(show, scores.use){
    score_bit <- switch(as.character(scores.use==0),
        "TRUE" = "Final Scores",
        "FALSE" = paste0("(Ref #", scores.use, ")"))
    paste(show, score_bit, sep = ", ")
}

# Sets the Title for the color legend based on which ref's pruning calls are shown
.prune_label <- function(pruned.use, calls.use){
    if (pruned.use!=calls.use) {
        which_pruned <- switch(as.character(pruned.use==0),
            "TRUE" = "(in Final)",
            "FALSE" = paste0("(in Ref #", pruned.use, ")"))
        paste("pruned", which_pruned)
    } else {
        "pruned"
    }
}

