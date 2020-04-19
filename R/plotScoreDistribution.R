#' Plot score distributions of labels.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param show String specifying whether to show the scores, the difference from the median or the difference from the next-best score.
#' @param labels.use String vector indicating one or more labels to show.
#' If \code{NULL}, all labels available in \code{results} are presented.
#' @param show.all.originals Logical which sets whether separate sets of plots for all original scores matrices shown be shown when \code{results} is the output of a multiple reference \code{\link{SingleR}} run or of the \code{\link{combineCommonResults}} or \code{\link{combineRecomputedResults}} functions.
#' @param calls.use,pruned.use Integer which sets which label calls or pruning calls to use when \code{results} is the output of a multiple reference \code{\link{SingleR}} run or of the \code{\link{combineCommonResults}} or \code{\link{combineRecomputedResults}} functions
#' A value of 0 points to the overall results, while any other integer indicates the index of the individual output that should be targetted.
#'
#' @param scores.use Integer which sets which scores to use when \code{results} is the output of a multiple reference \code{\link{SingleR}} run or of the \code{\link{combineCommonResults}} or \code{\link{combineRecomputedResults}} functions
#' A value of 0 points to the overall results, while any other integer indicates the index of the individual output that should be targetted.
#'
#' Alternatively, \code{scores.use} can be an integer vector pointing to multiple original results dataframes.
#'
#' Overwritten when \code{show.all.originals = TRUE}. Set \code{show.all.originals = FALSE} to adjust.
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
#' @param grid.vars named list of extra variables to pass to \code{\link[gridExtra]{grid.arrange}}
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
#' plotScoreDistribution(results = pred, labels.use = "B")
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
#' ### Multi-Reference Compatibility ###
#'
#' # When SingleR is run with multiple references, default output will contain
#' # separate plots for each original reference, as well as for the the combined
#' # set when 'show' = "scores".
#' example(combineRecomputedResults, echo = FALSE)
#' plotScoreDistribution(results = combined, show = "scores")
#' plotScoreDistribution(results = combined, show = "delta.med")
#' plotScoreDistribution(results = combined, show = "delta.next")
#'
#' # To color and group cells by non-final label and pruned calls,
#' # use 'calls.use' and 'pruned.use'
#' plotScoreDistribution(results = combined, show = "delta.med",
#'     calls.use = 1, pruned.use = 1)
#' plotScoreDistribution(results = combined, show = "delta.next",
#'     calls.use = 1, pruned.use = 1)
#'
#' # To instead target only final or a certain reference's scores,
#' # use 'scores.use' and add 'show.all.originals = FALSE'
#' plotScoreDistribution(results = combined, show = "scores",
#'     show.all.originals = FALSE, scores.use = 1)
#'
#' plotScoreDistribution(results = combined, show = "delta.med",
#'     show.all.originals = FALSE, scores.use = 1,
#'     show.nmads = 3,
#'     show.min.diff = 0.03)
#' plotScoreDistribution(results = combined, show = "delta.next",
#'     show.all.originals = FALSE, scores.use = 1,
#'     show.min.diff = 0.03)
#'
#' @export
plotScoreDistribution <- function(
    results,
    show = c("delta.med", "delta.next", "scores"),
    labels.use = colnames(results$scores),
    show.all.originals = TRUE,
    scores.use = 0,
    calls.use = 0,
    pruned.use = calls.use,
    size = 0.5,
    ncol = 5,
    dots.on.top = TRUE,
    this.color = "#F0E442",
    pruned.color = "#E69F00",
    other.color = "gray60",
    show.nmads = NULL,
    show.min.diff = NULL,
    grid.vars = list())
{

    ### For multi-ref results
    if (!is.null(results$orig.results)) {
        if (show.all.originals) {
            scores.use <- c(seq_along(results$orig.results))
            if (show == "scores") {
                scores.use <- c(0,scores.use)
            }
        }
        # When scores.use is a vector, recursively make calls for each and return plots in a grid.
        if (length(scores.use)>1) {
            plots <- lapply(
                scores.use,
                function(this) {
                    plotScoreDistribution(results, show, labels.use, FALSE, this,
                        calls.use, pruned.use, size, ncol, dots.on.top,
                        this.color, pruned.color, other.color,
                        show.nmads, show.min.diff, list())
                })
            grid.vars <- c(grid.vars, grobs = plots)

            return(do.call(gridExtra::grid.arrange, grid.vars))
        }
    }

    show <- match.arg(show)
    if (show == "delta.next" && calls.use != scores.use && !is.null(results$orig.results)) {
        message("'calls.use', and 'scores.use' must be the same for 'show=\"delta.next\"'. 'calls.use' was updated to ", scores.use, ".")
        calls.use <- scores.use
    }

    ## Extract data
    score.res <- .grab_results(results, scores.use)
    scores <- score.res$scores
    tuning.scores <- score.res$tuning.scores
    scores.title <- .values_title(results, scores.use, show)

    calls <- .grab_results(results, calls.use)$labels
    calls.title <- .calls_title(results, calls.use, "Calls", show, scores.use)

    prune.calls <- .grab_results(results, calls.use)$pruned.labels
    prune.label <- .prune_label(results, pruned.use, calls.use)

    # Make & trim data frame
    if (show!="delta.next") {
        df <- .scores_df_gather(
            show, scores, calls, prune.calls, prune.label)
    } else {
        df <- .next_df_gather(
            tuning.scores, calls, prune.calls, prune.label)
    }
    labels.use <- labels.use[labels.use %in% df$label]
    df <- df[df$label %in% labels.use,]

    ### Creating the plot
    p <- .plot_score_distribution(
        df, calls.title, scores.title, prune.label,
        this.color, pruned.color, other.color, size, ncol, dots.on.top)
    if (grepl("delta",show) && !(is.null(show.nmads)) || !(is.null(show.min.diff))) {
        p <- .add_cutoff_lines(
            p, score.res, df, show, show.nmads, show.min.diff)
    }

    p
}

#' @importFrom DelayedMatrixStats rowMedians
.scores_df_gather <- function(
    show, values, calls, prune.calls, prune.label) {

    if (show=="delta.med") {
        values <- values - DelayedMatrixStats::rowMedians(DelayedArray(values), na.rm = TRUE)
    }

    # Create a dataframe with separate rows for each score in values.
    df <- data.frame(
        label = rep(colnames(values), nrow(values)),
        values = as.numeric(t(values)),
        stringsAsFactors = FALSE)

    # Add whether this label is the final label given to each cell.
    df$cell.calls <- rep("other", nrow(df)) # rep() protects when nrow(df)=0.
    is.called <- df$label == rep(calls, each=ncol(values))
    df$cell.calls[is.called] <- "assigned"
    # Replace cell.call if cell was pruned.
    if (!is.null(prune.calls)) {
        is.pruned <- rep(is.na(prune.calls), each=ncol(values))
        df$cell.calls[is.pruned & is.called] <- prune.label
    }

    df
}

.next_df_gather <- function(
    tuning.scores, calls, prune.calls, prune.label) {

    if (is.null(tuning.scores)) {
        stop("Target 'results' lacks fine-tuning diagnostics for 'show=\"delta.next\"'")
    }

    df <- data.frame(
        values = tuning.scores$first - tuning.scores$second,
        label = calls,
        cell.calls = rep("assigned", nrow(tuning.scores)), # don't unrep, protects when nrow(tuning.scores)=0.
        stringsAsFactors = FALSE)

    if (!is.null(prune.calls)) {
        df$cell.calls[is.na(prune.calls)] <- prune.label
    }

    df
}

.plot_score_distribution <- function(
    df,
    calls.title, scores.title, prune.label,
    this.color, pruned.color, other.color, size, ncol, dots.on.top) {

    p <- ggplot2::ggplot(data = df,
            ggplot2::aes_string(x = "cell.calls", y = "values", fill = "cell.calls")) +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(
            name = calls.title,
            breaks = c("assigned", prune.label, "other"),
            values = c(this.color, pruned.color, other.color)) +
        ggplot2::facet_wrap(facets = ~label, ncol = ncol) +
        ggplot2::ylab(scores.title)

    if (length(levels(as.factor(df$label))) == 1) {
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

    if (dots.on.top) {
        p <- p + ggplot2::geom_jitter(
            height = 0, width = 0.3, color = "black", shape = 16,size = size,
            na.rm = TRUE)
    }

    p
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

# Sets the Title for the color legend based on which ref's pruning calls are shown
.prune_label <- function(results, pruned.use, calls.use){
    if (pruned.use!=calls.use) {
        which_pruned <- switch(as.character(pruned.use==0),
            "TRUE" = ifelse(is.null(results$orig.results), "", "(in Final)"),
            "FALSE" = paste0("(in Ref #", pruned.use, ")"))
        paste("pruned", which_pruned)
    } else {
        "pruned"
    }
}

