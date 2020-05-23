#' Plot delta distributions 
#'
#' Plot the distribution of deltas 
#' (i.e., the gap between the assignment score for the assigned label and those of the remaining labels)
#' across cells assigned to each reference label.
#'
#' @inheritParams plotScoreDistribution
#' @param deltas.use Integer scalar or vector specifying the individual annotation result from which to take deltas.
#' This is only relevant for combined results, see Details.
#' @param show String specifying whether to show the difference from the median (\code{"delta.med"}) 
#' or the difference from the next-best score (\code{"delta.next"}).
#'
#' @return
#' If \code{deltas.use} specifies a single set of deltas,
#' a \link[ggplot2]{ggplot} object is returned showing the deltas in violin plots.
#'
#' If \code{deltas.use} specifies multiple deltas for a combined result,
#' multiple ggplot objects are generated in a grid on the current graphics device.
#' 
#' If \code{delta.use} specifies multiple deltas and \code{grid.vars} is set to \code{NULL},
#' a list is returned containing the ggplot objects for manual display.
#'
#' @details
#' This function creates jitter and violin plots showing the deltas for all cells across one or more labels.
#' The idea is to provide a visual diagnostic for the confidence of assignment of each cell to its label.
#' The \code{show} argument determines what values to show on the y-axis:
#' \itemize{
#' \item \code{"delta.med"}, the difference between the score of the assigned label and the median of all scores for each cell.
#' \item \code{"delta.next"}, the difference between best and second-best scores of each cell at the last round of fine-tuning.
#' }
#'
#' If any fine-tuning was performed, the highest scoring label for an individual cell may not be its final label.
#' This may manifest as negative values when \code{show="delta.med"}.
#' \code{show="delta.next"} is guaranteed to be positive but may be overly stringent for references involving very similar labels.
#'
#' Pruned calls are identified as \code{NA}s in the \code{pruned.labels} field in \code{results}.
#' Points corresponding to cells with pruned calls are colored by \code{pruned.color};
#' this can be disabled by setting \code{pruned.color=NA}.
#'
#' For combined results (see \code{?"\link{combine-predictions}"}),
#' this function can show both the combined and individual deltas or labels.
#' This is done using the \code{deltas.use} and \code{labels.use} arguments,
#' entries of which refer to columns of \code{results$orig.results}.
#' \code{labels.use} is also allowed to contain zero, in which case the combined labels are used.
#' For example:
#' \itemize{
#' \item If we set \code{deltas.use=2} and \code{labels.use=1},
#' we will plot the deltas from the second individual reference against the labels faceted from the first reference.
#' \item If we set \code{deltas.use=1:2} and \code{labels.use=0},
#' we will plot the deltas from first and second references (in separate plots) faceted by the combined labels.
#' \item By default, the function will create a plot of deltas each individual reference.
#' In each plot, the deltas are faceted by the combined labels. 
#' }
#'
#' @seealso
#' \code{\link{SingleR}}, to generate scores.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores, and to see more about the quailty cutoffs.
#'
#' \code{\link[gridExtra]{grid.arrange}}, for tweaks to the how plots are arranged when multiple are output together.
#'
#' @author Daniel Bunis and Aaron Lun
#' @examples
#' example(SingleR, echo=FALSE)
#'
#' # Showing the delta to the median:
#' plotDeltaDistribution(pred)
#'
#' # Showing the delta to the next-highest score:
#' plotDeltaDistribution(pred, show = "delta.next")
#'
#' # Multi-reference compatibility:
#' example(combineRecomputedResults, echo = FALSE)
#' plotDeltaDistribution(results = combined)
#'
#' # To have plots output in a grid rather than as separate pages, provide,
#' # a list of inputs for gridExtra::grid.arrange() to 'grids.vars'.
#' plotDeltaDistribution(combined, grid.vars = list(ncol = 1))
#'
#' # An empty list will use grid.arrange defaluts
#' plotDeltaDistribution(combined, grid.vars = list())
#'
#' @export
plotDeltaDistribution <- function(
    results,
    show = c("delta.med", "delta.next"),
    labels.use = colnames(results$scores),
    deltas.use = NULL,
    calls.use = 0,
    size = 2,
    ncol = 5,
    dots.on.top = TRUE,
    this.color = "#F0E442",
    pruned.color = "#E69F00",
    grid.vars = list())
{
    results <- .ensure_named(results)
    show <- match.arg(show)
    is.combined <- !is.null(results$orig.results)
    ref.names <- colnames(results$orig.results)

    if (is.null(deltas.use)) {
        if (is.combined) {
            # Combined 'delta.med'/'delta.next' are not worth showing.
            deltas.use <- seq_along(results$orig.results) 
        } else {
            deltas.use <- 0L
        }
    }

    calls.use <- rep(calls.use, length.out=length(deltas.use))
    if (is.combined && show == "delta.next" && any(calls.use != deltas.use)) {
        warning("updating 'calls.use' to be the same as 'deltas.use' for 'show=\"delta.next\"'")
        calls.use <- deltas.use
    }

    plots <- vector("list", length(deltas.use))
    for (i in seq_along(plots)) {

        # Pulling out the scores to use in this iteration.
        chosen.scores <- deltas.use[i]
        if (is.combined) {
            if (chosen.scores==0L) {
                stop("deltas cannot be shown for combined results.")
            }
            score.results <- results$orig.results[[chosen.scores]]
        } else {
            score.results <- results
        }
        scores.title <- .values_title(is.combined, chosen.scores, ref.names, show)

        # Computing the values that we want to show.
        if (show=="delta.med") {
            values <- getDeltaFromMedian(score.results)
        } else if (show=="delta.next") {
            values <- score.results$tuning.scores$first - score.results$tuning.scores$second
        }

        # Pulling out the labels to use in this iteration.
        chosen.calls <- calls.use[i]
        if (chosen.calls==0L) {
            call.results <- results
        } else {
            call.results <- results$orig.results[[chosen.calls]]
        }

        labels <- call.results$labels
        labels.title <- .values_title(is.combined, chosen.calls, ref.names, "Labels")

        # Checking if we need pruning.
        pruned <- NULL
        if (is.na(pruned.color)) {
            pruned <- is.na(call.results$pruned.labels)
        }

        # Actually creating the plot
        plots[[i]] <- .plot_delta_distribution(
            values=values, labels=labels, pruned=pruned, labels.use=labels.use,
            labels.title=labels.title, scores.title=scores.title, 
            this.color=this.color, pruned.color=pruned.color,
            size=size, ncol=ncol, dots.on.top=dots.on.top)
    }

    if (length(plots)==1L) {
        # Doing this to be consistent with raw ggplot output.
        plots[[1]]
    } else {
        if (!is.null(grid.vars) && length(deltas.use) > 1L) {
            do.call(gridExtra::grid.arrange, c(plots, grid.vars))
        } else {
            plots
        }
    }
}

.plot_delta_distribution <- function(
    values, labels, pruned, labels.use,
    labels.title, scores.title, 
    this.color, pruned.color,
    size, ncol, dots.on.top)
{
    df <- data.frame(values=values, 
        label=labels, 
        x=character(length(values)))

    aes <- list(y="values", x="x")
    if (is.null(pruned)) {
        df$pruned <- pruned
        aes$color_by <- "pruned"
    }

    # Trim dataframe by labels:
    if (any(keep <- labels.use %in% df$label)) {
        labels.use <- labels.use[keep]
        df <- df[df$label %in% labels.use,]
    } else {
        warning("ignoring 'labels.use' as it has no values in ", scores.title)
    }

    # Making the violin plots.
    p <- ggplot2::ggplot(data = df, do.call(ggplot2::aes_string, aes)) + 
        ggplot2::xlab("") + 
        ggplot2::scale_color_manual(
            name = labels.title,
            breaks = c("FALSE", "TRUE"),
            values = c(this.color, pruned.color))

    p <- .pretty_violins(p, df=df, ncol=ncol, scores.title=scores.title, 
        size=size, dots.on.top=dots.on.top, fill="grey")

    p
}
