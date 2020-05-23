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
#' @param show.nmads Numeric scalar that shows the threshold that would be used for pruning with \code{\link{pruneScores}}.
#' Only used when \code{show="delta.med"}.
#' @param show.min.diff Numeric scalar that shows the threshold that would be used for pruning with \code{\link{pruneScores}}.
#' Only used when \code{show="delta.med"} or \code{"delta.next"}.
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
#' Note that \code{\link{pruneScores}} trims based on the \code{min.diff.med} and \code{min.diff.next} cutoffs first,
#' before calculating the first-labels' delta medians.
#' Thus, the actual \code{nmads} cut-off used in \code{\link{pruneScores}} may vary from the one portrayed in the plot.
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
#' # Adjusting the default 'nmads' for the pruning threshold:
#' plotDeltaDistribution(pred, show.nmads = 2)
#'
#' # Alternatively using a 'min.diff' cutoff:
#' plotDeltaDistribution(pred, show.min.diff = 0.03)
#'
#' plotDeltaDistribution(pred, show = "delta.next", show.min.diff = 0.03)
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
    pruned.use = 0,
    size = 0.5,
    ncol = 5,
    dots.on.top = TRUE,
    show.nmads = 3,
    show.min.diff = NULL,
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
    pruned.use <- rep(pruned.use, length.out=length(deltas.use))

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

        # Calculate nmad.cutoff values for 'show.nmads'
        if (show == "delta.med" && !is.null(show.nmads)) {
            nmad.vals <- pruneScores(score.results, get.thresholds=TRUE, nmads = show.nmads)
        }

        # Actually creating the plot
        plots[[i]] <- .plot_delta_distribution(
            values=values, labels=labels, labels.use=labels.use,
            labels.title=labels.title, scores.title=scores.title, 
            this.color=this.color, size=size, ncol=ncol, dots.on.top=dots.on.top,
            show=show, show.nmads=show.nmads, show.min.diff=show.min.diff, nmad.vals=nmad.vals)
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
    values, labels, labels.use,
    labels.title, scores.title, 
    this.color, size, ncol, dots.on.top,
    show, show.nmads, show.min.diff, nmad.vals) 
{
    df <- data.frame(values=values, label=labels, x=character(length(values)))

    # Trim dataframe by labels:
    if (any(keep <- labels.use %in% df$label)) {
        labels.use <- labels.use[keep]
        df <- df[df$label %in% labels.use,]
    } else {
        warning("ignoring 'labels.use' as it has no values in ", scores.title)
    }

    # Making the violin plots.
    p <- ggplot2::ggplot(data = df, ggplot2::aes_string(y = "values", x = "x")) +
        ggplot2::xlab("")
    p <- .pretty_violins(p, df=df, ncol=ncol, scores.title=scores.title, 
        size=size, dots.on.top=dots.on.top)

    # Add cutoff lines:
    if (!is.null(show.nmads) || !(is.null(show.min.diff))) {
        p <- .add_cutoff_lines(p, show=show, show.nmads=show.nmads, 
            show.min.diff=show.min.diff, nmad.vals=nmad.vals, labels.use=labels.use)
    }

    p
}

.add_cutoff_lines <- function(p, show, show.nmads, show.min.diff, nmad.vals, labels.use) {
    # Add cutoff line for nmads cutoff
    if (show == "delta.med" && !is.null(show.nmads)) {
        p <- .add_nmads_lines(p, nmad.vals, labels.use)
    }

    # Add cutoff line for min.diff cutoff
    if (!is.null(show.min.diff)) {
        df_min.diff <- data.frame(color = "2.min.diff", show.min.diff = show.min.diff)
        p <- p + ggplot2::geom_hline(
            data = df_min.diff, na.rm = TRUE, size = 1.1,
            ggplot2::aes_string(yintercept = "show.min.diff", color = "color"))
    }

    # Set the colors and a=labels for combined legend.
    p + ggplot2::scale_color_manual(
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
}

.add_nmads_lines <- function(p, nmad.vals, labels.use) {
    df <- data.frame(
        label=names(nmad.vals),
        nmads.cutoff=nmad.vals,
        x = character(length(nmad.vals)),
        values=numeric(length(nmad.vals)), # to keep the inherited aes_string() happy.
        mad=rep("1.nmad", length(nmad.vals)) # Add dummy var for setting color later.
    ) 

    df <- df[df$label %in% labels.use, , drop=FALSE]

    p + ggplot2::geom_errorbar(data = df, na.rm=TRUE,
            ggplot2::aes_string(ymin= "nmads.cutoff", ymax= "nmads.cutoff", color = "mad"),
            width = 1, size = 1.1, show.legend = c(color = TRUE, fill = FALSE))
}
