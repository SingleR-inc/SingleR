#' Plot delta distributions 
#'
#' Plot the distribution of deltas 
#' (i.e., the gap between the assignment score for the assigned label and those of the remaining labels)
#' across cells assigned to each reference label.
#'
#' @inheritParams plotScoreDistribution
#' @param chosen.only Logical scalar indicating whether to only show deltas 
#' for individual labels that were chosen as the final label in a combined result.
#' @param show String specifying whether to show the difference from the median (\code{"delta.med"}) 
#' or the difference from the next-best score (\code{"delta.next"}).
#'
#' @return
#' If \code{references} specifies a single set of deltas,
#' a \link[ggplot2]{ggplot} object is returned showing the deltas in violin plots.
#'
#' If \code{references} specifies multiple deltas for a combined result,
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
#' this function will show the deltas faceted by the assigned label within each individual reference.
#' The references to show in this manner can be specified using the \code{references} argument, 
#' entries of which refer to columns of \code{results$orig.results}.
#'
#' By default, a separate plot is created for each individual reference in a combined \code{results}.
#' Deltas are only shown in each plot if the label in the corresponding reference 
#' was chosen as the overall best label in the combined results.
#' However, this can be changed to show all deltas for an individual reference by setting \code{chosen.only=FALSE}.
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
#'
#' plotDeltaDistribution(results = combined)
#'
#' plotDeltaDistribution(results = combined, chosen.only=FALSE)
#'
#' # Tweaking the grid controls:
#' plotDeltaDistribution(combined, grid.vars = list(ncol = 2))
#'
#' @export
plotDeltaDistribution <- function(
    results,
    show = c("delta.med", "delta.next"),
    labels.use = colnames(results$scores),
    references = NULL,
    chosen.only = TRUE,
    size = 2,
    ncol = 5,
    dots.on.top = TRUE,
    this.color = "#000000",
    pruned.color = "#E69F00",
    grid.vars = list())
{
    results <- .ensure_named(results)
    show <- match.arg(show)
    is.combined <- !is.null(results$orig.results)
    ref.names <- colnames(results$orig.results)

    if (is.null(references)) {
        if (is.combined) {
            # Combined 'delta.med'/'delta.next' are not worth showing.
            references <- seq_along(results$orig.results) 
        } else {
            references <- 0L
        }
    }

    plots <- vector("list", length(references))
    for (i in seq_along(plots)) {

        # Pulling out the scores to use in this iteration.
        chosen <- references[i]
        if (is.combined) {
            if (chosen==0L) {
                stop("deltas cannot be shown for combined results")
            }
            current.results <- results$orig.results[[chosen]]
            if (chosen.only) {
                current.results <- current.results[chosen==results$reference,]
            }
        } else {
            current.results <- results
        }
        scores.title <- .values_title(is.combined, chosen, ref.names, show)

        # Computing the values that we want to show.
        if (show=="delta.med") {
            values <- getDeltaFromMedian(current.results)
        } else if (show=="delta.next") {
            tuned <- current.results$tuning.scores
            values <- tuned$first - tuned$second
        }

        # Pulling out the labels to use in this iteration.
        labels <- current.results$labels
        labels.title <- .values_title(is.combined, chosen, ref.names, "Labels")

        # Checking if we need pruning.
        pruned <- NULL
        if (!is.na(pruned.color)) {
            pruned <- is.na(current.results$pruned.labels)
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
        if (!is.null(grid.vars) && length(references) > 1L) {
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

    aes.jit <- NULL
    if (!is.null(pruned)) {
        df$pruned <- pruned
        aes.jit <- ggplot2::aes_string(color = "pruned")
    }

    # Trim dataframe by labels:
    if (any(keep <- labels.use %in% df$label)) {
        labels.use <- labels.use[keep]
        df <- df[df$label %in% labels.use,]
    } else {
        warning("ignoring 'labels.use' as it has no values in ", scores.title)
    }

    # Making the violin plots.
    p <- ggplot2::ggplot(data = df, ggplot2::aes_string(x="x", y="values")) + 
        ggplot2::xlab("")

    if (!is.null(pruned)) {
        p <- p + ggplot2::scale_color_manual(
            name = "Pruned",
            breaks = c("FALSE", "TRUE"),
            values = c(this.color, pruned.color))
    }

    jit <- ggplot2::geom_jitter(mapping = aes.jit, height = 0, width = 0.3, 
        shape = 16, size = size, na.rm = TRUE)

    .pretty_violins(p, df=df, ncol=ncol, scores.title=scores.title, 
        size=size, dots.on.top=dots.on.top, jitter=jit, fill="grey")
}
