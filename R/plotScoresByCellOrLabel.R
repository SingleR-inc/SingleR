#' Plot SingleR scoring on a per label or per cell basis.
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param cell.id Integer specifying which cell to show for \code{plotScoresSingleCell}
#' @param labels.use String vector indicating what labels to show in \code{plotScoresSingleCell} and \code{plotScoresMultiLabel}
#' If \code{labels.use} is left \code{NULL}, all labels available in \code{results} are presented.
#' @param label String indicating which individual label to plot in \code{plotScoresSingleLabel}
#' @param prune.calls Logical vector, the output of \code{\link{pruneScores}}.
#' This input will be unnecessary once \code{\link{pruneScores}}'s output is added to the results DataFrame
#' @param dots.on.top Logical which sets whether cell dots are plotted on top of, versus behind, the violin plots in \code{plotScoresSingleLabel} and \code{plotScoresMultiLabel}
#' @param colors String vector that sets the colors.  Order of colors should be: `this label`, `this label - pruned`, `other label`, `other label - pruned`.
#' @param size the size of the dots
#' @name plotScoresByCellOrLabel
#' @return Each function returns a \link{ggplot} object showing SingleR scores in a dot and/or violin plot representation.
#' 
#' @details
#' NOTE: these functions show initial scores only.
#' Aside from use of post-fine-tuning labels, they are completely independent of the fine-tuning step.
#' 
#' The \code{plotScoresSingleCell} function creates a dot plot showing the scores of a single cell accross many labels,
#' with a dotted line added to indicate the median score accross all labels.
#' It can be used to assess how an individual cell scored versus all/many labels of the reference dataset,
#' and may be useful for visualizing and tuning the \code{min.diff.med} per-cell cutoff of the \code{\link{pruneScores}} function.
#' 
#' The \code{plotScoresSingleLabel} and \code{plotScoresMultiLabel} functions creates jitter and violin plots showing
#' the scores of all cells accross either a single label or multiple labels, respectively.
#' These functions can be used to assess the distribution of scores of all cells for individual labels,
#' and may be useful for visualizing and tuning the \code{nmads} per-label cutoff of the \code{\link{pruneScores}} function.
#' 
#' Scores are grouped and colored by whether they were the final calls for a cell or not.
#' If pruneScores has been run and prune.calls are provided, scores are also separated and colored based on whether cells
#' were pruned versus not.
#' 
#' @author Daniel Bunis
#' @examples
#' example(SingleR, echo=FALSE)
#' prune <- pruneScores(pred)
#' 
#' plotScoresSingleCell(results = pred, prune.calls = prune, cell.id = 1)
#' plotScoresSingleLabel(results = pred, prune.calls = prune, label = "B",
#'     dots.on.top = TRUE)
#' plotScoresMultiLabels(results = pred,
#'     dots.on.top = TRUE, size = 0.5)
#' 
NULL

#' @describeIn plotScoresByCellOrLabel Plot scores accross labels of an individual cells
#' @export
#' @importFrom stats median
plotScoresSingleCell <- function(results, cell.id, prune.calls = NULL,
    labels.use = levels(as.factor(results$labels)), size = 2,
    colors = c("#F0E442", "#56B4E9", "gray70", "gray40")) {

    if (length(colors)<4) {
        stop("4 colors are expected. Order = 'this label', 'this label - pruned',
            'other label', 'other label - pruned'")
    }
    # Name the colors
    if (is.null(names(colors))) {
        names(colors) <- 
            c('this label', 'this label - pruned',
            'other label', 'other label - pruned')
    }

    # Add rownames to the results, which will be used for trimming scores data
    #   to the target cell later on
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }

    # Get the scores data for all cells
    df <- .data_gather(results, prune.calls)
    # Trim to just the data for the target cell
    df <- df[df$id == rownames(results)[cell.id],]
    # Calculate the cell's median score based on all labels
    scores.median <- median(df$score)
    # Trim to just the data for the target labels
    df <- df[df$label %in% labels.use,]

    # Make the plot
    p <- ggplot(
            data = df,
            aes(x = label, y = score, fill = called.this)) +
        theme_classic() +
        # Set labels names to be rotated 60 degrees.
        theme(axis.text.x= element_text(
            angle=60, hjust = 1, vjust = 1, size=12)) +
        # Remove "label" label from the x-axis
        xlab(NULL) +
        # Add median score line
        geom_hline(yintercept = scores.median, color = "gray",
            linetype = "dashed") +
        # Add scores data points
        geom_point(color = "black", shape = 21, size = size, alpha = 1) +
        # Set the colors
        scale_fill_manual(values = colors)

    p
}

#' @describeIn plotScoresByCellOrLabel Plot scores accross labels of an individual cells
#' @export
plotScoresSingleLabel <- function(results, prune.calls = NULL, label, size = 0.5, dots.on.top = FALSE, df = NULL,
    colors = c("#F0E442", "#56B4E9", "gray70", "gray40")){

    if (length(colors)<4) {
        stop("4 colors are expected. Order = 'this label', 'this label - pruned',
            'other label', 'other label - pruned'")
    }
    # Name the colors
    if (is.null(names(colors))) {
        names(colors) <- 
            c('this label', 'this label - pruned',
            'other label', 'other label - pruned')
    }

    # Add rownames to the results, which will be used for trimming scores data
    #   to the target cell later on
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }

    # Get the scores data
    df <- .data_gather(results, prune.calls, label)
    # Trim to the target label
    df <- df[df$label == label,]
    
    # Make the plot
    p <- ggplot(
            data = df,
            aes(x = called.this, y = score, fill = called.this)) + 
        theme_classic() +
        # Set the colors
        scale_fill_manual(values = colors) + 
        # Remove x-axis labels for the groupings (already in the legend),
        #   but do show the name of the target `label` as the axis title.
        scale_x_discrete(name = label, labels = NULL)
    # Add the data as jitter and violin.
    #   Violin plot added first iff `dots.on.top = TRUE`
    if (dots.on.top) {
        p <- p+ geom_violin()
    }
    p <- p + geom_jitter(
        height = 0, width = 0.3, color = "black", shape = 16,size = size)
    if (!dots.on.top) {
        p <- p + geom_violin()
    }
    
    p
}

#' @describeIn plotScoresByCellOrLabel Plot scores accross labels of an individual cells
#' @export
plotScoresMultiLabels <- function(results, prune.calls = NULL, size = 0.2, dots.on.top = FALSE,
    labels.use = levels(as.factor(results$labels)), ncol = 5,
    colors = c("#F0E442", "#56B4E9", "gray70", "gray40"), ...){
    if(length(colors)<4){stop("4 colors are expected.")}
    if (is.null(names(colors))){names(colors) <- c('this label', 'this label - pruned', 'other label', 'other label - pruned')}
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    df <- .data_gather(results, prune.calls, labels.use)
    max <- max(df$score)
    plots <- lapply(labels.use, function(X) {
        plotScoresSingleLabel(results, prune.calls, label = X, size, dots.on.top, df = df, colors) +
            theme(legend.position = "none", axis.ticks.x=element_blank()) +
            coord_cartesian(ylim = c(0,max)) + ylab(NULL)
    })
    plots <- c(plots,
               list(cowplot::ggdraw(cowplot::get_legend(
        plotScoresSingleLabel(results, prune.calls = prune.calls, label = labels.use[1])))))
    gridExtra::grid.arrange(grobs=plots, ncol = ncol, ...)
}

.data_gather <- function(results, prune.calls = NULL, labels.use = levels(as.factor(results$labels)))
{
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    labels.use <- labels.use[labels.use %in% colnames(results$scores)]
    scores <- results$scores[,colnames(results$scores) %in% labels.use]

    df <- data.frame(
        id = c(sapply(rownames(results), function(X) rep(X, length(labels.use)))),
        called = c(sapply(results$labels, function(X) rep(X, length(labels.use)))),
        label = rep(colnames(results$scores)[colnames(results$scores) %in% labels.use], nrow(results)),
        score = as.numeric(t(scores)),
        stringsAsFactors = FALSE)
    df$called.this <- "other label"
    df$called.this[df$label == df$called] <- "this label"
    if (!is.null(prune.calls)){
        prune.string <- as.character(factor(prune.calls, labels = c(""," - pruned")))
        df$called.this <- paste0(df$called.this,
                                 c(sapply(prune.string, function(X) rep(X, length(labels.use)))))
        df$called.this <- factor(df$called.this, levels = c('this label', 'this label - pruned', 'other label', 'other label - pruned'))
    }
    df
}
