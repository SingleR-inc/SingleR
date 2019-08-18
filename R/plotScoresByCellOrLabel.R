#' Plot SingleR scoring on a per label or per cell basis.
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param cell.id Integer specifying which cell to show for \code{plotScoresSingleCell}
#' @param labels.use String vector indicating what labels to show in \code{plotScoresSingleCell} and \code{plotScoresMultiLabel}
#' If \code{labels.use} is left \code{NULL}, all labels available in \code{results} are presented.
#' @param label String indicating which individual label to plot in \code{plotScoresSingleLabel}
#' @param prune.calls Logical vector, the output of \code{\link{pruneScores}}.
#' This input will be unnecessary once \code{\link{pruneScores}}'s output is added to the SingleR.results dataframe
#' @param dots.on.top Logical which sets whether cell dots are plotted on top of, versus behind, the violin plots in \code{plotScoresSingleLabel} and \code{plotScoresMultiLabel}
#' @param colors String vector that sets the colors.  Order of colors should be: `this label`, `this label - pruned`, `other label`, `other label - pruned`.
#' @param size the size of the dots
#' @name plotScoresByCellOrLabel
#' @return Dot and violin plot of SingleR scores separated by whether this was the final call or not.
#' If prune.calls are provided, will also separate and color based on whether cells' final calls were pruned versus not. 
#'
#' @details
#' The \code{plotScoresSingleCell} function creates a jitter plot showing the scores of a single cell accross many labels.
#' 
#' The \code{plotScoresSingleLabel} function creates a jitter and violin plot showing the scores of all cells accross a single label.
#' 
#' The \code{plotScoresMultiLabel} function creates \code{plotScoresSingleLabel} plots for multiple individual labels. 
#' 
#' @author Daniel Bunis
#' @examples
#' example(SingleR, echo=FALSE)
#' prune <- pruneScores(pred$scores)
#' 
#' plotScoresSingleCell(results = pred, prune.calls = prune, cell.id = 1)
#' plotScoresSingleLabel(results = pred, prune.calls = prune, label = "B",
#'     dots.on.top = TRUE)
#' plotScoresMultiLabels(results = pred,
#'     dots.on.top = TRUE, size = 0.5)
#' 
#' NULL
NULL


#' @describeIn plotScoresByCellOrLabel Plot scores accross labels of an individual cells
#' @export
plotScoresSingleCell <- function(results, cell.id, prune.calls = NULL,
    labels.use = levels(as.factor(results$labels)), size = 2,
    colors = c("#F0E442", "#56B4E9", "gray70", "gray40")){
    if(length(colors)<4){stop("4 colors are expected.")}
    if (is.null(names(colors))){names(colors) <- c('this label', 'this label - pruned', 'other label', 'other label - pruned')}
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    df <- plotScores.dataGather(results, prune.calls, labels.use)
    p <- ggplot(data = df[df$id == rownames(results)[cell.id],],
                aes(x = label, y = score, fill = called.this)) + 
        theme_classic() + geom_point(
                color = "black",
                shape = 21,
                size = size,
                alpha = 1) + scale_fill_manual(values = colors) +
        theme(axis.text.x= element_text(angle=60, hjust = 1, vjust = 1, size=12)) +
        xlab(NULL)
    p
}

#' @describeIn plotScoresByCellOrLabel Plot scores accross labels of an individual cells
#' @export
plotScoresSingleLabel <- function(results, prune.calls = NULL, label, size = 0.5, dots.on.top = FALSE, df = NULL,
    colors = c("#F0E442", "#56B4E9", "gray70", "gray40")){
    if(length(colors)<4){stop("4 colors are expected.")}
    if (is.null(names(colors))){names(colors) <- c('this label', 'this label - pruned', 'other label', 'other label - pruned')}
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    if (is.null(df)) {df <- plotScores.dataGather(results, prune.calls, label)}
    p <- ggplot(data = df[df$label == label,],
                aes(x = called.this, y = score, fill = called.this)) + 
        theme_classic() +
        scale_fill_manual(values = colors) + 
        scale_x_discrete(name = label, labels = NULL)
    if(dots.on.top){ p <- p+ geom_violin()}
    p <- p + geom_jitter(height = 0, width = 0.3, color = "black", shape = 16, size = size)
    if(!dots.on.top){ p <- p + geom_violin()}
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
    df <- plotScores.dataGather(results, prune.calls, labels.use)
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

plotScores.dataGather <- function(results, prune.calls = NULL, labels.use = levels(as.factor(results$labels)))
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