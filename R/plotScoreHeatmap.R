#' Plot a score heatmap
#'
#' Create a heatmap of the \code{\link{SingleR}} assignment scores across all cell-label combinations.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}},
#' \code{\link{classifySingleR}}, \code{\link{combineCommonResults}}, or \code{\link{combineRecomputedResults}}.
#' @param cells.use Integer or string vector specifying the single cells (i.e., rows of \code{results}) to show.
#' If \code{NULL}, all cells are shown.
#' @param labels.use Character vector specifying the labels to show in the heatmap rows.
#' Defaults to all labels in \code{results}.
#' @param scores.use Integer scalar or vector specifying the individual annotation result from which to take scores.
#' This is only relevant for combined results, see Details.
#' @param calls.use Integer scalar or vector specifying the individual annotation result from which to take labels,
#' for use in the annotation bar when \code{show.labels=TRUE}.
#' This is only relevant for combined results, see Details.
#' @param clusters String vector or factor containing cell cluster assignments, to be shown as an annotation bar in the heatmap.
#' @param show.labels Logical indicating whether the assigned labels should be shown as an annotation bar.
#' @param show.pruned Logical indicating whether the pruning status of the cell labels,
#' as defined by \code{\link{pruneScores}}, should be shown as an annotation bar.
#' @param max.labels Integer scalar specifying the maximum number of labels to show.
#' @param normalize Logical specifying whether correlations should be normalized to lie in [0, 1].
#' @param order.by String providing the annotation to be used for cells/columns ordering.
#' Can be "labels" (default) or "clusters" (when provided).
#' Ignored if \code{cells.order} or \code{cluster_cols} are specified.
#' @param cells.order Integer or String vector specifying how to order the cells/columns of the heatmap.
#' Regardless of \code{cells.use}, this input should be the the same length as the total number of cells.
#' Ignored if \code{cluster_cols} is set.
#' @param rows.order String vector specifying how to order rows of the heatmap.
#' Contents should be the reference-labels in the order you would like them to appear, from top-to-bottom.
#' For combined results, include labels for all plots in a single vector and labels relevant to each plot will be extracted.
#' @param na.color String specifying the color for non-calculated scores of combined \code{results}.
#' @param annotation_col,cluster_cols,show_colnames,color,silent,...
#' Additional parameters for heatmap control passed to \code{\link[pheatmap]{pheatmap}}.
#' @param grid.vars A named list of extra variables to pass to \code{\link[gridExtra]{grid.arrange}},
#' used to arrange the multiple plots generated when \code{scores.use} is of length greater than 1.
#'
#' @return
#' If \code{scores.use} specifies a single set of scores,
#' the output of \code{\link[pheatmap]{pheatmap}} is returned showing the heatmap on the current graphics device.
#'
#' If \code{scores.use} specifies multiple scores for a combined result,
#' multiple heatmaps are generated in a grid on the current graphics device.
#'
#' If \code{scores.use} specifies multiple scores and \code{grid.vars} is set to \code{NULL},
#' a list is returned containing the \code{\link[pheatmap]{pheatmap}} globs for manual display.
#'
#' @details
#' This function creates a heatmap containing the \code{\link{SingleR}} initial assignment scores
#' for each cell (columns) to each reference label (rows).
#' Users can then easily identify the high-scoring labels associated with each cell and/or cluster of cells.
#'
#' If \code{show.labels=TRUE}, an annotation bar will be added to the heatmap showing the label assigned to each cell.
#' This is also used to order the columns for a more organized visualization when \code{order.by="label"}.
#' Note that scores shown in the heatmap are initial scores prior to the fine-tuning step,
#' so the reported labels may not match up to the visual maximum for each cell in the heatmap.
#'
#' If \code{max.labels} is less than the total number of unique labels, only the top labels are shown in the plot.
#' Labels that were called most frequently are prioritized.
#' The remaining labels are then selected based on:
#' \itemize{
#' \item Labels with max z-scores after per-cell centering and scaling of the scores matrix,
#' if \code{results} does not contain combined scores.
#' \item Labels which were suggested most frequently by individual references,
#' if \code{results} contains combined scores.
#' }
#'
#' @section Working with combined results:
#' For combined results (see \code{?\link{combineRecomputedResults}}),
#' this function can show both the combined and individual scores or labels.
#' This is done using the \code{scores.use} and \code{calls.use} arguments,
#' entries of which refer to columns of \code{results$orig.results} if positive or to the combined results if zero.
#' For example:
#' \itemize{
#' \item If we set \code{scores.use=2} and \code{calls.use=1},
#' we will plot the scores from the second individual reference
#' with the annotation bar containing label assignments from the first reference.
#' \item If we set \code{scores.use=1:2} and \code{calls.use=1:2},
#' we will plot the scores from first and second references (in separate plots)
#' with the annotation bar in each plot containing the corresponding label assignments.
#' \item By default, the function will create a separate plot the combined scores and each individual reference.
#' In each plot, the annotation bar contains the combined labels;
#' this is equivalent to \code{scores.use=0:N} and \code{calls.use=0} for \code{N} individual references.
#' }
#'
#' @section Tweaking the output:
#' Additional arguments can be passed to \code{\link[pheatmap]{pheatmap}} for further tweaking of the heatmap.
#' Particularly useful parameters are \code{show_colnames}, which can be used to display cell/cluster names;
#' \code{treeheight_row}, which sets the width of the clustering tree;
#' and \code{annotation_col}, which can be used to add extra annotation layers.
#' Clustering, pruning and label annotations are automatically generated and appended to \code{annotation_col} when available.
#'
#' @section Normalization of colors:
#' If \code{normalize=TRUE}, scores will be linearly adjusted for each cell
#' so that the smallest score is 0 and the largest score is 1.
#' This is followed by cubing of the adjusted scores to improve dynamic range near 1.
#' Visually, the color scheme is changed to a blue-green-yellow scale.
#'
#' The adjustment is intended to inflate differences between scores within a given cell for easier visualization.
#' This is because the scores are often systematically shifted between cells,
#' making the raw values difficult to directly compare.
#' However, it may be somewhat misleading;
#' fine-tuning may appear to assign a cell to a label with much lower score whereas the actual scores are much closer.
#' It is for this reason that the color bar values are not shown as the absolute values of the score have little meaning.
#'
#' Note that this transformation is not dependent on the choice of the top \code{max.labels} labels.
#' Altering \code{max.labels} will not change the normalized values, only the labels that are shown.
#' However, the transformation will respond to \code{labels.use}.
#'
#' @seealso
#' \code{\link{SingleR}}, to generate \code{scores}.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores.
#'
#' \code{\link[pheatmap]{pheatmap}}, for additional tweaks to the heatmap.
#'
#' \code{\link[gridExtra]{grid.arrange}}, for tweaks to the how heatmaps are arranged when multiple are output together.
#'
#' @author Daniel Bunis, based on code by Dvir Aran.
#'
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' # Grab the original identities of the cells as mock clusters
#' clusts <- test$label
#'
#' # Creating a heatmap with just the labels.
#' plotScoreHeatmap(pred)
#'
#' # Creating a heatmap with clusters also displayed.
#' plotScoreHeatmap(pred,
#'     clusters=clusts)
#'
#' # Creating a heatmap with whether cells were pruned displayed.
#' plotScoreHeatmap(pred,
#'     show.pruned = TRUE)
#'
#' # We can also turn off the normalization with Normalize = FALSE
#' plotScoreHeatmap(pred, clusters=clusts,
#'     normalize = FALSE)
#'
#' # To only show certain labels, you can use labels.use or max.labels
#' plotScoreHeatmap(pred, clusters=clusts,
#'     labels.use = c("A","B","D"))
#' plotScoreHeatmap(pred, clusters=clusts,
#'     max.labels = 4)
#'
#' # We can pass extra tweaks the heatmap as well
#' plotScoreHeatmap(pred, clusters=clusts,
#'     fontsize_row = 20)
#' plotScoreHeatmap(pred, clusters=clusts,
#'     treeheight_row = 15)
#' plotScoreHeatmap(pred, clusters=clusts, cluster_col = TRUE,
#'     cutree_cols = 5)
#'
#' ### Multi-Reference Compatibility ###
#'
#' example(combineRecomputedResults, echo = FALSE)
#' plotScoreHeatmap(combined)
#'
#' # 'scores.use' sets which particular run's scores to show, and can be multiple
#' plotScoreHeatmap(combined,
#'     scores.use = 1)
#' plotScoreHeatmap(combined,
#'     scores.use = c(0,2))
#'
#' # 'calls.use' adjusts which run's labels and pruning calls to display.
#' plotScoreHeatmap(combined,
#'     calls.use = 1)
#'
#' # To have plots output in a grid rather than as separate pages, provide,
#' # a list of inputs for gridExtra::grid.arrange() to 'grids.vars'.
#' plotScoreHeatmap(combined,
#'     grid.vars = list(ncol = 1))
#'
#' # An empty list will use grid.arrange defaluts
#' plotScoreHeatmap(combined,
#'     grid.vars = list())
#'
#' @export
#' @importFrom utils head
#' @importFrom DelayedArray rowMaxs rowMins
plotScoreHeatmap <- function(results, cells.use = NULL, labels.use = NULL,
    clusters = NULL, show.labels = TRUE, show.pruned = FALSE,
    max.labels = 40, normalize = TRUE,
    cells.order = NULL, order.by = c("labels","clusters"), rows.order = NULL,
    scores.use = NULL, calls.use = 0, na.color = "gray30",
    cluster_cols = FALSE,
    annotation_col = NULL, show_colnames = FALSE,
    color = grDevices::colorRampPalette(c("#D1147E", "white", "#00A44B"))(100),
    silent = FALSE, ..., grid.vars = list())
{
    results <- .ensure_named(results)
    is.combined <- !is.null(results$orig.results)
    ref.names <- colnames(results$orig.results)

    if (is.null(scores.use)) {
        scores.use <- c(0L, seq_along(results$orig.results)) # seq_along(NULL) is nothing.
    }
    calls.use <- rep(calls.use, length.out=length(scores.use))

    # Delaying the plotting to a single grid.arrange call,
    # even if it's not requested to be silent.
    use.grid <- !is.null(grid.vars) && length(scores.use) > 1L

    plots <- vector("list", length(scores.use))
    for (i in seq_along(plots)) {

        # Pulling out the scores to use in this iteration.
        chosen.scores <- scores.use[i]
        if (chosen.scores==0L) {
            score.results <- results
        } else {
            score.results <- results$orig.results[[chosen.scores]]
        }

        scores <- score.results$scores
        rownames(scores) <- rownames(results)
        scores.title <- .values_title(is.combined, chosen.scores, ref.names, "Scores")
        scores.labels <- score.results$labels

        # Pulling out the labels to use in this iteration.
        chosen.calls <- calls.use[i]
        if (chosen.calls==0L) {
            call.results <- results
        } else {
            call.results <- results$orig.results[[chosen.calls]]
        }

        labels <- call.results$labels
        prune.calls <- call.results$pruned.labels
        names(labels) <- names(prune.calls) <- rownames(scores)
        labels.title <- .values_title(is.combined, chosen.calls, ref.names, "Labels")

        # Actually creating the heatmap.
        output <- .plot_score_heatmap(
            scores=scores,
            labels=labels,
            prune.calls=prune.calls,
            cells.use=cells.use,
            labels.use=labels.use,
            max.labels=max.labels,
            clusters=clusters,
            cells.order,
            order.by=order.by,
            rows.order=rows.order,
            show.labels=show.labels,
            show.pruned=show.pruned,
            scores.title=scores.title,
            labels.title=labels.title,
            show_colnames=show_colnames,
            cluster_cols=cluster_cols,
            annotation_col=annotation_col,
            silent=silent || use.grid,
            color=color,
            na.color=na.color,
            normalize=normalize,
            scores.labels=scores.labels,
            ...)

        if (use.grid) {
            plots[[i]] <- output[[4]]
        } else {
            plots[[i]] <- output
        }
    }

    if (length(plots)==1L) {
        # Doing this to be consistent with raw pheatmap() output.
        plots[[1]]
    } else {
        if (use.grid) {
            do.call(gridExtra::grid.arrange, c(plots, grid.vars))
        } else {
            plots
        }
    }
}

.plot_score_heatmap <- function(
    scores, labels, prune.calls,
    cells.use, labels.use, max.labels,
    clusters, cells.order, order.by, rows.order,
    show.labels, show.pruned,
    scores.title, labels.title,
    show_colnames, cluster_cols, annotation_col, silent,
    color, na.color, normalize, scores.labels, ...)
{
    # 'scores' is guaranteed to be named by this point.
    clusters <- .name_unless_NULL(clusters, rownames(scores))
    cells.order <- .name_unless_NULL(cells.order, rownames(scores))

    # Adjust data
    scores <- .trim_normalize_reorder_scores(
        scores=scores,
        scores.title=scores.title,
        labels.use=labels.use,
        max.labels=max.labels,
        cells.use=cells.use,
        normalize=normalize,
        cluster_cols=cluster_cols,
        order.by=order.by,
        cells.order=cells.order,
        rows.order=rows.order,
        labels=labels,
        clusters=clusters,
        scores.labels)

    # Compile annotations
    annotation_col <- .make_annotation_col(
        annotation_col=annotation_col,
        show.labels=show.labels,
        labels=labels,
        labels.title=labels.title,
        show.pruned=show.pruned,
        prune.calls=prune.calls,
        clusters=clusters)

    ### Create base args list for making the heatmap
    args <- list(border_color = NA, show_colnames = show_colnames,
        clustering_method = 'ward.D2', cluster_cols = cluster_cols,
        silent = silent, annotation_col = annotation_col,
        ...)

    if (is.null(args$cluster_rows)) {
        args$cluster_rows <- is.null(rows.order) && ncol(scores)>1
    }
    if (is.null(args$main)) {
        args$main <- scores.title
    }

    # Add annotation colors
    if (is.null(args$annotation_colors)) {
        args$annotation_colors <-
            .make_heatmap_annotation_colors(args, show.pruned)
    }

    # Add scores & score colors
    ## Set score colors and legend display
    if (normalize && ncol(scores) > 1) {
        color <- viridis::viridis(100)
        args$breaks <- seq(0, 1, length.out = 101)
        args$legend_breaks <- c(0,1)
        args$legend_labels <- c("Lower", "Higher")
    } else {
        abs.max <- max(abs(range(scores, na.rm = TRUE)))
        breaks.len <- length(color)+1
        args$breaks <- seq(-abs.max, abs.max, length.out = breaks.len)
        args$legend_breaks <- c(-abs.max, abs.max, length.out = 3)
        args$legend_labels <- round(args$legend_breaks, 3)
    }
    args$color <- color

    # Replace NAs and add na.color
    if (any(is.na(scores))) {
        # value should be 10% distance above current max
        NA_val <- max(args$breaks) + 0.1*diff(range(args$breaks))
        scores[is.na(scores)] <- NA_val
        args$color <- c(args$color, na.color)
        args$breaks <- c(args$breaks, NA_val)
        args$legend_breaks <- c(args$legend_breaks, NA_val)
        args$legend_labels <- c(args$legend_labels, "NA")
    }
    args$mat <- t(scores)

    # Troubleshooting and testing purposes
    if (!is.null(args$return.data) && args$return.data) {
        return(args)
    }

    do.call(pheatmap::pheatmap, args)
}

.make_annotation_col <- function(
    annotation_col = NULL, show.labels, labels, labels.title,
    show.pruned, prune.calls, clusters = NULL) {

    if (is.null(annotation_col)) {
        annotation_col <- data.frame(row.names = names(labels))
    }
    if (show.pruned && !is.null(prune.calls)) {
        # Pruned calls added this way to ensure they come first for coloring purposes.
        Pruned <- as.character(is.na(prune.calls)[rownames(annotation_col)])
        annotation_col <- cbind(Pruned,annotation_col)
    }
    if (show.labels) {
        annotation_col$Labels <- labels[rownames(annotation_col)]
        annot.titles <- colnames(annotation_col)
        annot.titles[annot.titles == "Labels"] <- labels.title
        names(annotation_col) <- annot.titles
    }
    if (!is.null(clusters)) {
        annotation_col$Clusters <- clusters[rownames(annotation_col)]
    }

    if (!ncol(annotation_col)>0) {
        return(NULL)
    }
    annotation_col
}

.trim_normalize_reorder_scores <- function(
    scores, scores.title,
    labels.use, max.labels, cells.use, normalize,
    cluster_cols, order.by, cells.order, rows.order,
    labels, clusters, scores.labels)
{
    scores <- .trim_byLabel_and_normalize_scores(
        scores, labels.use, max.labels, normalize, scores.title, scores.labels)

    if (!is.null(cells.use)) {
        # Trim by cell
        scores <- scores[cells.use,,drop=FALSE]

        # Trim potential ordering vars
        clusters <- clusters[cells.use]
        labels <- labels[cells.use]
        cells.order <- cells.order[cells.use]
    }

    if (!cluster_cols) {
        # Order: priority = 'cells.order', then 'order.by' which can be labels or clusters.
        scores <- .order_score_matrix_cells(
            scores, cluster_cols, order.by, cells.order, labels, clusters)
    }

    if (!is.null(rows.order)) {
        if (any(!colnames(scores) %in% rows.order)) {
            missing <- colnames(scores)[!colnames(scores) %in% rows.order]
            warning("Label(s) of ", scores.title, " missing from 'rows.order' will not be plotted: ",
                    paste0(missing, collapse = ", "))
        }
        scores <- scores[,rows.order[rows.order %in% colnames(scores)]]
    }

    scores
}

.trim_byLabel_and_normalize_scores <- function(
    scores, labels.use, max.labels, normalize, scores.title, scores.labels) {

    # Trim by labels (remove any with no scores)
    all.na <- apply(scores, 2, FUN = function(x) all(is.na(x)))
    scores <- scores[,!all.na, drop = FALSE]

    # Trim by labels (labels.use)
    if (!is.null(labels.use)) {
        labels.use <- labels.use[labels.use %in% colnames(scores)]
        if (length(labels.use)>0){
            scores <- scores[,labels.use,drop=FALSE]
        } else {
            warning("ignoring 'labels.use' without any values in ", scores.title)
        }
    }

    # Trim by labels (max.labels), using primarily the most frequent labels.
    times.best <- table(factor(scores.labels, levels = unique(colnames(scores))))[colnames(scores)]
    if (!any(is.na(scores))) {
        # To break ties, we sort by the scaled maximum if there are no NAs.
        # This is done _before_ within-cell normalization of the scores,
        # after which it makes little sense to compare scores between cells.
        secondary <- rowMaxs(scale(t(scores)), na.rm = TRUE)
    } else {
        # If there are NAs - usually from combineRecomputedResults -
        # we sort by the frequency of non-NA occurrences.
        secondary <- apply(scores, 2, FUN = function(x) sum(!is.na(x)))
    }
    to.keep <- order(times.best, secondary, decreasing=TRUE)
    to.keep <- head(to.keep, max.labels)

    # Normalize the scores between [0, 1] and cube to create more separation.
    if (normalize) {
        if (ncol(scores) > 1L) {
            mmax <- rowMaxs(scores, na.rm = TRUE)
            mmin <- rowMins(scores, na.rm = TRUE)
            scores <- (scores-mmin)/pmax(mmax-mmin, 1e-8)
            scores <- scores^3
        } else {
            warning("disabling normalization with only one label in ", scores.title)
        }
    }

    # Drop labels exceeding 'max.labels'.
    scores[,to.keep,drop=FALSE]
}

.order_score_matrix_cells <- function(
    scores, cluster_cols, order.by = c("labels","clusters"),
    cells.order, labels, clusters) {
    # By: cells.order, if provided, else by 'order.by' which = "labels" by default, or "clusters".

    if (!is.null(cells.order)) {
        order <- order(cells.order)
    } else {
        order.stat <- switch(match.arg(order.by),
            labels=labels,
            clusters=clusters
        )
        if (is.null(order.stat)) {
            stop("'clusters' input is required when 'order.by=\"clusters\"'")
        }
        order <- order(order.stat)
    }

    scores[order,,drop=FALSE]
}

.make_heatmap_annotation_colors <- function(args, show.pruned) {
    # Create pheatmap annotations_colors dataframe
        # list of character vectors, all named.
            # vector names = annotation titles
            # vector members' (colors') names = annotation identities

    # Extract a default color-set
    annotation.colors.d <- .make_heatmap_colors_discrete(show.pruned)
    annotation.colors.n <- .make_heatmap_colors_numeric()

    # Initiate variables
    next.color.index.discrete <- 1
    next.color.index.numeric <- 1
    col_colors <- NULL
    row_colors <- NULL

    # Columns First (if there)
    if (!is.null(args$annotation_col)) {
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_col,
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        col_colors <- dfcolors_out$df_colors
        next.color.index.discrete <- dfcolors_out$next.color.index.discrete
        next.color.index.numeric <- dfcolors_out$next.color.index.numeric
    }

    # Rows Second (if there)
    if (!is.null(args$annotation_row)) {
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_row,
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        row_colors <- dfcolors_out$df_colors
    }

    c(col_colors, row_colors)
}

.make_heatmap_colors_discrete <- function(show.pruned) {
    # Creates a default vector of colors with 40*10 (overkill) options.
    annotation.colors <- rep(
        # DittoSeq-v1.4 Colors (based on Okabe-Ito colors)
        c(
            "#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#666666",
            "#AD7700", "#1C91D4", "#007756", "#D5C711",
            "#005685", "#A04700", "#B14380", "#4D4D4D",
            "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
            "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
            "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
            "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3",
            "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
            "#00446B", "#803800", "#8D3666", "#3D3D3D"),
        10)
    if (show.pruned) {
        annotation.colors <- c("white", annotation.colors)
    }
    annotation.colors
}

.make_heatmap_colors_numeric <- function() {
    # Creates a default vector of colors with 8*3 (overkill) options.
    # These represent max.colors for discrete color scales.
    rep(
        # DittoSeq-v0.2.10 Colors, distinct order, (based on Okabe-Ito colors)
        c(
            "#B14380", "#A04700", "#005685", "#D5C711", "#007756",
            "#1C91D4", "#AD7700", "#4D4D4D", "#CC79A7", "#D55E00",
            "#0072B2", "#F0E442", "#009E73", "#56B4E9", "#E69F00",
            "#666666"),
        3)
}

# Interpret annotations dataframe,
# Pick, name, and add colors.
.pick_colors_for_df <- function(
    annotation_df,
    next.color.index.discrete, next.color.index.numeric,
    annotation.colors.d, annotation.colors.n
    ) {
    df_colors <- NULL
    for (i in seq_len(ncol(annotation_df))){

        # Determine the distinct contents of the first annotation
        in.this.annot <- levels(as.factor(annotation_df[,i]))

        # Make new colors
        if(!is.numeric(annotation_df[,i])){
            # Take colors for each, and name them.
            new.colors <- annotation.colors.d[
                seq_along(in.this.annot) + next.color.index.discrete - 1
                ]
            names(new.colors) <- in.this.annot

            next.color.index.discrete <-
                next.color.index.discrete + length(in.this.annot)
        } else {
            # Make a 100 color split as in pheatmap code.
            a <- cut(
                annotation_df[order(annotation_df[,i]),i],
                breaks = 100)
            # Assign to colors.
            this.ramp <- annotation.colors.n[next.color.index.numeric]
            new.colors <-
                grDevices::colorRampPalette(c("white",this.ramp))(100)[a]

            next.color.index.numeric <- next.color.index.numeric + 1
        }
        # Add the new colors as the list
        df_colors <- c(
            df_colors,
            list(new.colors))
    }
    names(df_colors) <- names(annotation_df)
    list(df_colors = df_colors,
        next.color.index.discrete = next.color.index.discrete,
        next.color.index.numeric = next.color.index.numeric)
}
