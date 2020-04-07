#' Plot a score heatmap
#'
#' Create a heatmap of the \code{\link{SingleR}} assignment scores across all cell-label combinations.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param cells.use Integer or string vector specifying the single cells to show.
#' If \code{NULL}, all cells are presented.
#' @param labels.use String vector indicating what labels to show.
#' If \code{NULL}, all labels available in \code{results} are presented.
#' @param scores.use,calls.use Integer which sets which scores or final labels and pruning calls to use when \code{results} is the output of a multiple reference \code{\link{SingleR}} run or of the \code{\link{combineCommonResults}} or \code{\link{combineRecomputedResults}} functions.
#' A value of 0 points to the overall results, while any other integer indicates the index of the individual output that should be targetted.
#' @param clusters String vector or factor containing cell cluster assignments, to be shown as an annotation bar in the heatmap.
#' @param show.labels Logical indicating whether the final labels of cells should be shown as an annotation bar.
#' @param show.pruned Logical indicating whether the pruning status of the labels should be shown as an annotation bar, as defined by \code{\link{pruneScores}}.
#' @param max.labels Integer scalar specifying the maximum number of labels to show.
#' @param normalize Logical specifying whether correlations should be normalized to lie in [0, 1].
#' @param order.by String providing the annotation to be used for clustering.
#' Can be "labels" (default) or "clusters" (when provided).
#' Subordinate to \code{cells.order} and \code{cluster_cols}.
#' @param cells.order Integer vector specifying the ordering of cells/columns of the heatmap.
#' Regardless of \code{cells.use}, this input should be the the same length as the total number of cells.
#' Subordinate to \code{cluster_cols}.
#' @param annotation_col,cluster_cols,show_colnames,color,... Additional parameters for heatmap control passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return A heatmap of assignment scores is generated on the current graphics device using \pkg{pheatmap}.
#'
#' @details
#' This function creates a heatmap containing the \code{\link{SingleR}} initial assignment scores for each cell (columns) to each reference label (rows).
#' Users can then easily identify the high-scoring labels associated with each cell and/or cluster of cells.
#'
#' If \code{show.labels=TRUE}, an annotation bar will be added to the heatmap indicating final labels assigned to the cells.
#' Note that scores shown in the heatmap are initial scores prior to the fine-tuning step, so the reported labels may not match up to the visual maximum for each cell in the heatmap.
#'
#' If \code{max.labels} is less than the total number of unique labels, only the labels with the largest maximum scores in \code{results} are shown in the plot.
#' Specifically, the set of scores for each cell is centred and scaled, and the maximum transformed score for each label is used to choose the labels to retain.
#'
#' Additional arguments can be passed to \code{\link[pheatmap]{pheatmap}} for further tweaking of the heatmap.
#' Particularly useful parameters are \code{show_colnames}, which can be used to display cell/cluster names;
#' and \code{annotation_col}, which can be used to add extra annotation layers.
#' Clustering, pruning and label annotations are automatically generated and appended to \code{annotation_col} when available.
#'
#' @section Normalization of colors:
#' If \code{normalize=TRUE}, scores will be linearly adjusted for each cell so that the smallest score is 0 and the largest score is 1.
#' This is followed by cubing of the adjusted scores to improve dynamic range near 1.
#' Visually, the color scheme is changed to a blue-green-yellow scale.
#'
#' The adjustment is intended to inflate differences between scores within a given cell for easier visualization.
#' This is because the scores are often systematically shifted between cells, making the raw values difficult to directly compare.
#' However, it may be somewhat misleading;
#' fine-tuning may appear to assign a cell to a label with much lower score whereas the actual scores are much closer.
#' It is for this reason that the color bar values are not shown as the absolute values of the score have little meaning.
#'
#' Also note that this transformation is done \emph{after} the choice of the top \code{max.labels} labels.
#'
#' @seealso
#' \code{\link{SingleR}}, to generate \code{scores}.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores.
#'
#' \code{\link[pheatmap]{pheatmap}}, for additional tweaks to the heatmap.
#'
#' @author Daniel Bunis, based on code by Dvir Aran.
#'
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#' # Grab the original identities of the cells as mock clusters
#' clusts <- g
#'
#' # Creating a heatmap with just the labels.
#' plotScoreHeatmap(pred)
#'
#' # Creating a heatmap with clusters also displayed.
#' plotScoreHeatmap(pred, clusters=clusts)
#'
#' # Creating a heatmap with whether cells were pruned displayed.
#' plotScoreHeatmap(pred, show.pruned = TRUE)
#'
#' # We can also turn off the normalization with Normalize = FALSE
#' plotScoreHeatmap(pred, clusters=clusts, normalize = FALSE)
#'
#' # To only show certain labels, you can use labels.use or max.labels
#' plotScoreHeatmap(pred, clusters=clusts, labels.use = c("A","B","D"))
#' plotScoreHeatmap(pred, clusters=clusts, max.labels = 4)
#'
#' # We can pass extra tweaks the heatmap as well
#' plotScoreHeatmap(pred, clusters=clusts, fontsize.row = 9)
#' plotScoreHeatmap(pred, clusters=clusts, cutree_col = 3)
#'
#' ### Multi-Reference Compatibility ###
#'
#' # When SingleR is run with multiple references, investigation of how scores
#' # and calls from particular references contributed to the final outcome can
#' # be achieved through utilization of scores.use and calls.use.
#' #
#' # NOTE: Final scoring in such runs are performed for a subset of labels for
#' # each cell, thus heatmaps cannot be made for final scores.
#' # Thus, default output shows final Labels, but targets scores of 1st reference
#' example(combineRecomputedResults, echo = FALSE)
#' plotScoreHeatmap(pred)
#'
#' # scores.use sets which run's scores to show
#' plotScoreHeatmap(pred,
#'     scores.use = 2)
#'
#' # calls.use sets which run's labels and pruning calls to display
#' plotScoreHeatmap(pred,
#'     scores.use = 1,
#'     calls.use = 1)
#'
#' @export
#' @importFrom utils head
#' @importFrom DelayedArray rowMaxs rowMins
plotScoreHeatmap <- function(results, cells.use = NULL, labels.use = NULL,
    scores.use = 0, calls.use = 0,
    clusters = NULL, show.labels = TRUE, show.pruned = FALSE,
    max.labels = 40, normalize = TRUE,
    cells.order = NULL, order.by = c("labels","clusters"), cluster_cols = FALSE,
    annotation_col = NULL, show_colnames = FALSE, color = NULL, ...)
{

    if (scores.use == 0 && !is.null(results$orig.results)) {
        message("Heatmap cannot be made from sparse final scores. Scores from run for 1st reference are shown instead.")
        scores.use <- 1
    }

    orig.res <- results
    results <- .ensure_named(.grab_results(orig.res, scores.use))

    if (is.null(annotation_col)) {
        annotation_col <- data.frame(row.names = rownames(results))
    }
    if (show.pruned) {
        prune.calls <- .grab_results(orig.res, calls.use)$pruned.labels
        if (!is.null(prune.calls)) {
            pruned <- as.character(is.na(prune.calls))
            names(pruned) <- rownames(results)
            annotation_col$Pruned <- pruned[rownames(annotation_col)]
        }
    }
    if (show.labels) {
        labels <- .grab_results(orig.res, calls.use)$labels
        names(labels) <- rownames(results)
        annotation_col$Labels <- labels[rownames(annotation_col)]
    }
    if (!is.null(clusters)) {
        names(clusters) <- rownames(results)
        annotation_col$Clusters <- clusters[rownames(annotation_col)]
    }

    scores <- results$scores
    rownames(scores) <- rownames(results)

    # Trim the scores and potential ordering vars by requested cells or labels
    labels <- results$labels
    if (!is.null(cells.use)) {
        scores <- scores[cells.use,,drop=FALSE]
        clusters <- clusters[cells.use]
        labels <- labels[cells.use]
        cells.order <- cells.order[cells.use]
    }
    if (!is.null(labels.use)) {
        scores <- scores[,labels.use,drop=FALSE]
    }

    # Determining how to order the cells.
    #   Simply use current order if clustering turned on.
    #   Otherwise, by cells.order, if provided, else by order.by(=Labels, default).
    if (cluster_cols) {
        order <- seq_len(nrow(scores))
    } else {
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
    }

    # Determine labels to show based on 'max.labels' with the highest
    # pre-normalized scores (not removed until later, as we still need these
    # values when normalizing the scores).
    m <- rowMaxs(scale(t(scores)))
    to.keep <- head(order(m,decreasing=TRUE), max.labels)

    # Normalize the scores between [0, 1] and cube to create more separation.
    if (normalize) {
        mmax <- rowMaxs(scores)
        mmin <- rowMins(scores)
        scores <- (scores-mmin)/(mmax-mmin)
        scores <- scores^3
    }

    scores <- scores[,seq_len(ncol(scores)) %in% to.keep,drop=FALSE]
    scores <- t(scores)

    # Create args list for making the heatmap
    args <- list(mat = scores[,order,drop=FALSE], border_color = NA,
        show_colnames = show_colnames, clustering_method = 'ward.D2',
        cluster_cols = cluster_cols, ...)

    if (normalize) {
        args$color <- viridis::viridis(100)
        args$legend_breaks <- c(0,1)
        args$legend_labels <- c("Lower", "Higher")
    } else {
        abs.max <- max(abs(range(scores)))
        breaks.len <- ifelse(is.null(color), 101, length(color)+1)
        args$breaks <- seq(-abs.max, abs.max, length.out = breaks.len)
    }

    if (ncol(annotation_col)>0) {
        args$annotation_col <- annotation_col
    }
    if (is.null(args$annotation_colors)) {
        args <- .make_heatmap_annotation_colors(args, show.pruned)
    }

    if (!is.null(color)) {
        args$color <- color
    }

    do.call(pheatmap::pheatmap, args)
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

    args$annotation_colors <- c(col_colors, row_colors)
    args
}

.make_heatmap_colors_discrete <- function(show.pruned) {
    # Creates a default vector of colors with 24*10 (overkill) options.
    annotation.colors <- rep(
        c(  # DittoSeq-v0.2.10 Colors
            "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
            "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
            "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
            "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
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
        c(  # DittoSeq-v0.2.10 Colors, distinct order
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
