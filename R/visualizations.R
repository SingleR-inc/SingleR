#' Plot a cell versus a reference 
#'
#' Plot a single cell's expression profile against that of a reference sample across all genes. 
#' 
#' @param test A numeric matrix of single-cell expression values where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param test.id Integer scalar or string specifying the index or name of the target cell to use.
#' @param ref A numeric matrix of reference expression values (usually log-transformed, see \code{\link{trainSingleR}}), where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param ref.id Integer scalar or string specifying the index or name of the reference cell/sample to use.
#' @param assay.type.test Integer scalar or string specifying the assay of \code{test} containing the relevant expression data.  
#' Used if \code{test} is a \linkS4class{SummarizedExperiment}.
#' @param assay.type.ref Integer scalar or string specifying the assay of \code{ref} containing the relevant expression data.  
#' Used if \code{ref} is a \linkS4class{SummarizedExperiment}.
#'
#' @return A \link[ggplot2]{ggplot} object containing a scatter plot of the cell against a reference.
#'
#' @details 
#' This function allows a user to manually check how an individual cell compares to reference cells of the same or different type.
#' It generates a scatter plot of one cell's expression profile against a chosen reference where each point is an individual gene.
#' It also displays a linear regression fit and the value of Spearman's rank correlation coefficient.
#' 
#' @author Daniel Bunis, based on code by Dvir Aran.
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#' test <- scater::logNormCounts(test)
#' 
#' # Checking if cell#1 resembles reference-set cells of its assigned label.
#' pred$labels[1]
#' same.type <- grep(pred$labels[1], sce$label)
#' 
#' # Compare expression of target cell to a reference cell of the same type
#' plotCellVsReference(test, test.id = 1, ref = sce, ref.id = same.type[1])
#' 
#' # Compare expression of target cell to reference cells of a different type
#' diff.type <- seq_along(pred$labels)[-same.type]
#' plotCellVsReference(test, test.id = 1, ref = sce, ref.id = diff.type[1])
#' 
#' @export
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
#' @importFrom stats cor
plotCellVsReference <- function(test, test.id, ref, ref.id, assay.type.test = 'logcounts', assay.type.ref = 'logcounts') {
    if (is(test, "SummarizedExperiment")) {
        test <- assay(test, assay.type.test)
    }
    if (is(ref, "SummarizedExperiment")) {
        ref <- assay(ref, assay.type.ref)
    }

    rownames(test) <- tolower(rownames(test))
    rownames(ref) <- tolower(rownames(ref))
    A <- intersect(rownames(test),rownames(ref))

    df <- data.frame(x = test[A,test.id], y = ref[A,ref.id])
    ggplot2::ggplot(df, ggplot2::aes_string(x="x", y="y")) + 
        ggplot2::geom_point(size=0.5,alpha=0.5,color='blue') +
        ggplot2::geom_smooth(method='lm',color='red') +
        ggplot2::theme(legend.position="none") + 
        ggplot2::xlab('Single cell') + 
        ggplot2::ylab('Reference sample') +
        ggplot2::ggtitle(paste('R =', format(round(cor(df$x,df$y,method='spearman',use='pairwise'), 3), 3))) + 
        ggplot2::theme_classic()
}

#' Plot a score heatmap
#'
#' Create a heatmap of the \code{\link{SingleR}} assignment scores across all cell-label combinations.
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param cells.use Integer or string vector specifying the single cells to show.
#' If \code{NULL}, all cells are presented.
#' @param labels.use String vector indicating what labels to show.
#' If \code{NULL}, all labels available in \code{results} are presented.
#' @param clusters String vector or factor containing cell cluster assignments, to be shown as annotation in the heatmap.
#' @param show.pruned Logical indicated whether prune.calls should be added as a column annotation
#' @param prune.calls Logical vector, the output of \code{\link{pruneScores}}.  This input will be unnecessary once \code{\link{pruneScores}}'s output is added to the SingleR.results dataframe
#' @param max.labels Integer scalar specifying the maximum number of labels to show.
#' @param normalize Logical specifying whether correlations should be normalized to lie in [0, 1].
#' @param order.by.clusters Logical scalar specifying if cells should be ordered by \code{clusters} and not by scores.
#' If set, this takes precedence over \code{cells.order} input.
#' @param cells.order Integer vector specifying the ordering of cells/columns of the heatmap. 
#' If set, turns off clustering of columns based on scoring.
#' Note: When used alongside \code{cells.use}, both arguments should be the same length. 
#' @param annotation_col Data.frame containing data for additional/alternative column annotations (clutering and prune-calls annotations are added automatically internally)
#' Format = row.names should be the names of the cells, and columns should be named with the title to be displayed for each annotation bar.
#' @param ... Additional parameters for heatmap control passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return A heatmap of assignment scores is generated on the current graphics device using \pkg{pheatmap}.
#'
#' @details
#' This function creates a heatmap containing the \code{\link{SingleR}} assignment scores for each cell (columns) to each reference label (rows).
#' Users can then easily identify the high-scoring labels associated with each cell and/or cluster of cells.
#'
#' If \code{max.labels} is less than the total number of unique labels, only the labels with the largest maximum scores in \code{results} are shown in the plot.
#' Specifically, the set of scores for each cell is centred and scaled, and the maximum transformed score for each label is used to choose the labels to retain.
#'
#' If \code{normalize=TRUE}, scores will be linearly adjusted for each cell so that the smallest score is 0 and the largest score is 1.
#' This is followed by cubing of the adjusted scores to improve dynamic range near 1.
#' Note that this transformation is done \emph{after} the choice of the top \code{max.labels} labels.
#'
#' @author Daniel Bunis, based on code by Dvir Aran.
#'
#' @examples
#' # Running the SingleR() example.
#' example(SingleR, echo=FALSE)
#'
#' # Creating a heatmap that shows cells showed.
#' plotScoreHeatmap(pred)
#'
#' # Creating a heatmap with clusters.
#' plotScoreHeatmap(pred, clusters=test$label)
#'
#' # We can also turn off the normalization with Normalize = FALSE
#' plotScoreHeatmap(pred, clusters=test$label, normalize = FALSE)
#' 
#' # To only show certain labels, you can use labels.use or max.labels
#' plotScoreHeatmap(pred, clusters=test$label, labels.use = c("A","B","D"))
#' plotScoreHeatmap(pred, clusters=test$label, max.labels = 4)
#' 
#' # We can pass extra tweaks the heatmap as well
#' plotScoreHeatmap(pred, clusters=test$label, fontsize.row = 9)
#' plotScoreHeatmap(pred, clusters=test$label, cutree_col = 3)
#' 
#' @export
#' @importFrom utils head
#' @importFrom DelayedArray rowMaxs rowMins
plotScoreHeatmap <- function(results, cells.use = NULL, labels.use = NULL,
    clusters = NULL, show.pruned = FALSE, prune.calls,
    max.labels = 40, normalize = TRUE,
    cells.order=NULL, order.by.clusters=FALSE, 
    annotation_col = NULL,
    ...) {
    # Add rownames (cell names) to the results DataFrame if not already there.
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    # Initiate an annotation dataframe, unless provided by the user.
    if (is.null(annotation_col)) {
        annotation_col <- data.frame(row.names = rownames(results))
    }
    # Retrieve prune.calls, add names, and add to annotation dataframe.
    if (show.pruned) {
        # prune.calls <- results$prune.calls    # UNCOMMENT after prune.calls added to results.
        names(prune.calls) <- rownames(results)
        annotation_col$Pruned <- as.character(prune.calls[rownames(annotation_col)])
    }
    # Add names to clusters and add to annotation dataframe.
    if (!is.null(clusters)) {
        names(clusters) <- rownames(results)
        annotation_col$Clusters <- clusters[rownames(annotation_col)]
    }

    # Retrieve scores and add names
    scores <- results$scores
    rownames(scores) <- rownames(results)
    
    # Trim the scores by requested cells or labels
    if (!is.null(cells.use)) {
        scores <- scores[cells.use,]
    }
    if (!is.null(labels.use)) {
        scores <- scores[,labels.use]
    }
    
    # Determining how to order the cells and create `order` vector of indices.
    cluster_cols <- FALSE
    if (order.by.clusters) {
        # Make `order` based on cluster identities (in requested cells)
        if (!is.null(cells.use)) {
            order <- order(clusters[cells.use])
        } else {
            order <- order(clusters)
        }
    } else if (!is.null(cells.order)) {
        # Make `order` based on requested order
        order <- cells.order
    } else {
        # If here, then no ordering was requested.
        #   Make `order` contain all indices, then set clustering to happen.
        order <- seq_len(ncol(scores))
        cluster_cols <- TRUE
    }
    
    # Determine labels to show based on max.labels (not removed until later)
    m <- rowMaxs(scale(t(scores)))
    to.keep <- head(order(m,decreasing=TRUE), max.labels)
    # Normalize the scores between [0, 1] and cube to create more separation.
    if (normalize) {
        mmax <- rowMaxs(scores)
        mmin <- rowMins(scores)
        scores <- (scores-mmin)/(mmax-mmin)
        scores <- scores^3
    }
    # Remove extra labels, then transpose the scores matrix.
    scores <- scores[,seq_len(ncol(scores)) %in% to.keep,drop=FALSE]
    scores <- t(scores)
    
    # Create args list for making the heatmap
    args <- list(mat = scores[,order,drop=FALSE], border_color = NA,
        show_colnames = FALSE, clustering_method = 'ward.D2',
        cluster_cols = cluster_cols, ...)
    # Add annotations dataframe to args list if it contains annotations
    if (ncol(annotation_col)>0) {
        args$annotation_col <- annotation_col
    }
    # Make the heatmap
    do.call(pheatmap::pheatmap, args)
}
