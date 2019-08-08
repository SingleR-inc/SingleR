#' Plot expression of a cell versus a reference cell/sample
#'
#' @param sc.data Numeric matrix of single-cell expression values (usually log-transformed
#' or otherwise variance-stabilized), where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param sc.id Integer scalar specifying the index of the target cell to use.
#' Alternatively, a string containing the name of the cell.
#' @param train.data Numeric matrix of reference dataset expression values (usually log-transformed
#' or otherwise variance-stabilized), where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param train.id Integer scalar specifying the reference cell/sample to use.
#' @param assay.type.sc Integer scalar or string specifying the assay of \code{sc.data} containing the relevant expression data.  
#' Used if \code{sc.data} is a \linkS4class{SingleCellExperiment}.
#' @param assay.type.train Integer scalar or string specifying the assay of \code{train.data} containing the relevant expression data.  
#' Used if provided \code{train.data} is a \linkS4class{SingleCellExperiment}.
#'
#' @return A \link{ggplot} object containing a scatter plot of the cell against a reference.
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
#' plotCellVsReference(test, sc.id = 1, train.data = sce, train.id = same.type[1])
#' 
#' # Compare expression of target cell to reference cells of a different type
#' diff.type <- seq_along(pred$labels)[-ref.sameType.as1]
#' plotCellVsReference(test, sc.id = 1, train.data = sce, train.id = diff.type[1])
#' 
#' @export
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
plotCellVsReference <- function(sc.data, sc.id, train.data,train.id, assay.type.sc = 'logcounts', assay.type.train = 'logcounts') {
    if (is(sc.data, "SingleCellExperiment")) {
        sc.data <- assay(sc.data, assay.type.sc)
    }
    if (is(train.data, "SingleCellExperiment")) {
        train.data <- assay(train.data, assay.type.train)
    }

    rownames(sc.data) <- tolower(rownames(sc.data))
    rownames(train.data) <- tolower(rownames(train.data))
    A <- intersect(rownames(sc.data),rownames(train.data))

    df <- data.frame(x = sc.data[A,sc.id], y = train.data[A,train.id])
    ggplot2::ggplot(df, ggplot2::aes_string(x="x", y="y")) + 
        ggplot2::geom_point(size=0.5,alpha=0.5,color='blue') +
        ggplot2::geom_smooth(method='lm',color='red') +
        ggplot2::theme(legend.position="none") + 
        ggplot2::xlab('Single cell') + 
        ggplot2::ylab('Reference sample') +
        ggplot2::ggtitle(paste('R =', format(round(cor(df$x,df$y,method='spearman',use='pairwise'), 3), 3))) + 
        ggplot2::theme_classic()
}

#' Plot a heatmap of the scoring of cells
#'
#' @param results A \linkS4class{DataFrame} containing the output from \code{\link{SingleR}} or \code{\link{classifySingleR}}.
#' @param cells.use Integer or character vector specifying the single cells to show.
#' If \code{NULL}, all cells are presented.
#' @param labels.use Character vector indicating what labels to show.
#' If \code{NULL}, all cell types are presented.
#' @param clusters Character vector or factor containing cell cluster assignments, to be shown as annotation in the heatmap.
#' @param max.labels Integer scalar specifying the maximum number of labels to show.
#' @param normalize Logical specifying whether correlations should be normalized to lie in [0, 1].
#' @param order.by.clusters Logical scalar specifying if cells should be ordered by \code{clusters} and not by scores.
#' If set, this takes precedence over \code{cells.order} input.
#' @param cells.order Integer vector specifying the ordering of cells/columns of the heatmap. 
#' If set, turns off clustering of columns based on scoring.
#' @param silent Logical scalar that specifying whether the plot drawn.
#' @param ... Additional parameters for heatmap control passed to \code{pheatmap()}
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
#' @export
#' @importFrom utils head
#' @importFrom DelayedArray rowMaxs rowMins
plotScoreHeatmap <- function(results, cells.use = NULL, labels.use = NULL,
    clusters=NULL, max.labels=40, normalize=TRUE,
    cells.order=NULL, order.by.clusters=FALSE, 
    fontsize.row=9, ...)
{
    scores <- results$scores
    rownames(scores) <- rownames(results)
    if(!is.null(cells.use)){
        scores <- scores[cells.use,]
    }
    if (!is.null(labels.use)) {
        scores <- scores[,labels.use]
    }

    # NOTE: calculation of 'm' before normalization is deliberate.
    m <- apply(t(scale(t(scores))),2,max)
    to.keep <- head(order(m,decreasing=TRUE), max.labels)

    if (normalize) {
        mmax <- rowMaxs(scores)
        mmin <- rowMins(scores)
        scores <- (scores-mmin)/(mmax-mmin)
        scores <- scores^3
    }

    scores <- scores[,seq_len(ncol(scores)) %in% to.keep,drop=FALSE]
    scores <- t(scores)

    # Determining how to order the cells. 
    if (!is.null(clusters)) {
        names(clusters) <- rownames(results)
        clusters <- scores.frame(Clusters = clusters[colnames(scores)], row.names = colnames(scores))
    }
    cluster_cols <- FALSE
    if (order.by.clusters && !is.null(clusters)) {
        order <- order(clusters$Clusters)
    } else if (!is.null(cells.order)){
        order <- cells.order
    } else {
        order <- seq_len(ncol(scores))
        cluster_cols <- TRUE
    }

    args <- list(mat = scores[,order,drop=FALSE], border_color = NA, show_colnames = FALSE,
        clustering_method = 'ward.D2', fontsize_row = fontsize.row,
        cluster_cols = cluster_cols, ...)
    if (!is.null(clusters)) {
        args$annotation_col <- clusters[order,,drop=FALSE]
    }
    do.call(pheatmap::pheatmap, args)
}
