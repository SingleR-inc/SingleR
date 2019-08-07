#' Plot expression of a cell versus a reference cell/sample
#'
#' @param sc.data Numeric matrix of single-cell expression values (usually log-transformed
#' or otherwise variance-stabilized), where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param sc.id Integer scalar specifying the index of the target cell to use.
#' Alternatively, the name of the cell.
#' @param train.data Numeric matrix of reference dataset expression values (usually log-transformed
#' or otherwise variance-stabilized), where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param train.id Integer scalar specifying the reference cell/sample to use.
#' @param assay.type.sc Integer scalar or Character specifying the assay of \code{sc.data} containing the relevant assay if provided as a \linkS4class{SingleCellExperiment}.
#' @param assay.type.train Integer scalar or Character specifying the assay of \code{train.data} containing the relevant assay if provided as a \linkS4class{SingleCellExperiment}.
#' if \code{sc.data} or \code{train.data} is a \linkS4class{SingleCellExperiment} object.
#' @return a gpplot scatterplot of the expression of indidivual genes, with linear regression and Spearman pariwaise correlation coefficient
#' @details This function allows a user to manually check how cells compare to reference cells of the same type or of a different type.
#' @author Daniel Bunis, based on code by Dvir Aran.
#' @examples
#' ###########################################
#' ## Mocking up some example training data ##
#' ###########################################
#'
#' Ngroups <- 5
#' Ngenes <- 1000
#' means <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
#' means[1:900,] <- 0
#' colnames(means) <- LETTERS[1:5]
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' train.data <- SingleCellExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' rownames(train.data) <- sprintf("GENE_%s", seq_len(nrow(train.data)))
#'
#' ##################################################
#' ## Mocking up some test data for classification ##
#' ##################################################
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' sc.data <- SingleCellExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' rownames(sc.data) <- sprintf("GENE_%s", seq_len(nrow(sc.data)))
#' 
#' #####################
#' ## Running SingleR ##
#' #####################
#' 
#' pred <- SingleR(sc.data, train.data, labels=train.data$label)
#' table(predicted=pred$labels, truth=g)
#' 
#' ##### Check if cell#1 resembles reference-set cells of the same type
#' # Identify reference cells of the same type as cell#1
#'   #This cell is type...
#'   pred$labels[1]
#' ref.sameType.as1 <- grep(pred$labels[1], train.data$label)
#' 
#' # Compare expression of target cell to reference cells of the same type
#' plotCellVsReference(sc.data, sc.id = 1,
#'                     train.data, train.id = ref.sameType.as1[1])
#' 
#' # Compare expression of target cell to reference cells of a different type
#' ref.diffType.as1 <- seq_along(train.data$label)[-ref.sameType.as1]
#' plotCellVsReference(sc.data, sc.id = 1,
#'                     train.data, train.id = ref.diffType.as1[1])
#' 
#' @export
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme xlab ylab ggtitle theme_classic
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
plotCellVsReference <- function(sc.data, sc.id, train.data,train.id, assay.type.sc = 1, assay.type.train = 1) {
    if (is(sc.data, "SingleCellExperiment")) {
        sc.data <- assay(sc.data, assay.type.sc)
    }
    if (is(train.data, "SingleCellExperiment")) {
        train.data <- assay(train.data, assay.type.train)
    }
    if (is(train.data, "list")) {
        train.data <- train.data$data
    }
    rownames(sc.data) <- tolower(rownames(sc.data))
    rownames(train.data) <- tolower(rownames(train.data))
    A <- intersect(rownames(sc.data),rownames(train.data))
    df <- data.frame(sc.data[A,sc.id],train.data[A,train.id])
    colnames(df) <- c('x','y')
    ggplot(df,aes(x=x, y=y)) + geom_point(size=0.5,alpha=0.5,color='blue') +
        geom_smooth(method='lm',color='red')+
        theme(legend.position="none") + xlab('Single cell') + ylab('Reference sample') +
        ggtitle(paste('R =', round(1000*cor(df$x,df$y,method='spearman',use='pairwise'))/1000)) + 
        theme_classic()
}

#' Plot a heatmap of the scoring of cells
#'
#' @param SingleR.results DataFrame output from a run of \code{SingleR()} or \code{classifySingleR()}
#' @param cells.use Integer or Character vector of single cells to present. If NULL, all cells are presented.
#' @param labels.use Character vector indicating what cell types to present. If NULL, all cell types are presented.
#' @param clustering Character vector of Factor representing cell clustering to be present as annotation in the heatmap.
#' @param max.labels Integer scalar representing the number of cell types to present. Default is 40. This can have an effect on the clustering which is performed only on the cell types presented.
#' @param normalize Logical that sets if scores are normalized to a 0-1 scale.  Default is TRUE.
#' @param order.by.clusters Logical that sets if cells are ordered by the input clusters and not clustered based on scoring.  Default is FALSE. Takes precedence over \code{cells.order} input.
#' @param cells.order Integer vector with length equal to the number of cells in the heatmap. Sets the ordering of cells/columns of the heatmap. Turns off clustering of columns based on scoring.
#' @param silent Logical that sets whether the plot drawn.
#' @param ... Additional parameters for heatmap control passed to \code{pheatmap()}
#' @export
#' @importFrom pheatmap pheatmap
plotScoreHeatmap <- function(SingleR.results, cells.use = NULL, labels.use = NULL,
                             clusters=NULL, max.labels=40, normalize=TRUE,
                             cells.order=NULL, order.by.clusters=FALSE, silent=FALSE,
                             fontsize.row=9, ...) {
    scores <- SingleR.results$scores
    if (!is.null(cells.use)) {
        scores <- scores[cells.use,]
    }
    if (!is.null(labels.use)) {
        scores <- scores[,labels.use]
    }
    m <- apply(t(scale(t(scores))),2,max)
    thres <- sort(m,decreasing=TRUE)[min(max.labels,length(m))]
    data <- as.matrix(scores)
    if (normalize==TRUE) {
        mmax <- rowMaxs(data)
        mmin <- rowMins(data)
        data <- (data-mmin)/(mmax-mmin)
        data <- data^3
    }
    data <- data[,m>(thres-1e-6)]
    data <- t(data)
    if (!is.null(clusters)) {
        clusters <- as.data.frame(clusters)
        colnames(clusters) <- 'Clusters'
        rownames(clusters) <- colnames(data)
    }
    cluster_cols <- FALSE
    if (order.by.clusters==TRUE) {
        order <- order(clusters$Clusters)
    } else if (!is.null(cells.order)){
        order <- cells.order
    } else {
        order <- seq_len(ncol(data))
        cluster_cols <- TRUE
    }
    args <- list(mat = data[,order], border_color = NA, show_colnames = FALSE,
                 clustering_method = 'ward.D2', fontsize_row = fontsize.row,
                 silent = silent, cluster_cols = cluster_cols, ...)
    if (!is.null(clusters)) {
        args$annotation_col <- clusters[order,,drop=FALSE]
    }
    do.call(pheatmap::pheatmap, args)
}
