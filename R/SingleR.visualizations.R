#' Plot a scatter plot of a single cell vs. a single reference cell/sample
#'
#' @param sc_data Similar to SinlgeR test input, a numeric matrix of single-cell expression values (usually
#' log-transformed or otherwise variance-stabilized), where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param cell_id a number indicating which single cell to use.
#' Alternatively, the name of the cell.
#' @param refdata Similar to SinlgeR training input, a numeric matrix of expression values. 
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' This should have the same rows as or a subset of the rows in \code{test}.
#' @param sample_id a number of the sample to use
#' @param assay.type An integer scalar or string specifying the assay of \code{sc_data} or \code{refdata} containing the relevant expression matrix,
#' if \code{sc_data} or \code{refdata} is a \linkS4class{SingleCellExperiment} object.
#'
#' @return a gpplot scatterplot of the expression of indidivual genes, with linear regression and Spearman pariwaise correlation coefficient
#' 
#' @export
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme xlab ylab ggtitle theme_classic
SingleR.DrawScatter = function(sc_data, cell_id, refdata,sample_id, assay.type = 1) {
  if (is(sc_data, "SingleCellExperiment")) {
    sc_data <- assay(sc_data, assay.type)
  }
  if (is(refdata, "SingleCellExperiment")) {
    refdata <- assay(refdata, assay.type)
  }
  if (is(refdata, "list")) {
    refdata <- refdata$data
  }
  rownames(sc_data) = tolower(rownames(sc_data))
  rownames(refdata) = tolower(rownames(refdata))
  A = intersect(rownames(sc_data),rownames(refdata))
  df = data.frame(sc_data[A,cell_id],refdata[A,sample_id])
  colnames(df) = c('x','y')
  ggplot(df,aes(x=x, y=y)) + geom_point(size=0.5,alpha=0.5,color='blue') +
    geom_smooth(method='lm',color='red')+
    theme(legend.position="none") + xlab('Single cell') + ylab('Reference sample') +
    ggtitle(paste('R =', round(1000*cor(df$x,df$y,method='spearman',use='pairwise'))/1000)) + 
    theme_classic()
}

#' Plot a heatmap of the scores for all the single cells
#'
#' @param SingleR.results the output from a run of the SingleR/classifySingleR
#' @param cells.use single cells to present, if NULL all single cells presented
#' @param types.use cell types to present, if NULL all cell types presented
#' @param clusters a clustering to present as annotation in the heatmap
#' @param top.n number of cell types to presents. Default is 40. This can have an effect on the clustering which is performed only on the cell types presented.
#' @param normalize if TRUE scores are normalized to a 0-1 scale.
#' @param order.by.clusters if TRUE columns are ordered by the input clusters, and are not clustered again
#' @param cells_order an input order for the column
#' @param silent if TRUE do not draw the plot
#' @importFrom pheatmap pheatmap
SingleR.DrawHeatmap = function(SingleR.results,cells.use = NULL, types.use = NULL,
                               clusters=NULL,top.n=40,normalize=TRUE,
                               order.by.clusters=FALSE,cells_order=NULL,silent=FALSE,
                               fontsize_row=9,...) {
  scores = SingleR.results$scores
  if (!is.null(cells.use)) {
    scores = scores[cells.use,]
  }
  if (!is.null(types.use)) {
    scores = scores[,types.use]
  }
  
  m = apply(t(scale(t(scores))),2,max)
  
  thres = sort(m,decreasing=TRUE)[min(top.n,length(m))]
  
  data = as.matrix(scores)
  
  if (normalize==TRUE) {
    mmax = rowMaxs(data)
    mmin = rowMins(data)
    data = (data-mmin)/(mmax-mmin)
    data = data^3
  }
  data = data[,m>(thres-1e-6)]
  
  data = t(data)
  
  if (!is.null(clusters)) {
    clusters = as.data.frame(clusters)
    colnames(clusters) = 'Clusters'
    rownames(clusters) = colnames(data)
    
  }
  additional_params = list(...)
  if (is.null(additional_params$annotation_colors)) {
    annotation_colors = NA
  } else {
    annotation_colors = additional_params$annotation_colors
  }
  clustering_method = 'ward.D2'
  if (order.by.clusters==TRUE) {
    data = data[,order(clusters$Clusters)]
    clusters = clusters[order(clusters$Clusters),,drop=FALSE]
    pheatmap::pheatmap(data,border_color=NA,show_colnames=FALSE,
                       clustering_method=clustering_method,fontsize_row=fontsize_row,
                       annotation_col = clusters,cluster_cols = FALSE,silent=silent, 
                       annotation_colors=annotation_colors)
  } else if (!is.null(cells_order)) {
    data = data[,cells_order]
    clusters = clusters[cells_order,,drop=FALSE]
    pheatmap::pheatmap(data,border_color=NA,show_colnames=FALSE,
                       clustering_method=clustering_method,fontsize_row=fontsize_row,
                       annotation_col = clusters,cluster_cols = FALSE,silent=silent, 
                       annotation_colors=annotation_colors)
  } else {
    if (!is.null(clusters)) {
      pheatmap::pheatmap(data,border_color=NA,show_colnames=FALSE,
                         clustering_method=clustering_method,fontsize_row=fontsize_row,
                         annotation_col = clusters,silent=silent, 
                         annotation_colors=annotation_colors)
    } else {
      pheatmap::pheatmap(data[,sample(ncol(data))],border_color=NA,show_colnames=FALSE,
                         clustering_method=clustering_method,fontsize_row=fontsize_row,
                         silent=silent, annotation_colors=annotation_colors)
      
    }
  }
}
