#' Aggregate reference samples
#'
#' Aggregate reference samples for a given label by averaging their count profiles.
#' This can be done with varying degrees of resolution to preserve the within-label heterogeneity.
#'
#' @param ref A numeric matrix of reference expression values, usually containing log-expression values.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param labels A character vector or factor of known labels for all cells in \code{ref}.
#' @param ncenters Integer scalar specifying the maximum number of aggregated profiles to produce for each label.
#' @param power Numeric scalar between 0 and 1 indicating how much aggregation should be performed, see Details.
#' @param rank Integer scalar specfiying the number of principal components to use during clustering.
#' @param assay.type An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object.
#' @param subset.row Integer, character or logical vector indicating the rows of \code{ref} to use for k-means clustering. 
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how parallelization should be performed.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object indicating which SVD algorithm should be used in \code{\link{runPCA}}.
#' 
#' @details
#' With single-cell reference datasets, it is often useful to aggregate individual cells into pseudo-bulk samples to serve as a reference.
#' This improves speed in downstream assignment with \code{\link{classifySingleR}} or \code{\link{SingleR}}.
#' The most obvious aggregation is to simply average all counts for all cells in a label to obtain a single pseudo-bulk profile.
#' However, this discards information about the within-label heterogeneity (e.g., the \dQuote{shape} and spread of the population in expression space)
#' that may be informative for assignment, especially for closely related labels.
#'
#' The default approach in this function is to create a series of pseudo-bulk samples to represent each label.
#' This is achieved by performing vector quantization via k-means clustering on all cells in a particular label.
#' Cells in each cluster are subsequently averaged to create one pseudo-bulk sample that serves as a representative for that location in the expression space.
#' This reduces the number of separate observations (for speed) while preserving some level of population heterogeneity (for fidelity).
#' 
#' The number of pseudo-bulk samples per label is controlled by \code{ncenters}.
#' By default, we set the number of clusters to \code{X^power} where \code{X} is the number of cells for that label.
#' This ensures that labels with more cells have more resolved representatives.
#' If \code{ncenters} is greater than the number of samples for a label and/or \code{power=1}, no aggregation is performed.
#' Setting \code{power=0} will aggregate all cells of a label into a single pseudo-bulk profile.
#'
#' k-means clustering is actually performed on the first \code{rank} principal components as computed using \code{\link{runPCA}}.
#' The use of PCs compacts the data for more efficient operation of \code{\link{kmeans}};
#' it also removes some of the high-dimensional noise to highlight major factors of within-label heterogenity.
#' Note that the PCs are only used for clustering and the full expression profiles are still used for the final averaging.
#' Users can disable the PCA step by setting \code{rank=Inf}.
#'
#' Setting \code{subset.row} to an appropriate indexing vector will subset the matrix by gene prior to PCA.
#' However, again, the aggregation is performed on the full expression profiles.
#' This option may be useful for clustering based on known genes of interest but retaining all genes in the aggregated results.
#'
#' We use the average rather than the sum in order to be compatible with \code{\link{trainSingleR}}'s internal marker detection.
#' Moreover, unlike counts, the sum of transformed and normalized expression values generally has little meaning.
#' We do not use the median to avoid consistently obtaining zeros for lowly expressed genes.
#' 
#' @return
#' A \linkS4class{SummarizedExperiment} object with a \code{"logcounts"} assay containing a matrix of aggregated expression values,
#' and a \code{label} column metadata field specifying the label corresponding to each column.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' # Making up some labels for demonstration purposes:
#' labels <- sample(LETTERS, ncol(sce), replace=TRUE)
#'
#' # Aggregation at different resolutions:
#' (aggr <- aggregateReference(sce, labels, power=0.5))
#'
#' (aggr <- aggregateReference(sce, labels, power=0))
#'
#' # No aggregation:
#' (aggr <- aggregateReference(sce, labels, power=1))
#' @export
#' @importFrom stats kmeans
#' @importFrom Matrix rowSums t
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray sweep colsum DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular bsparam runPCA
aggregateReference <- function(ref, labels, ncenters=NULL, power=0.5, assay.type="logcounts", 
    rank=20, subset.row=NULL, check.missing=TRUE, BPPARAM=SerialParam(), BSPARAM=bsparam())
{
    output.vals <- output.labs <- list()
    ref <- .to_clean_matrix(ref, assay.type, check.missing, msg="ref")

    oldp <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldp), add=TRUE)

    for (u in unique(labels)) {
        chosen <- u==labels
        current <- ref[,chosen,drop=FALSE]

        cur.ncenters <- ncenters
        if (is.null(cur.ncenters)) {
            cur.ncenters <- floor(ncol(current)^power)
        }

        if (cur.ncenters==0) {
            val <- matrix(rowMeans(current), dimnames=list(rownames(current), NULL))
        } else if (cur.ncenters >= ncol(current)) {
            val <- current
        } else {
            to.use <- t(current)
            if (!is.null(subset.row)) {
                to.use <- to.use[,subset.row,drop=FALSE]
            }

            # Identifying the top PCs to avoid realizing the entire matrix.
            if (rank <= min(dim(to.use))-1L) {
                pcs <- runPCA(to.use, rank=rank, get.rotation=FALSE, BSPARAM=BSPARAM, BPPARAM=BPPARAM)$x
            } else {
                pcs <- as.matrix(to.use)
            }

            clustered <- kmeans(pcs, centers=cur.ncenters)
            val <- colsum(DelayedArray(current), clustered$cluster)
            tab <- table(clustered$cluster)[colnames(val)]
            val <- sweep(val, 2, tab, "/")
        }

        colnames(val) <- sprintf("%s.%s", u, seq_len(ncol(val)))
        output.vals[[u]] <- val
        output.labs[[u]] <- rep(u, ncol(val))
    }

    if (length(output.vals)==0L) {
        output.vals[[1]] <- matrix(0, nrow(ref), 0, dimnames=list(rownames(ref), NULL))
        output.labs[[1]] <- labels[0]
    }

    SummarizedExperiment(list(logcounts=do.call(cbind, output.vals)),
        colData=DataFrame(label=unlist(output.labs, use.names=FALSE)))
}
