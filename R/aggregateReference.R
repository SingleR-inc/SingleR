#' Aggregate reference samples
#'
#' Aggregate reference samples for a given label by averaging their count profiles.
#' This can be done with varying degrees of resolution to preserve the within-label heterogeneity.
#'
#' @param ref A numeric matrix of reference expression values, usually containing log-expression values.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param labels A character vector or factor of known labels for all cells in \code{ref}.
#' @param power Numeric scalar between 0 and 1 indicating how much aggregation should be performed, see Details.
#' @param assay.type An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object.
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' 
#' @details
#' With single-cell reference datasets, it is often useful to aggregate individual cells into pseudo-bulk samples to serve as a reference.
#' This improves speed (and to some extent, reduces noise) in downstream assignment with \code{\link{classifySingleR}} or \code{\link{SingleR}}.
#' 
#' The most obvious aggregation is to simply average all counts for all cells in a label to obtain a single pseudo-bulk profile.
#' This can be achieved by setting \code{power=0}.
#' However, this discards information about the within-label heterogeneity (e.g., the \dQuote{shape} and spread of the population in expression space)
#' that may be informative for assignment, especially for closely related labels.
#'
#' Instead, the default approach in this function is to create a series of pseudo-bulk samples to represent each label.
#' This is achieved by performing vector quantization using k-means clustering on all cells in a particular label.
#' Cells in each cluster are subsequently averaged to create one pseudo-bulk sample.
#' We set the number of clusters to be \code{ncol(ref)^power} so that labels with more cells have more resolved representatives.
#'
#' If \code{power=1}, no aggregation is performed.
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
#' library(scater)
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
#' @importFrom Matrix rowSums
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedArray sweep colsum DelayedArray
aggregateReference <- function(ref, labels, power=0.5, assay.type="logcounts", check.missing=TRUE) {
    output.vals <- output.labs <- list()
    ref <- .to_clean_matrix(ref, assay.type, check.missing, msg="ref")

    for (u in unique(labels)) {
        chosen <- u==labels
        current <- ref[,chosen,drop=FALSE]

        if (power==0) {
            val <- matrix(rowMeans(current), dimnames=list(rownames(current), NULL))
        } else if (power==1) {
            val <- current
        } else {
            ncenters <- ncol(current)^power
            clustered <- kmeans(as.matrix(t(current)), centers=ncenters)
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
