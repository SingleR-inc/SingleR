#' Annotate scRNA-seq data
#'
#' Returns the best annotation for each cell in a test dataset,
#' given a labelled reference dataset in the same feature space.
#'
#' @param test A numeric matrix of single-cell expression values where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @inheritParams trainSingleR
#' @param ref A numeric matrix of (usually log-transformed) expression values from a reference dataset,
#' or a \linkS4class{SummarizedExperiment} object containing such a matrix;
#' see \code{\link{trainSingleR}} for details.
#'
#' Alternatively, a list or \linkS4class{List} of SummarizedExperiment objects or numeric matrices containing multiple references.
#' Row names may be different across entries but only the intersection will be used, see Details.
#' @param method String specifying whether annotation should be performed on single cells in \code{test},
#' or whether they should be aggregated into cluster-level profiles prior to annotation.
#' @param clusters A character vector or factor of cluster identities for each cell in \code{test}.
#' Only used if \code{method="cluster"}.
#' @param genes,sd.thresh,de.method,de.n,de.args Arguments controlling the choice of marker genes used for annotation, see \code{\link{trainSingleR}}.
#' @param aggr.ref,aggr.args Arguments controlling the aggregation of the references prior to annotation, see \code{\link{trainSingleR}}.
#' @param quantile,fine.tune,tune.thresh,prune Further arguments to pass to \code{\link{classifySingleR}}.
#' @param assay.type.test An integer scalar or string specifying the assay of \code{test} containing the relevant expression matrix,
#' if \code{test} is a \linkS4class{SummarizedExperiment} object.
#' @param assay.type.ref An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object (or is a list that contains one or more such objects).
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for building nearest neighbor indices.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed, if any.
#'
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This is identical to the output of \code{\link{classifySingleR}}.
#'
#' @details
#' If \code{method="single"}, this function is effectively just a convenient wrapper around \code{\link{trainSingleR}} and \code{\link{classifySingleR}}.
#' 
#' If \code{method="cluster"}, per-cell profiles are summed to obtain per-cluster profiles and annotation is performed on these clusters.
#'
#' The function will automatically restrict the analysis to the intersection of the genes available in both \code{ref} and \code{test}.
#' If this intersection is empty (e.g., because the two datasets use different annotation in their row names), an error will be raised.
#'
#' \code{ref} can contain both single-cell or bulk data, but in the case of the former, read the Note in \code{?\link{trainSingleR}}.
#' 
#' @references
#' Aran D, Looney AP, Liu L et al. (2019).
#' Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
#' \emph{Nat. Immunology} 20, 163–172.
#'
#' @author Aaron Lun, based on code by Dvir Aran.
#' @examples
#' # Mocking up data with log-normalized expression values:
#' ref <- .mockRefData()
#' test <- .mockTestData(ref)
#'
#' ref <- scuttle::logNormCounts(ref)
#' test <- scuttle::logNormCounts(test)
#'
#' # Running the classification with different options:
#' pred <- SingleR(test, ref, labels=ref$label)
#' table(predicted=pred$labels, truth=test$label)
#'
#' k.out<- kmeans(t(assay(test, "logcounts")), center=5) # mock up a clustering
#' pred2 <- SingleR(test, ref, labels=ref$label, 
#'     method="cluster", clusters=k.out$cluster) 
#' table(predicted=pred2$labels, cluster=rownames(pred2))
#'
#' @export
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
#' @importFrom DelayedArray colsum DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom BiocParallel SerialParam
SingleR <- function(test, ref, 
    labels, method = c("single", "cluster"), clusters = NULL, 
    genes = "de", sd.thresh=1, de.method ="classic", de.n = NULL, de.args = list(),
    aggr.ref = FALSE, aggr.args = list(), recompute=TRUE, restrict=NULL,
    quantile = 0.8, fine.tune = TRUE, tune.thresh = 0.05, prune=TRUE, 
    assay.type.test = "logcounts", assay.type.ref="logcounts", 
    check.missing=TRUE, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    test <- .to_clean_matrix(test, assay.type.test, check.missing, msg="test")

    # Converting to a common list format for ease of data munging.
    if (single.ref <- !.is_list(ref)) {
        ref <- list(ref)
    }

    ref <- lapply(ref, FUN=.to_clean_matrix, assay.type=assay.type.ref, 
        check.missing=check.missing, msg="ref")
    refnames <- Reduce(intersect, lapply(ref, rownames))

    keep <- intersect(rownames(test), refnames)
    if (length(keep) == 0) {
        stop("no common genes between 'test' and 'ref'")
    }
    if (!identical(keep, rownames(test))) {
        test <- test[keep,]
    }
    for (i in seq_along(ref)) {
        if (!identical(keep, rownames(ref[[i]]))) {
            ref[[i]] <- ref[[i]][keep,,drop=FALSE]
        }
    }

    # Converting back.
    if (single.ref) {
        ref <- ref[[1]]
    }

    trained <- trainSingleR(ref, labels, genes = genes, sd.thresh = sd.thresh, 
        de.method = de.method, de.n = de.n, de.args = de.args,
        aggr.ref = aggr.ref, aggr.args = aggr.args, recompute=recompute,
        restrict = restrict, check.missing=FALSE, BNPARAM=BNPARAM)

    method <- match.arg(method)
    if (method=="cluster") {
        if (is.null(clusters)) {
            stop("'clusters' must be specified when 'method=\"cluster\"'")
        }

        oldp <- getAutoBPPARAM()
        setAutoBPPARAM(BPPARAM)
        on.exit(setAutoBPPARAM(oldp), add=TRUE)
        test <- colsum(DelayedArray(test), clusters)
    }

    # Do not set sd.thresh, use the value from 'trainSingleR'.
    classifySingleR(test, trained, quantile=quantile, fine.tune=fine.tune,
        tune.thresh=tune.thresh, prune=prune, check.missing=FALSE, BPPARAM=BPPARAM)
}
