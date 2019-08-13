#' Annotate scRNA-seq data
#'
#' Returns the best annotation for each cell in a test dataset,
#' given a labelled reference dataset in the same feature space.
#'
#' @param test A numeric matrix of single-cell expression values where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param ref A numeric matrix of reference expression values (usually log-transformed, see \code{\link{trainSingleR}}).
#' This should have the same rows as or a subset of the rows in \code{test}.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param labels A character vector or factor of known labels for all cells in \code{ref}.
#' @param method String specifying whether annotation should be performed on single cells in \code{test},
#' or whether they should be aggregated into cluster-level profiles prior to annotation.
#' @param clusters A character vector or factor of cluster identities for each cell in \code{test}.
#' Only used if \code{method="cluster"}.
#' @param genes,sd.thresh Arguments controlling the genes that are used for annotation, see \code{\link{trainSingleR}}.
#' @param quantile,fine.tune,tune.thresh Further arguments to pass to \code{\link{classifySingleR}}.
#' @param assay.type.test An integer scalar or string specifying the assay of \code{test} containing the relevant expression matrix,
#' if \code{test} is a \linkS4class{SummarizedExperiment} object.
#' @param assay.type.ref An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object.
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
#' @author Aaron Lun, based on code by Dvir Aran.
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
#' sce <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' 
#' rownames(sce) <- sprintf("GENE_%s", seq_len(nrow(sce)))
#' sce <- scater::logNormCounts(sce)
#'
#' ##################################################
#' ## Mocking up some test data for classification ##
#' ##################################################
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' test <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#'
#' rownames(test) <- sprintf("GENE_%s", seq_len(nrow(test)))
#' test <- scater::logNormCounts(test)
#' 
#' ###############################
#' ## Performing classification ##
#' ###############################
#' 
#' pred <- SingleR(test, sce, labels=sce$label)
#' table(predicted=pred$labels, truth=g)
#'
#' pred2 <- SingleR(test, sce, labels=sce$label, 
#'     method="cluster", clusters=test$label) 
#' table(predicted=pred2$labels, truth=rownames(pred2))
#'
#' @export
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
#' @importFrom DelayedArray colsum DelayedArray
#' @importFrom BiocParallel SerialParam
SingleR <- function(test, ref, labels, method = c("single", "cluster"),
    clusters = NULL, genes = "de", quantile = 0.8, fine.tune = TRUE, 
    tune.thresh = 0.05, sd.thresh = 1, assay.type.test = "logcounts", assay.type.ref="logcounts", 
    check.missing=TRUE, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    test <- .to_clean_matrix(test, assay.type.test, check.missing, msg="test")
    ref <- .to_clean_matrix(ref, assay.type.ref, check.missing, msg="ref")

    keep <- intersect(rownames(test), rownames(ref))
    if (length(keep) == 0) {
        stop("no common genes between 'test' and 'ref'")
    }
    if (!identical(keep, rownames(test))) {
        test <- test[keep,]
    } 
    if (!identical(keep, rownames(ref))) {
        ref <- ref[keep,]
    }

    trained <- trainSingleR(ref, labels, genes = genes, sd.thresh=sd.thresh, 
        check.missing=FALSE, BNPARAM=BNPARAM)

    method <- match.arg(method)
    if (method=="cluster") {
        if (is.null(clusters)) {
            stop("'clusters' must be specified when 'method=\"cluster\"'")
        }
        test <- colsum(DelayedArray(test), clusters)
    }

    # Do not set sd.thresh, use the value from 'trainSingleR'.
    classifySingleR(test, trained, quantile=quantile, fine.tune=fine.tune,
        tune.thresh=tune.thresh, check.missing=FALSE, BPPARAM=BPPARAM)
}
