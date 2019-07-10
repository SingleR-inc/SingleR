#' Annotate scRNA-seq data
#'
#' Returns the best annotation for each cell in a training dataset,
#' given a labelled.reference dataset in the same features space.
#'
#' @param test A numeric matrix of single-cell expression values (usually log-transformed or otherwise variance-stabilized),
#' where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param training A numeric matrix of single-cell expression values, like \code{test}.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' This should have the same rows as or a subset of the rows in \code{test}.
#' @param labels A character vector or factor of known labels for all cells in \code{training}.
#' @param method String specifying whether annotation should be performed on single cells in \code{test},
#' or whether they should be aggregated into cluster-level profiles prior to annotation.
#' @param clusters A character vector or factor of cluster identities for each cell in \code{test}.
#' Only used if \code{method="cluster"}.
#' @param genes,sd.thresh Arguments controlling the genes that are used for annotation, see \code{\link{trainSingleR}}.
#' @param quantile,fine.tune,tune.thresh Further arguments to pass to \code{\link{classifySingleR}}.
#' @param assay.type An integer scalar or string specifying the assay of \code{x} containing the relevant expression matrix,
#' if \code{x} is a \linkS4class{SingleCellExperiment} object.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for building nearest neighbor indices.
#'
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This is identical to the output of \code{\link{classifySingleR}}.
#'
#' @details
#' If \code{method="single"}, this function is effectively just a convenient wrapper around \code{\link{trainSingleR}} and \code{\link{classifySingleR}}.
#' 
#' If \code{method="cluster"}, per-cell profiles are summed to obtain per-cluster profiles and annotation is performed on these clusters.
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
#' sce <- SingleCellExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' rownames(sce) <- sprintf("GENE_%s", seq_len(nrow(sce)))
#'
#' ##################################################
#' ## Mocking up some test data for classification ##
#' ##################################################
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' test <- SingleCellExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' rownames(test) <- sprintf("GENE_%s", seq_len(nrow(test)))
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
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @importFrom DelayedArray colsum
SingleR <- function(test, training, labels, method = c("single", "cluster"),
    clusters = NULL, genes = "de", quantile = 0.8, fine.tune = TRUE, 
    tune.thresh = 0.05, sd.thresh = 1, assay.type = 1, BNPARAM=KmknnParam()) 
{
    trained <- trainSingleR(training, labels, genes = genes, sd.thresh=sd.thresh, 
        assay.type=assay.type, BNPARAM=BNPARAM)

    if (is(test, "SingleCellExperiment")) {
        test <- assay(test, assay.type)
    }

    method <- match.arg(method)
    if (method=="cluster") {
        if (is.null(clusters)) {
            stop("'clusters' must be specified when 'method=\"cluster\"'")
        }
        test <- colsum(test, clusters)
    }

    # Do not set sd.thresh, use the value from 'trainSingleR'.
    classifySingleR(test, trained, quantile=quantile, fine.tune=fine.tune,
        tune.thresh=tune.thresh, assay.type=assay.type)
}
