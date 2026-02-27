#' Annotate scRNA-seq data
#'
#' Returns the best annotation for each cell in a test dataset,
#' given a labelled reference dataset in the same feature space.
#'
#' @param test A numeric matrix of single-cell expression values where rows are genes and columns are cells.
#' Alternatively, a \link[SummarizedExperiment]{SummarizedExperiment} object containing such a matrix.
#' @inheritParams trainSingleR
#' @param ref A numeric matrix of (usually normalized and log-transformed) expression values from a reference dataset,
#' or a \link[SummarizedExperiment]{SummarizedExperiment} object containing such a matrix;
#' see \code{\link{trainSingleR}} for details.
#'
#' Alternatively, a list or \link[S4Vectors]{List} of SummarizedExperiment objects or numeric matrices containing multiple references.
#' Row names may be different across entries but only the intersection will be used, see Details.
#' @param method Deprecated.
#' @param clusters A character vector or factor of cluster identities for each cell in \code{test}.
#' If set, annotation is performed on the aggregated cluster profiles, otherwise it defaults to per-cell annotation.
#' @param genes,sd.thresh,de.method,de.n,de.args Arguments controlling the choice of marker genes used for annotation, see \code{\link{trainSingleR}}.
#' @param aggr.ref,aggr.args Arguments controlling the aggregation of the references prior to annotation, see \code{\link{trainSingleR}}.
#' @param quantile,fine.tune,tune.thresh,fine.tune.combined,prune Further arguments to pass to \code{\link{classifySingleR}}.
#' @param assay.type.test An integer scalar or string specifying the assay of \code{test} containing the relevant expression matrix,
#' if \code{test} is a \link[SummarizedExperiment]{SummarizedExperiment} object.
#' @param assay.type.ref An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \link[SummarizedExperiment]{SummarizedExperiment} object (or is a list that contains one or more such objects).
#' @param check.missing.test Logical scalar indicating whether rows of \code{test} should be checked for missing values (and if found, removed).
#' @param check.missing.ref Logical scalar indicating whether rows of \code{ref} should be checked for missing values (and if found, removed).
#' @param check.missing Deprecated, use \code{check.missing.test} and \code{check.missing.ref} instead.
#' @param num.threads Integer scalar specifying the number of threads to use for index building and classification.
#' @param BNPARAM Deprecated and ignored.
#' @param BPPARAM Deprecated, use \code{num.threads} instead.
#'
#' @return A \link[S4Vectors]{DataFrame} is returned containing the annotation statistics for each cell (one cell per row).
#' This is identical to the output of \code{\link{classifySingleR}}.
#'
#' @details
#' This function is just a convenient wrapper around \code{\link{trainSingleR}} and \code{\link{classifySingleR}}.
#' The function will automatically restrict the analysis to the intersection of the genes in both \code{ref} and \code{test}.
#' If this intersection is empty (e.g., because the two datasets use different gene annotations), an error will be raised.
#'
#' If \code{clusters} is specified, per-cell profiles are summed to obtain per-cluster profiles.
#' Annotation is then performed by running \code{\link{classifySingleR}} on these profiles.
#' This yields a DataFrame with one row per level of \code{clusters}.
#'
#' The default settings of this function are based on the assumption that \code{ref} contains or bulk data.
#' If it contains single-cell data, this usually requires a different \code{de.method} choice.
#' Read the Note in \code{?\link{trainSingleR}} for more details.
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
#' ref <- scrapper::normalizeRnaCounts.se(ref)
#' test <- scrapper::normalizeRnaCounts.se(test)
#'
#' # Running the classification with different options:
#' pred <- SingleR(test, ref, labels=ref$label)
#' table(predicted=pred$labels, truth=test$label)
#'
#' k.out<- kmeans(t(assay(test, "logcounts")), center=5) # mock up a clustering
#' pred2 <- SingleR(test, ref, labels=ref$label, clusters=k.out$cluster) 
#' table(predicted=pred2$labels, cluster=rownames(pred2))
#'
#' @export
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
#' @importFrom DelayedArray DelayedArray
SingleR <- function(
    test, 
    ref, 
    labels, 
    method = NULL, 
    clusters = NULL, 
    genes = "de", 
    sd.thresh=1, 
    de.method ="classic", 
    de.n = NULL, 
    de.args = list(),
    aggr.ref = FALSE, 
    aggr.args = list(), 
    recompute=TRUE, 
    restrict=NULL,
    hint.sce = TRUE,
    quantile = 0.8, 
    fine.tune = TRUE, 
    tune.thresh = 0.05,
    fine.tune.combined=fine.tune,
    prune=TRUE, 
    assay.type.test = "logcounts", 
    assay.type.ref="logcounts", 
    check.missing.test=FALSE, 
    check.missing.ref=check.missing, 
    check.missing=TRUE, 
    num.threads = 1,
    BNPARAM = NULL,
    BPPARAM = NULL
) {
    num.threads <- .get_num_threads(num.threads, BPPARAM)

    # We have to do all this row-subsetting at the start before trainSingleR,
    # otherwise 'test.genes' won't match up to the filtered 'test'.
    test <- .to_clean_matrix(test, assay.type.test, check.missing.test, msg="test", num.threads=num.threads)
    tmp.ref <- ref
    if (!is.list(tmp.ref) || is.data.frame(tmp.ref)) {
        tmp.ref <- list(ref)
    }
    for (rr in tmp.ref) {
        keep <- rownames(test) %in% rownames(rr)
        if (!all(keep)) {
            test <- DelayedArray(test)[keep,,drop=FALSE] # only keeping the intersection, for safety's sake - see ?combineRecomputedResults.
        }
    }
    if (nrow(test) == 0) {
        stop("no common genes between 'test' and 'ref")
    }

    trained <- trainSingleR(
        ref, 
        labels, 
        genes = genes, 
        sd.thresh = sd.thresh, 
        de.method = de.method, 
        de.n = de.n, 
        de.args = de.args,
        aggr.ref = aggr.ref, 
        aggr.args = aggr.args, 
        recompute=recompute,
        restrict = restrict, 
        test.genes=rownames(test),
        check.missing=check.missing.ref, 
        hint.sce=hint.sce,
        BNPARAM=BNPARAM, 
        num.threads = num.threads, 
        BPPARAM=BPPARAM
    )

    if (!is.null(method)) {
        .Deprecated(msg="'method=\"cluster\"' is no longer necessary when 'cluster=' is specified")
    }

    if (!is.null(clusters)) {
        agg <- scrapper::aggregateAcrossCells(test, list(clusters=clusters), num.threads=num.threads)
        test <- agg$sums
        colnames(test) <- agg$combinations$clusters
    }

    classifySingleR(
        test, 
        trained, 
        quantile=quantile, 
        fine.tune=fine.tune,
        tune.thresh=tune.thresh, 
        fine.tune.combined=fine.tune.combined,
        prune=prune, 
        check.missing=FALSE, 
        num.threads = num.threads, 
        BPPARAM=BPPARAM
    )
}
