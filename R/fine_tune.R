# COMMENTS FROM AARON:
# We decline to use nearest neighbors here, because there's no opportunity to build a single index.
# Moreover, the data are so reduced in both dimensions that the algorithmic complexity is likely
# offset by the reduction in overhead when just computing the correlations directly.
# One could possibly improve vectorization by grouping together test cells with the same
# combination of topLabels, but it adds a lot of complexity and additional overhead.

#' @importFrom BiocParallel SerialParam
#' @importFrom beachmat colBlockApply
.fine_tune_de <- function(exprs, scores, references, quantile, tune.thresh, de.info, BPPARAM=SerialParam()) {
    # Checking that all names are in sync.
    stopifnot(identical(names(references), colnames(scores)))
    stopifnot(identical(names(references), names(de.info)))
    for (markers in de.info) {
        stopifnot(identical(names(markers), names(de.info)))
    }

    # Scanning across all references and subsetting to the common genes.
    # This should reduce the amount of data that gets distributed,
    # as well as the number of cache misses.
    universe <- unique(unlist(lapply(de.info, unlist, use.names=FALSE), use.names=FALSE))
    references <- lapply(references, function(x) x[universe,,drop=FALSE])
    exprs <- exprs[universe,,drop=FALSE]

    # Converting character vectors into integer indices.
    de.info <- lapply(de.info, function(markers) {
        lapply(markers, function(x) match(x, universe) - 1L)
    })

    references <- .realize_references(references)

    # We assume that classifySingleR() has already set up the backend.
    bp.out <- colBlockApply(exprs, FUN=.fine_tune_de0, BPPARAM=BPPARAM,
        scores=t(scores), References=references, quantile=quantile, 
        tune_thresh=tune.thresh, marker_genes=de.info)

    do.call(mapply, c(bp.out, list(FUN=c, SIMPLIFY=FALSE, USE.NAMES=FALSE)))
}

#' @importFrom DelayedArray makeNindexFromArrayViewport
.fine_tune_de0 <- function(block, scores, ...) {
    vp <- attr(block, "from_grid")[[attr(block, "block_id")]]
    idx <- makeNindexFromArrayViewport(vp, expand.RangeNSBS = TRUE)[[2]]
    if (!is.null(idx)) {
        scores <- scores[,idx,drop=FALSE]
    }

    fine_tune_label_de(block, scores, ...)
}

#' @importFrom BiocParallel SerialParam
#' @importFrom beachmat colBlockApply
.fine_tune_sd <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh, BPPARAM=SerialParam()) {
    stopifnot(identical(names(references), colnames(scores)))

    references <- .realize_references(references)

    bp.out <- colBlockApply(exprs, FUN=.fine_tune_sd0, BPPARAM=BPPARAM, 
        scores=t(scores), References=references, quantile=quantile, 
        tune_thresh=tune.thresh, median_mat=t(median.mat), sd_thresh=sd.thresh)

    do.call(mapply, c(bp.out, list(FUN=c, SIMPLIFY=FALSE, USE.NAMES=FALSE)))
}

#' @importFrom DelayedArray makeNindexFromArrayViewport
.fine_tune_sd0 <- function(block, scores, ...) {
    vp <- attr(block, "from_grid")[[attr(block, "block_id")]]
    idx <- makeNindexFromArrayViewport(vp, expand.RangeNSBS = TRUE)[[2]]
    if (!is.null(idx)) {
        scores <- scores[,idx,drop=FALSE]
    }
    fine_tune_label_sd(block, scores, ...)
}
