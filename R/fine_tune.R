# COMMENTS FROM AARON:
# We decline to use nearest neighbors here, because there's no opportunity to build a single index.
# Moreover, the data are so reduced in both dimensions that the algorithmic complexity is likely
# offset by the reduction in overhead when just computing the correlations directly.
# One could possibly improve vectorization by grouping together test cells with the same
# combination of topLabels, but it adds a lot of complexity and additional overhead.

#' @importFrom BiocParallel bplapply bpmapply SerialParam
.fine_tune_de <- function(exprs, scores, references, quantile, tune.thresh, de.info, BPPARAM=SerialParam()) {
    # Converting character vectors into integer indices.
    # We assume that SingleR() has already set up the backend.
    de.info <- bplapply(de.info, function(markers, genes, labels) {
        for (j in seq_along(markers)) {
            markers[[j]] <- match(markers[[j]], genes) - 1L
        }
        markers[labels]
    }, genes=rownames(exprs), labels=colnames(scores), BPPARAM=BPPARAM)
    de.info <- de.info[colnames(scores)]

    M <- .prep_for_parallel(exprs, BPPARAM)
    S <- .prep_for_parallel(t(scores), BPPARAM)
    bp.out <- bpmapply(Exprs=M, scores=S, FUN=fine_tune_label_de, 
        MoreArgs=list(References=references, quantile=quantile, tune_thresh=tune.thresh, marker_genes=de.info), 
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    unlist(bp.out)
}

#' @importFrom BiocParallel bpmapply SerialParam
.fine_tune_sd <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh, BPPARAM=SerialParam()) {
    M <- .prep_for_parallel(exprs, BPPARAM)
    S <- .prep_for_parallel(t(scores), BPPARAM)
    bp.out <- bpmapply(Exprs=M, scores=S, FUN=fine_tune_label_sd, 
        MoreArgs=list(References=references, quantile=quantile, tune_thresh=tune.thresh, median_mat=t(median.mat), sd_thresh=sd.thresh),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    unlist(bp.out)
}

#' @importFrom BiocParallel bpnworkers
.prep_for_parallel <- function(mat, BPPARAM) {
    is.int <- !is.double(mat[0,])
    n_cores <- bpnworkers(BPPARAM)

    if (n_cores==1L) {
        # Can't be bothered to template it twice at the C++ level,
        # as we'd have to have both int/numeric versions for the test and reference.
        if (is.int) {
            mat <- mat + 0
        }
        return(list(mat))
    }

    # Split the matrix *before* parallelization,
    # otherwise the full matrix gets serialized to all workers.
    boundaries <- as.integer(seq(from = 1L, to = ncol(mat)+1L, length.out = n_cores + 1L)) 
    out <- vector("list", n_cores)

    for (i in seq_along(out)) {
        cur_start <- boundaries[i]
        cur_end <- boundaries[i+1]
        curmat <- mat[,(cur_start - 1L) + seq_len(cur_end - cur_start),drop=FALSE]
        if (is.int) {
            curmat <- curmat + 0
        }
        out[[i]] <- curmat
    }

    out
}
