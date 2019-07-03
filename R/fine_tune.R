# COMMENTS FROM AARON:
# We decline to use nearest neighbors here, because there's no opportunity to build a single index.
# Moreover, the data are so reduced in both dimensions that the algorithmic complexity is likely
# offset by the reduction in overhead when just computing the correlations directly.
# One could possibly improve vectorization by grouping together test cells with the same
# combination of topLabels, but it adds a lot of complexity and additional overhead.

#' @importFrom BiocParallel bplapply SerialParam
.fine_tune_de <- function(exprs, scores, references, quantile, tune.thresh, de.info, BPPARAM=SerialParam()) {
    de.info <- do.call(cbind, de.info)
    if (is.null(colnames(de.info)) || !identical(colnames(de.info), rownames(de.info))) {
        stop("marker list should be named during training")
    }

    out <- bplapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell, exprs=exprs, scores=scores, 
        references=references, quantile=quantile, tune.thresh=tune.thresh, 
        commonFUN=.fine_tune_de_genes, de.info=de.info, BPPARAM=BPPARAM)            
    unlist(out)
}

#' @importFrom BiocParallel bplapply SerialParam
.fine_tune_sd <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh, BPPARAM=SerialParam()) {
    out <- bplapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell, exprs=exprs, scores=scores, 
        references=references, quantile=quantile, tune.thresh=tune.thresh, 
        commonFUN=.fine_tune_sd_genes, median.mat, sd.thresh=sd.thresh, BPPARAM=BPPARAM)
    unlist(out)
}

.fine_tune_de_genes <- function(top.labels, de.info) {
    # Finding the subset of genes (assuming 'extras' is a matrix of lists).
    all.combos <- expand.grid(top.labels, top.labels)
    Reduce(union, de.info[as.matrix(all.combos)])
}

#' @importFrom DelayedMatrixStats rowSds
.fine_tune_sd_genes <- function(top.labels, extras, sd.thresh, sd.n=500) {
    sds <- rowSds(extras, col=match(top.labels, colnames(extras)))
    sd.n <- min(length(sds), sd.n)
    sd.thresh <- min(sd.thresh, -sort(-sds, partial=sd.n)[sd.n])
    rownames(extras)[sds > sd.thresh]
}

.fine_tune_cell <- function(i, exprs, scores, references, quantile, tune.thresh, commonFUN, ...) {
    cur.exprs <- exprs[,i]
    cur.scores <- scores[i,]
    top.labels <- names(cur.scores)[cur.scores > max(cur.scores) - tune.thresh]
    old.labels <- character(0)

    # Need to compare to old.labels, to avoid an infinite loop 
    # if the correlations are still close after fine tuning.
    while (length(top.labels) > 1L && !identical(top.labels, old.labels)) {
        common <- commonFUN(top.labels, ...)
        cur.scores <- .compute_label_scores_manual(common, top.labels, cur.exprs, references, quantile=quantile)
        old.labels <- top.labels
        cur.scores <- cur.scores[!is.na(cur.scores)]
        top.labels <- names(cur.scores)[cur.scores > max(cur.scores) - tune.thresh] 
    }

    if (length(top.labels)==1L) {
        top.labels
    } else if (length(top.labels)==0L) {
        NA_character_
    } else {
        names(cur.scores)[which.max(cur.scores)]
    }
}

#' @importFrom stats cor
.compute_label_scores_manual <- function(common, top.labels, cur.exprs, references, quantile) {
    cur.exprs <- cur.exprs[common]
    cur.scores <- numeric(length(top.labels))
    names(cur.scores) <- top.labels

    for (u in top.labels) {
        ref <- references[[u]]
        ref <- as.matrix(ref[common,,drop=FALSE]) # should be cheap with few 'common'.

        # We a 'k'-based method for selecting the quantile, for consistency with classifySingleR.
        k <- max(1, round((1-quantile) * ncol(ref)))
        cur.cor <- cor(cur.exprs, ref, method="spearman")
        cur.scores[u] <- -sort(-cur.cor, partial=k)[k]
    }
    cur.scores
}
