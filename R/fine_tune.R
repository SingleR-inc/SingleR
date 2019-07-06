# COMMENTS FROM AARON:
# We decline to use nearest neighbors here, because there's no opportunity to build a single index.
# Moreover, the data are so reduced in both dimensions that the algorithmic complexity is likely
# offset by the reduction in overhead when just computing the correlations directly.
# One could possibly improve vectorization by grouping together test cells with the same
# combination of topLabels, but it adds a lot of complexity and additional overhead.

#' @importFrom BiocParallel bplapply SerialParam
.fine_tune_de <- function(exprs, scores, references, quantile, tune.thresh, de.info, BPPARAM=SerialParam()) {
    genes <- rownames(exprs)
    for (i in seq_along(de.info)) {
        for (j in seq_along(de.info[[i]])) {
            de.info[[i]][[j]] <- match(de.info[[i]][[j]], genes) - 1L
        }
    }

    # Coercing it to double as a hack; can't be bothered to template the C++ code twice.
    if (!is.double(exprs[0,])) {
        exprs <- exprs + 0
    }
    fine_tune_label_de(exprs, t(scores), references, quantile, tune.thresh, de.info) 
}

#' @importFrom BiocParallel bplapply SerialParam
.fine_tune_sd <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh, BPPARAM=SerialParam()) {
    out <- bplapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell, exprs=exprs, scores=scores, 
        references=references, quantile=quantile, tune.thresh=tune.thresh, 
        commonFUN=.fine_tune_sd_genes, median.mat, sd.thresh=sd.thresh, BPPARAM=BPPARAM)
    unlist(out)
}


