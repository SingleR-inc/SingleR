# COMMENTS FROM AARON:
# We decline to use nearest neighbors here, because there's no opportunity to build a single index.
# Moreover, the data are so reduced in both dimensions that the algorithmic complexity is likely
# offset by the reduction in overhead when just computing the correlations directly.
# One could possibly improve vectorization by grouping together test cells with the same
# combination of topLabels, but it adds a lot of complexity and additional overhead.

#' @importFrom BiocParallel bplapply SerialParam
.fine_tune_data <- function(exprs, scores, genes, extras, references, quant.thresh, tune.thresh, sd.thresh, sd.n, BPPARAM=SerialParam()) {
    if (genes=="de") {
        extras <- do.call(cbind, extras)
        out <- bplapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell_de, exprs=exprs, 
            scores=scores, extras=extras, references=references, quant.thresh=quant.thresh, 
            tune.thresh=tune.thresh, BPPARAM=BPPARAM)            

    } else if (genes=="sd") {
        out <- bplapply(seq_len(ncol(exprs)), FUN=.fine_tune_cell_sd, exprs=exprs, 
            scores=scores, extras=extras, references=references, quant.thresh=quant.thresh, 
            tune.thresh=tune.thresh, sd.thresh=sd.thresh, sd.n=sd.n, BPPARAM=BPPARAM)
    }
    unlist(out)
}

.fine_tune_cell_de <- function(i, exprs, scores, extras, references, quant.thresh, tune.thresh) {
    .fine_tune_cell(i, exprs, scores, references, quant.thresh, tune.thresh, commonFUN=.fine_tune_de_genes, 
        extras=extras)
}

.fine_tune_de_genes <- function(top.labels, extras) {
    # Finding the subset of genes (assuming 'extras' is a matrix of lists).
    all.combos <- expand.grid(top.labels, top.labels)
    unlist(extras[as.matrix(all.combos)])
}

.fine_tune_cell_sd <- function(i, exprs, scores, extras, references, quant.thresh, tune.thresh, sd.thresh, sd.n) {
    .fine_tune_cell(i, exprs, scores, references, quant.thresh, tune.thresh, commonFUN=.fine_tune_sd_genes, 
        extras=extras, sd.thresh=sd.thresh, sd.n=sd.n)
}

#' @importFrom DelayedMatrixStats rowSds
.fine_tune_sd_genes <- function(top.labels, extras, sd.thresh, sd.n) {
    sds <- rowSds(extras[,top.labels])
    sd.n <- min(length(sds), sd.n)
    sd.thresh <- min(sd.thresh, sort(sds, partial=sd.n)[sd.n])
    rownames(extras)[sds > sd.thresh]
}

#' @importFrom stats cor quantile
.fine_tune_cell <- function(i, exprs, scores, references, quant.thresh, tune.thresh, commonFUN, ...) {
    cur.exprs <- exprs[,i]
    cur.scores <- scores[i,]
    top.labels <- names(cur.scores)[cur.scores > max(cur.scores) - tune.thresh]
    old.labels <- character(0)

    # Need to compare to old.labels, to avoid an infinite loop 
    # if the correlations are still close after fine tuning.
    while (length(top.labels) > 1L && !identical(top.labels, old.labels)) {
        common <- commonFUN(top.labels, ...)

        cur.scores <- numeric(length(top.labels))
        names(cur.scores) <- top.labels
        for (u in top.labels) {
            ref <- references[[u]]
            ref <- as.matrix(ref[common,,drop=FALSE]) # should be cheap with few 'common'.
            cur.scores[u] <- quantile(cor(cur.exprs[common], ref, method="spearman"), 
                na.rm=TRUE, p=quant.thresh)
        }

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
