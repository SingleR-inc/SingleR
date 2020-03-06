#' @importFrom BiocGenerics rbind 
#' @importFrom BiocParallel SerialParam bpiterate
#' @importFrom S4Vectors DataFrame selfmatch
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom DelayedArray rowMaxs
.fine_tune_loop <- function(exprs, scores, references, 
    quantile, tune.thresh, defineFUN, 
    BPPARAM=SerialParam(), BNPARAM=KmknnParam())
{
    indices <- list()
    present <- DataFrame(scores[0,])
    all.labels <- colnames(scores)

    survivors <- seq_len(nrow(scores))
    remaining <- exprs
    last.nhits <- integer(nrow(scores))

    final.labels <- character(nrow(scores))
    final.best <- final.second <- numeric(nrow(scores))

    repeat {
        above <- scores >= rowMaxs(scores) - tune.thresh
        nhits <- rowSums(above)
        finished <- nhits==1L | nhits==last.nhits[survivors]
        last.nhits[survivors] <- nhits

        # TODO: make this faster?
        for (i in which(finished)) {        
            curscores <- sort(scores[i,], decreasing=TRUE)
            original <- survivors[i]
            final.labels[original] <- names(curscores)[1]
            final.best[original] <- curscores[1]
            final.second[original] <- curscores[2]
        }

        if (all(finished)) {
           break
        } else if (any(finished)) {
            survivors <- survivors[!finished]
            above <- above[!finished,,drop=FALSE]
            remaining <- remaining[,!finished,drop=FALSE]
        }

        encoding <- DataFrame(above, check.names=FALSE)
        use <- selfmatch(encoding)

        ITER <- .fine_tune_iterator(use, encoding=encoding, references=references, 
            defineFUN=defineFUN, exprs=remaining)

        output <- bpiterate(ITER=ITER, FUN=.fine_tune_searcher, quantile=quantile,
            all.labels=all.labels, BPPARAM=BPPARAM, BNPARAM=BNPARAM)

        scores <- do.call(rbind, output)
        scores[order(use),] <- scores # need to reorder as iterator splits up by 'use'.
    }

    list(final.labels, final.best, final.second)
}

.fine_tune_iterator <- function(use, encoding, references, defineFUN, exprs) {
    u.use <- sort(unique(use))
    counter <- 1L

    function() {
        if (counter > length(u.use)) {
            return(NULL)
        }

        in.use <- u.use[counter]
        used <- use==in.use
        
        chosen <- unlist(encoding[which(used)[1],], use.names=FALSE)
        chosen <- colnames(encoding)[chosen]
        common <- defineFUN(chosen)

        counter <<- counter + 1L
        list(
            exprs=exprs[common,used,drop=FALSE],
            references=lapply(references[chosen], function(x) x[common,,drop=FALSE])
        )
    }
}

#' @importFrom BiocNeighbors buildIndex
.fine_tune_searcher <- function(data, quantile, all.labels, BNPARAM) {
    ranked <- .scaled_colranks_safe(data$exprs)

    overlord <- matrix(-Inf, 
        nrow=nrow(ranked), ncol=length(all.labels),
        dimnames=list(rownames(ranked), all.labels))

    for (i in names(data$references)) {
        rout <- .scaled_colranks_safe(data$references[[i]])
        idx <- buildIndex(rout, BNPARAM=BNPARAM)
        overlord[,i] <- .find_nearest_quantile(ranked, idx, quantile=quantile)
    }
    overlord
}

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
.fine_tune_de <- function(exprs, scores, references, quantile, tune.thresh, de.info, 
    BPPARAM=SerialParam(), BNPARAM=KmknnParam()) 
{
    deFUN <- function(labels) {
        targets <- lapply(de.info[labels], function(y) unique(unlist(y[labels], use.names=FALSE)))
        Reduce(union, targets)
    }

    .fine_tune_loop(exprs, scores=scores, references=references,
        quantile=quantile, tune.thresh=tune.thresh, defineFUN=deFUN,
        BPPARAM=BPPARAM, BNPARAM=BNPARAM)
}

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom DelayedMatrixStats rowSds
.fine_tune_sd <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh, 
    BPPARAM=SerialParam(), BNPARAM=KmknnParam()) 
{
    sdFUN <- function(labels) {
        curmat <- median.mat[,labels,drop=FALSE]
        rownames(median.mat)[rowSds(mat) >= sd.thresh]
    }

    .fine_tune_loop(exprs, scores=scores, references=references,
        quantile=quantile, tune.thresh=tune.thresh, defineFUN=sdFUN,
        BPPARAM=BPPARAM, BNPARAM=BNPARAM)
}
