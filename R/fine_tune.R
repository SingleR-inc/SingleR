#' @importFrom BiocGenerics rbind 
#' @importFrom BiocParallel SerialParam bpiterate
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom DelayedArray rowMaxs
.fine_tune_loop <- function(exprs, scores, references, 
    quantile, tune.thresh, defineFUN, 
    BPPARAM=SerialParam(), BNPARAM=KmknnParam())
{
    indices <- list()
    present <- DataFrame(scores[0,])
    all.labels <- colnames(scores)

    final <- character(nrow(scores))
    survivors <- seq_len(nrow(scores))
    remaining <- exprs
    last.nhits <- integer(nrow(scores))

    repeat {
        above <- scores >= rowMaxs(scores) - tune.thresh
        nhits <- rowSums(above)
        finished <- nhits==1L | nhits==last.nhits[survivors]
        last.nhits[survivors] <- nhits

        final[survivors[finished]] <- "X" # TODO: fix this.

        if (all(finished)) {
           break
        } else if (any(finished)) {
            survivors <- survivors[!finished]
            above <- above[!finished,,drop=FALSE]
            remaining <- remaining[,!finished,drop=FALSE]
        }

        i.out <- .build_indices(above, present=present, indices=indices, 
            references=references, defineFUN=defineFUN, BNPARAM=BNPARAM)
        present <- i.out$present
        indices <- i.out$indices
        use <- i.out$use

        ITER <- .fine_tune_iterator(use, present=present, indices=indices, exprs=remaining)
        output <- bpiterate(ITER=ITER, FUN=.fine_tune_searcher, BPPARAM=BPPARAM, quantile=quantile, all.labels=all.labels)
        scores <- do.call(rbind, output)
        scores[order(use),] <- scores # need to reorder as iterator splits up by 'use'.
    }

    final
}

#' @importFrom BiocNeighbors buildIndex KmknnParam
#' @importFrom S4Vectors DataFrame selfmatch
#' @importFrom BiocGenerics match
.build_indices <- function(encoding, present, indices, references, defineFUN, BNPARAM=KmknnParam()) {
    encoding2 <- DataFrame(encoding)
    smatched <- selfmatch(encoding2)
    smatched <- as.integer(factor(smatched))

    representatives <- encoding2[!duplicated(smatched),,drop=FALSE]
    imatched <- match(representatives, present)
    absent <- which(is.na(imatched))

    absent.indices <- vector("list", length(absent))
    for (i in seq_along(absent)) {
        chosen <- encoding[absent[i],]
        chosen <- colnames(encoding)[chosen]
        common <- defineFUN(chosen)

        idxs <- lapply(references[chosen], function(x) {
            rout <- .scaled_colranks_safe(x[common,,drop=FALSE])
            buildIndex(rout, BNPARAM=BNPARAM)
        })
        absent.indices[[i]] <- list(trained=idxs, genes=common)
    }

    imatched[absent] <- length(indices) + seq_along(absent.indices)
    indices <- c(indices, absent.indices)
    present <- rbind(present, representatives[absent,])

    list(
        present=present,
        indices=indices,
        use=imatched[smatched]
    )
}

.fine_tune_iterator <- function(use, present, indices, exprs) {
    u.use <- sort(unique(use))
    counter <- 1L

    function() {
        if (counter > length(u.use)) {
            return(NULL)
        }

        in.use <- u.use[counter]
        used <- use==in.use
        curdex <- indices[[in.use]]

        counter <<- counter + 1L
        list(
            exprs=exprs[curdex$genes,used,drop=FALSE],
            trained=curdex$trained
        )
    }
}

.fine_tune_searcher <- function(data, quantile, all.labels) {
    exprs <- data$exprs
    trained <- data$trained
    ranked <- .scaled_colranks_safe(exprs)
    output <- lapply(trained, FUN=.find_nearest_quantile, ranked=ranked, quantile=quantile)

    overlord <- matrix(-Inf, 
        nrow=nrow(ranked), ncol=length(all.labels),
        dimnames=list(rownames(ranked), all.labels))

    for (i in names(output)) {
        overlord[,i] <- output[[i]]
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
