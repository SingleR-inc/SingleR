#' @importFrom BiocGenerics rbind 
#' @importFrom BiocParallel SerialParam bpiterate
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocNeighbors KmknnParam
.fine_tune_loop <- function(exprs, scores, references, 
    quantile, tune.thresh, defineFUN, 
    BPPARAM=SerialParam(), BNPARAM=KmknnParam())
{
    indices <- list()
    present <- DataFrame(scores[0,])
    final <- character(nrow(scores))
    all.labels <- colnames(scores)
    survivors <- seq_len(nrow(scores))

    repeat {
        max.score <- scores[cbind(seq_len(nrow(scores)), max.col(scores))]
        limit <- max.score - tune.thresh
        above <- scores >= limit
        encoding <- DataFrame(overlord)
        nhits <- Reduce("+", encoding)
        finished <- nhits==1L

        final[survivors[finished]] <- "X" # TODO: fix this.
        survivors <- survivors[!finished]
        encoding <- encoding[!finished,,drop=FALSE]
        if (length(survivors)==0L) {
            break
        }

        i.out <- .build_indices(encoding, present=present, indices=indices, defineFUN=defineFUN, BNPARAM=BNPARAM)
        present <- i.out$present
        indices <- i.out$indices
        use <- i.out$use

        ITER <- .fine_tune_iterator(use, present=present, indices=indices, exprs=exprs) 
        output <- bpiterate(ITER=ITER, FUN=.fine_tune_searcher, BPPARAM=BPPARAM, quantile=quantile, all.labels=all.labels)
        scores <- do.call(rbind, output)
        scores[order(use),] <- scores # need to reorder as iterator splits up by 'use'.
    }

    final
}

#' @importFrom BiocNeighbors buildIndex KmknnParam
#' @importFrom S4Vectors selfmatch
.build_indices <- function(encoding, present, indices, defineFUN, BNPARAM=KmknnParam()) {
    smatched <- selfmatch(encoding)
    representatives <- encoding[!duplicated(matched),,drop=FALSE]
    imatched <- match(representatives, present)
    absent <- which(is.na(imatched))

    absent.indices <- vector("list", length(absent))
    for (i in seq_along(absent)) {
        chosen <- representatives[absent[i],]
        chosen <- colnames(scores)[unlist(chosen)]
        common <- defineFUN(chosen)

        idxs <- lapply(references[chosen], function(x) {
            rout <- .scaled_colranks_safe(x[common,,drop=FALSE])
            buildIndex(rout, BNPARAM=BNPARAM)
        })
        absent.indices[[i]] <- list(trained=idxs, genes=common)
    }

    indices <- c(indices, absent.indices)
    present <- rbind(present, representatives[absent,])
    imatched[absent] <- length(indices) + seq_along(absent.indices)

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

    overlord <- matrix(-Inf, nrow=nrow(ranked), ncol=length(all.labels),
        dimnames=list(rownames(ranked), all.labels))
    for (i in names(output)) {
        overlord[,i] <- output[[i]]
    }
    overlord
}

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

    do.call(mapply, c(bp.out, list(FUN=c, SIMPLIFY=FALSE, USE.NAMES=FALSE)))
}

#' @importFrom BiocParallel bpmapply SerialParam
.fine_tune_sd <- function(exprs, scores, references, quantile, tune.thresh, median.mat, sd.thresh, BPPARAM=SerialParam()) {
    M <- .prep_for_parallel(exprs, BPPARAM)
    S <- .prep_for_parallel(t(scores), BPPARAM)
    bp.out <- bpmapply(Exprs=M, scores=S, FUN=fine_tune_label_sd, 
        MoreArgs=list(References=references, quantile=quantile, tune_thresh=tune.thresh, median_mat=t(median.mat), sd_thresh=sd.thresh),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    do.call(mapply, c(bp.out, list(FUN=c, SIMPLIFY=FALSE, USE.NAMES=FALSE)))
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
