#' Get classic markers
#'
#' Find markers between pairs of labels using the \dQuote{classic} approach,
#' i.e., based on the log-fold change between the medians of labels.
#'
#' @inheritParams trainSingleR
#' @param de.n An integer scalar specifying the number of DE genes to use.
#' Defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#'
#' @return
#' A list of lists of character vectors, 
#' where both the outer and inner lists have names equal to the unique levels of \code{labels}.
#' The character vector contains the names of the top \code{de.n} genes with the largest positive log-fold changes
#' in one label (entry of the outer list) against another label (entry of the inner list).
#'
#' @details
#' This function implements the classic mode of marker detection in \pkg{SingleR},
#' based only on the magnitude of the log-fold change between labels.
#' In many respects, this approach may be suboptimal as it does not consider the variance within each label
#' and has limited precision when the expression values are highly discrete.
#' Nonetheless, it is often the only possible approach when dealing with reference datasets
#' that lack replication and thus cannot be used with more advanced marker detection methods.
#'
#' If multiple references are supplied, 
#' ranking is performed based on the average of the log-fold changes within each reference.
#' This avoids comparison of expression values across references that can be distorted by batch effects.
#' If a pair of labels does not co-occur in at least one reference,
#' no attempt is made to perform the comparison and the corresponding character vector is left empty in the output.
#'
#' The character vector corresponding to the comparison of a label to itself is always empty.
#' 
#' @author Aaron Lun, based on the original \code{SingleR} code by Dvir Aran.
#' @seealso
#' \code{\link{trainSingleR}} and \code{\link{SingleR}},
#' where this function is used when \code{genes="de"} and \code{de.method="classic"}.
#' 
#' @examples
#' ref <- .mockRefData()
#' ref <- scuttle::logNormCounts(ref)
#' out <- getClassicMarkers(ref, labels=ref$label)
#' str(out)
#'
#' # Works with multiple references:
#' ref2 <- .mockRefData()
#' ref2 <- scuttle::logNormCounts(ref2)
#' out2 <- getClassicMarkers(list(ref, ref2), labels=list(ref$label, ref2$label))
#' str(out2)
#'
#' @export
#' @importFrom S4Vectors selfmatch DataFrame
#' @importFrom BiocParallel SerialParam bpnworkers
getClassicMarkers <- function(ref, labels, assay.type="logcounts", check.missing=TRUE, de.n=NULL, num.threads=bpnworkers(BPPARAM), BPPARAM=SerialParam()) { 
    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (!.is_list(ref)) { 
        ref <- list(ref)
        labels <- list(labels)
    }

    # Setting up references.
    common <- Reduce(intersect, lapply(ref, rownames))
    if (length(common)==0L && any(vapply(ref, nrow, 0L) > 0L)) {
        stop("no common row names across 'ref'")
    }

    for (i in seq_along(ref)) {
        current <- ref[[i]][common,,drop=FALSE]
        current <- .to_clean_matrix(current, assay.type, check.missing, msg="ref", BPPARAM=BPPARAM)
        ref[[i]] <- .median_by_label(current, labels[[i]], BPPARAM=BPPARAM)
    }
    
    ulabels <- .get_levels(unlist(lapply(ref, colnames)))
    available <- list()
    for (i in seq_along(ref)) {
        m <- match(ulabels, colnames(ref[[i]]))
        available[[i]] <- !is.na(m)
        ref[[i]] <- ref[[i]][,m,drop=FALSE]
        colnames(ref[[i]]) <- ulabels
    }

    # Identify all label combinations within each reference.
    collated <- list()
    for (i in seq_along(ref)) {
        curavail <- which(available[[i]])
        pairs <- expand.grid(first=curavail, second=curavail)
        pairs <- pairs[pairs$first!=pairs$second,]
        collated[[i]] <- DataFrame(pairs)
    }
    collated <- do.call(rbind, collated)
    choices <- unique(collated)

    # Identify top hits based on the average (or sum, it doesn't matter)
    # of the log-fold changes between labels across references.
    if (is.null(de.n)) {
        de.n <- round(500*(2/3)^log2(length(ulabels)))
    }
    stopifnot(de.n > 0)

    out <- find_classic_markers(nlabels=length(ulabels), 
        ngenes=length(common), 
        left=choices$first - 1L, 
        right=choices$second - 1L, 
        ref=ref, 
        available=available, 
        de_n=de.n,
        nthreads=num.threads
    )

    names(out) <- ulabels
    for (i in seq_along(out)) {
        names(out[[i]]) <- ulabels
    }

    relist(common[unlist(out)], out)
}

#' @importFrom DelayedMatrixStats rowMedians
#' @importFrom DelayedArray DelayedArray
.median_by_label <- function(mat, labels, BPPARAM=SerialParam()) {
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    ulabels <- .get_levels(labels)
    output <- matrix(0, nrow(mat), length(ulabels))
    rownames(output) <- rownames(mat)
    colnames(output) <- ulabels

    for (u in ulabels) {
        # Disambiguate from Biobase::rowMedians.
        output[,u] <- DelayedMatrixStats::rowMedians(DelayedArray(mat), cols=u==labels)
    }
    output
}
