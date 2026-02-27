#' Get classic markers
#'
#' Find markers between pairs of labels using the \dQuote{classic} approach,
#' i.e., based on the log-fold change between the medians of labels.
#'
#' @inheritParams trainSingleR
#' @param de.n An integer scalar specifying the number of DE genes to use.
#' Defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param BPPARAM Deprecated, use \code{num.threads} instead.
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
#' ref <- scrapper::normalizeRnaCounts.se(ref)
#' out <- getClassicMarkers(ref, labels=ref$label)
#' str(out)
#'
#' # Works with multiple references:
#' ref2 <- .mockRefData()
#' ref2 <- scrapper::normalizeRnaCounts.se(ref2)
#' out2 <- getClassicMarkers(list(ref, ref2), labels=list(ref$label, ref2$label))
#' str(out2)
#'
#' @export
#' @importFrom S4Vectors selfmatch DataFrame
#' @importFrom BiocGenerics cbind
#' @importFrom utils relist
#' @importFrom beachmat initializeCpp
getClassicMarkers <- function(ref, labels, assay.type="logcounts", check.missing=TRUE, de.n=NULL, num.threads=1, BPPARAM=NULL) { 
    num.threads <- .get_num_threads(num.threads, BPPARAM)

    if (!.is_list(ref)) { 
        ref <- list(ref)
    } else {
        labels <- unlist(labels)
    }

    for (i in seq_along(ref)) {
        ref[[i]] <- .to_clean_matrix(ref[[i]], assay.type, check.missing, msg="ref", num.threads=num.threads)
    }

    common <- Reduce(intersect, lapply(ref, rownames))
    if (length(common)==0L && any(vapply(ref, nrow, 0L) > 0L)) {
        stop("no common row names across 'ref'")
    }
    common <- as.character(common) # avoid problems with NULL rownames for zero-row inputs.
    for (i in seq_along(ref)) {
        # Use match() as this works with zero-row matrices that aren't allowed to have rownames,
        # see discussion at https://mailman.stat.ethz.ch/pipermail/r-devel/2006-August/038893.html.
        curmat <- ref[[i]]
        ref[[i]] <- DelayedArray(curmat)[match(common, rownames(curmat)),,drop=FALSE]
    }

    blocks <- NULL
    if (length(ref) > 1L) {
        blocks <- rep(seq_along(ref) - 1L, vapply(ref, ncol, FUN.VALUE=0L))
    }

    ref <- do.call(cbind, ref)
    ulabels <- .get_levels(labels)

    out <- find_classic_markers(
        initializeCpp(ref),
        length(ulabels),
        match(labels, ulabels) - 1L,
        blocks,
        de_n=de.n,
        nthreads=num.threads
    )

    names(out) <- ulabels
    for (i in seq_along(out)) {
        names(out[[i]]) <- ulabels
    }

    relist(common[unlist(out)], out)
}
