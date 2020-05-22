#' Get classic markers
#'
#' Find markers between pairs of labels using the \dQuote{classic} approach,
#' i.e., based on the log-fold change between the medians of labels.
#'
#' @inheritParams trainSingleR
#' @param de.n An integer scalar specifying the number of DE genes to use.
#' Defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
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
#' ref <- scater::logNormCounts(ref)
#' out <- getClassicMarkers(ref, labels=ref$label)
#' str(out)
#'
#' # Works with multiple references:
#' ref2 <- .mockRefData()
#' ref2 <- scater::logNormCounts(ref2)
#' out2 <- getClassicMarkers(list(ref, ref2), labels=list(ref$label, ref2$label))
#' str(out2)
#'
#' @export
#' @importFrom S4Vectors selfmatch
getClassicMarkers <- function(ref, labels, assay.type="logcounts", check.missing=TRUE, de.n=NULL) { 

    if (!.is_list(ref)) {
        ulabels <- .get_levels(labels)
        stats <- .pairwise_median_logfc(ref, labels, assay.type=assay.type, check.missing=check.missing)
        choices <- stats$choices
        lfc <- stats$lfc
    } else {
        ulabels <- .get_levels(unlist(labels))
        stats <- mapply(FUN=.pairwise_median_logfc, ref=ref, labels=labels, 
            assay.type=assay.type, check.missing=check.missing, SIMPLIFY=FALSE)

        all.choices <- lapply(stats, "[[", i="choices")
        all.choices <- do.call(rbind, all.choices)
        all.lfc <- lapply(stats, "[[", i="lfc")
        all.lfc <- unlist(all.lfc, recursive=FALSE)

        m <- selfmatch(all.choices)
        by.m <- split(all.lfc, m)
        lfc <- lapply(by.m, Reduce, f="+")
        choices <- all.choices[match(names(lfc), as.character(m)),]
    }

    .choose_top_markers(ulabels, choices=choices, lfc=lfc, de.n=de.n)
}

.pairwise_median_logfc <- function(ref, labels, assay.type, check.missing) {
    ref <- .to_clean_matrix(ref, assay.type, check.missing, msg="ref")
    mat <- .median_by_label(ref, labels)

    choices <- expand.grid(first=colnames(mat), second=colnames(mat), stringsAsFactors=FALSE)
    choices <- choices[choices$first!=choices$second,]
    choices <- DataFrame(choices)

    lfc <- vector("list", nrow(choices))
    for (i in seq_along(lfc)) {
        left <- choices$first[i]
        right <- choices$second[i]
        lfc[[i]] <- mat[,left] - mat[,right]
    }

    list(choices=choices, lfc=lfc)
}

#' @importFrom utils head
.choose_top_markers <- function(ulabels, choices, lfc, de.n) {
    if (is.null(de.n)) {
        de.n <- round(500*(2/3)^log2(length(ulabels)))
    }

    output <- list()
    for (i in ulabels) {
        subout <- list()
        for (j in ulabels) {
            subout[[j]] <- character(0)
        }
        output[[i]] <- subout
    }

    for (i in seq_len(nrow(choices))) {
        left <- choices$first[i]
        right <- choices$second[i]
        chosen <- lfc[[i]]

        # FYI, the as.character() accounts for the edge case of NULL names,
        # which apparently happens due to incorrect matrix zero-subsetting.
        output[[left]][[right]] <- as.character(names(chosen)[head(order(chosen, decreasing=TRUE), de.n)])
    }

    output
}
