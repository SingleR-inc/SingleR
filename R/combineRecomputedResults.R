#' Combine SingleR results with recomputation
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} (usually against different references) into a single \linkS4class{DataFrame}.
#' The label from the results with the highest score for each cell is retained.
#' Unlike \code{\link{combineCommonResults}}, this does not assume that each run of \code{\link{classifySingleR}} was performed using the same set of common genes, instead recomputing the scores for comparison across references.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} when run on each reference separately.
#' @inheritParams SingleR
#' @param trained A list of \linkS4class{List}s containing the trained outputs of multiple references,
#' equivalent to either (i) the output of \code{\link{trainSingleR}} on multiple references with \code{recompute=TRUE},
#' or (ii) running \code{trainSingleR} on each reference separately and manually making a list of the trained outputs.
#'
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This mimics the output of \code{\link{classifySingleR}} and contains the following fields:
#' \itemize{
#' \item \code{scores}, a numeric matrix of correlations containing the \emph{recomputed} scores.
#' For any given cell, entries of this matrix are only non-\code{NA} for the assigned label in each reference;
#' scores are not recomputed for the other labels.
#' \item \code{labels}, a character vector containing the per-cell combined label across references.
#' \item \code{references}, an integer vector specifying the reference from which the combined label was derived.
#' \item \code{orig.results}, a DataFrame containing \code{results}.
#' }
#' It may also contain \code{first.labels} and \code{pruned.labels} if these were also present in \code{results}.
#'
#' The \code{\link{metadata}} contains \code{label.origin}, 
#' a DataFrame specifying the reference of origin for each label in \code{scores}.
#' Note that, unlike \code{\link{combineCommonResults}}, no \code{common.genes} is reported
#' as this function does not use a common set of genes across all references.
#'
#' @details
#' Here, the strategy is to performed classification separately within each reference, 
#' then collating the results to choose the label with the highest score across references.
#' For a given cell in \code{test}, we extract its assigned label from \code{results} for each reference.
#' We also retrieve the marker genes associated with that label and take the union of markers across all references.
#' This defines a common feature space in which the score for each reference's assigned label is recomputed using \code{ref};
#' the label from the reference with the top recomputed score is then reported as the combined annotation for that cell.
#' 
#' Unlike \code{\link{combineCommonResults}}, the union of markers is not used for the within-reference calls.
#' This avoids the inclusion of noise from irrelevant genes in the within-reference assignments.
#' Obviously, \code{combineRecomputedResults} is slower as it does require recomputation of the scores,
#' but the within-reference calls are faster as there are fewer genes in the union of markers for assigned labels
#' (compared to the union of markers across all labels, as required by \code{\link{combineCommonResults}}),
#' so it is likely that the net compute time should be lower.
#'
#' It is strongly recommended that the universe of genes be the same across all references.
#' The intersection of genes across all \code{ref} and \code{test} is used when recomputing scores,
#' and differences in the availability of genes between references may have unpredictable effects.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{SingleR}} and \code{\link{classifySingleR}}, for generating predictions to use in \code{results}.
#'
#' \code{\link{combineCommonResults}}, for another approach to combining predictions.
#'
#' @references
#' Lun A, Bunis D, Andrews J (2020).
#' Thoughts on a more scalable algorithm for multiple references.
#' \url{https://github.com/LTLA/SingleR/issues/94}
#'
#' @examples
#' # Making up data.
#' ref <- .mockRefData(nreps=8)
#' ref1 <- ref[,1:2%%2==0]
#' ref2 <- ref[,1:2%%2==1]
#' ref2$label <- tolower(ref2$label)
#'
#' test <- .mockTestData(ref)
#'
#' # Performing classification within each reference.
#' test <- scuttle::logNormCounts(test)
#'
#' ref1 <- scuttle::logNormCounts(ref1)
#' train1 <- trainSingleR(ref1, labels=ref1$label)
#' pred1 <- classifySingleR(test, train1)
#'
#' ref2 <- scuttle::logNormCounts(ref2)
#' train2 <- trainSingleR(ref2, labels=ref2$label)
#' pred2 <- classifySingleR(test, train2)
#'
#' # Combining results with recomputation of scores.
#' combined <- combineRecomputedResults(
#'     results=list(pred1, pred2), 
#'     test=test,
#'     trained=list(train1, train2))
#' 
#' combined[,1:5]
#'
#' @export
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam buildIndex
#' @importFrom beachmat colBlockApply
combineRecomputedResults <- function(results, test, trained, quantile=0.8, 
    assay.type.test="logcounts", check.missing=TRUE, allow.lost=FALSE, warn.lost=TRUE,
    BNPARAM=KmknnParam(), BPPARAM=SerialParam())
{
    all.names <- c(list(colnames(test)), lapply(results, rownames))
    if (length(unique(all.names)) != 1) {
        stop("cell/cluster names in 'results' are not identical")
    }
    all.nrow <- c(ncol(test), vapply(results, nrow, 0L))
    if (length(unique(all.nrow)) != 1) {
        stop("numbers of cells/clusters in 'results' are not identical")
    }

    ##############################################

    # Preparing genes (part 1).
    markers <- vector("list", length(trained)) 
    for (i in seq_along(trained)) {
        current <- trained[[i]]
        if (current$search$mode=="de") {
            tmp <- lapply(current$search$extra, unlist, use.names=FALSE)
            tmp <- lapply(tmp, unique)
            stopifnot(identical(names(tmp), names(current$original.exprs)))
            markers[[i]] <- tmp
        } else {
            nlabs <- length(current$original.exprs)
            tmp <- current$common.genes
            markers[[i]] <- rep(list(tmp), nlabs)
        }
    }

    universe <- unique(unlist(markers))
    if (!all(universe %in% rownames(test))) {
        stop("all markers stored in 'results' should be present in 'test'")
    }

    refnames <- lapply(trained, function(x) rownames(x$original.exprs[[1]]))
    refnames <- Reduce(intersect, refnames)
    if (has.lost <- !all(universe %in% refnames)) {
        if (warn.lost) {
            stop("not all markers are present across all 'results'")
        }
        if (!allow.lost) {
            universe <- intersect(universe, refnames)
            for (i in seq_along(markers)) {
                markers[[i]] <- lapply(markers[[i]], intersect, y=refnames)
            }
        }
    }

    ##############################################

    # Preparing expression values: subsetting them down to 
    # improve cache efficiency later on.
    test <- .to_clean_matrix(test, assay.type=assay.type.test, check.missing=check.missing, msg="test")
    test <- test[universe,,drop=FALSE]

    all.ref <- vector("list", length(trained))
    for (i in seq_along(all.ref)) {
        current <- trained[[i]]$original.exprs
        current <- lapply(current, FUN=.subset_with_NAs, universe=universe)
        all.ref[[i]] <- current
    }

    # Preparing genes (part 2).
    for (i in seq_along(markers)) {
        current <- markers[[i]]
        for (j in seq_along(current)) {
            current[[j]] <- match(current[[j]], universe) - 1L
        }
        markers[[i]] <- current
    }

    # Preparing labels.
    collated <- vector("list", length(results))
    for (i in seq_along(collated)) {
        collated[[i]] <- match(results[[i]]$labels, names(all.ref[[i]]))
    }
    all.labels <- do.call(rbind, collated)
    stopifnot(!any(is.na(all.labels)))

    ##############################################

    bp.out <- colBlockApply(test, FUN=.nonred_recompute_scores, BPPARAM=BPPARAM,
        labels=all.labels, all.ref=all.ref, markers=markers, quantile=quantile, has.lost=has.lost)
    scores <- do.call(rbind, bp.out)

    base.scores <- vector("list", length(results))
    for (r in seq_along(base.scores)) {
        mat <- results[[r]]$scores   
        mat[] <- NA_real_
        idx <- cbind(seq_len(nrow(mat)), all.labels[r,]) 
        mat[idx] <- scores[,r]
        base.scores[[r]] <- mat
    }

    all.scores <- do.call(cbind, base.scores)
    output <- DataFrame(scores = I(all.scores), row.names=rownames(results[[1]]))
    metadata(output)$label.origin <- .create_label_origin(base.scores)

    chosen <- max.col(scores, ties.method="first")
    cbind(output, .combine_result_frames(chosen, results))
}

#' @importFrom S4Vectors DataFrame selfmatch
#' @importFrom BiocNeighbors buildIndex KmknnParam
#' @importFrom DelayedArray currentViewport makeNindexFromArrayViewport
.nonred_recompute_scores <- function(block, labels, all.ref, markers, quantile, has.lost) {
    vp <- currentViewport()
    idx <- makeNindexFromArrayViewport(vp, expand.RangeNSBS = TRUE)[[2]]
    if (!is.null(idx)) {
        labels <- labels[,idx,drop=FALSE]
    }

    # Finding groups of cells with the same combination of per-reference assigned labels.
    collated <- DataFrame(t(labels))
    groups <- selfmatch(collated)
    by.group <- split(seq_along(groups) - 1L, groups)

    if (has.lost) {
        RECOMPFUN <- recompute_scores_with_na
    } else {
        RECOMPFUN <- recompute_scores
    }

    scores <- RECOMPFUN(
        Groups=by.group,
        Exprs=block,
        Labels=labels - 1L,
        References=all.ref,
        Genes=markers,
        quantile=quantile
    )

    t(scores)
}

.subset_with_NAs <- function(x, universe) {
    m <- match(universe, rownames(x))
    if (!anyNA(m)) {
        .realize_reference(x[m,,drop=FALSE])
    } else if (nrow(x)==0) {
        # Impossible in practice, but just in case...
        matrix(NA_real_, length(universe), ncol(x), dimnames=list(universe, colnames(x)))
    } else {
        # A little song and dance because some matrix formats don't take well to NA indices.
        lost <- is.na(m)
        m[lost] <- 1L
        output <- x[m,,drop=FALSE]
        output <- .realize_reference(output)
        output[lost,] <- NA_real_
        output
    }
}
