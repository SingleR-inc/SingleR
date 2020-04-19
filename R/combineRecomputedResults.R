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
#' @details
#' This function implements a variant of Option 3 described in \code{?"\link{combine-predictions}"}.
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
#' test <- scater::logNormCounts(test)
#'
#' ref1 <- scater::logNormCounts(ref1)
#' train1 <- trainSingleR(ref1, labels=ref1$label)
#' pred1 <- classifySingleR(test, train1)
#'
#' ref2 <- scater::logNormCounts(ref2)
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
#' @importFrom S4Vectors DataFrame 
#' @importFrom BiocParallel bpiterate SerialParam
#' @importFrom BiocNeighbors KmknnParam buildIndex
combineRecomputedResults <- function(results, test, trained, quantile=0.8, 
    assay.type.test="logcounts", check.missing=TRUE,
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
    refnames <- lapply(trained, function(x) rownames(x$original.exprs[[1]]))
    available <- c(list(rownames(test)), refnames)
    if (length(unique(available)) != 1) {
        warning("'test' and 'trained' differ in the universe of available genes")
    }
    available <- Reduce(intersect, available)

    markers <- vector("list", length(trained)) 
    for (i in seq_along(trained)) {
        current <- trained[[i]]

        if (current$search$mode=="de") {
            tmp <- lapply(current$search$extra, unlist, use.names=FALSE)
            tmp <- lapply(tmp, unique)
            stopifnot(identical(names(tmp), names(current$original.exprs)))
            tmp <- lapply(tmp, intersect, y=available)
            markers[[i]] <- tmp

        } else {
            nlabs <- length(current$original.exprs)
            tmp <- intersect(current$common.genes, available)
            markers[[i]] <- rep(list(tmp), nlabs)
        }
    }

    # Preparing expression values: subsetting them down to 
    # improve cache efficiency later on.
    test <- .to_clean_matrix(test, assay.type=assay.type.test, check.missing=check.missing, msg="test")
    all.ref <- lapply(trained, function(x) as.list(x$original.exprs))

    universe <- unique(unlist(markers))
    test <- test[universe,,drop=FALSE]
    for (i in seq_along(all.ref)) {
        current <- all.ref[[i]]
        for (j in seq_along(current)) {
            current[[j]] <- current[[j]][universe,,drop=FALSE]
        }
        all.ref[[i]] <- current
    }

    # Preparing genes (part 2).
    for (i in seq_along(markers)) {
        current <- markers[[i]]
        for (j in seq_along(current)) {
            current[[j]] <- match(current[[j]], universe)
        }
        markers[[i]] <- current
    }

    # Preparing labels.
    collated <- mapply(results, all.ref, FUN=function(x, y) {
        match(x$labels, names(y))
    }, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    all.labels <- do.call(rbind, collated)
    stopifnot(!any(is.na(all.labels)))

    ##############################################

    M <- .prep_for_parallel(test, BPPARAM, use.grid=is(test, "DelayedMatrix"))
    S <- .cofragment_matrix(M, all.labels)

    bp.out <- bpmapply(exprs=M, labels=S, FUN=.nonred_recompute_scores,
        MoreArgs=list(all.ref=all.ref, markers=markers, quantile=quantile),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
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
    chosen <- max.col(scores, ties.method="first")
    cbind(output, .combine_result_frames(chosen, results))
}

#' @importFrom S4Vectors DataFrame selfmatch
#' @importFrom BiocNeighbors buildIndex KmknnParam
.nonred_recompute_scores <- function(exprs, labels, all.ref, markers, quantile, BNPARAM=KmknnParam()) 
# This is a non-redundant calculator of the recomputed scores, where all label
# combinations that appear >= minimum times are used in a loop to avoid
# recalculation of ranked references in the corresponding C++ code. However, if
# it only appears once, then it is used directly in the C++ code for fast
# looping. Note that the realization of DA's requires sufficiently small 'exprs'.
{
    if (is(exprs, "DelayedMatrix")) {
        exprs <- as.matrix(exprs)
    }

    # Finding groups of cells with the same combination of per-reference assigned labels.
    collated <- DataFrame(t(labels))
    groups <- selfmatch(collated)
    by.group <- split(seq_along(groups), groups)
    above.min <- lengths(by.group) >= getOption("SingleR.recompute.minimum", 2L)
    scores <- matrix(0, ncol(exprs), ncol(collated))

    for (i in which(above.min)) {
        curgroup <- by.group[[i]]
        curlabels <- collated[curgroup[1],,drop=FALSE]

        curmarkers <- integer(0)
        for (i in seq_along(curlabels)) {
            curlab <- curlabels[[i]]
            curmarkers <- union(curmarkers, markers[[i]][[curlab]])
        }

        curtest <- exprs[curmarkers,curgroup,drop=FALSE]
        curtest <- .scaled_colranks_safe(curtest)

        for (i in seq_along(curlabels)) {
            curlab <- curlabels[[i]]
            curorig <- all.ref[[i]][[curlab]]
            curref <- curorig[curmarkers,,drop=FALSE]
            curref <- .scaled_colranks_safe(curref)
            curdex <- buildIndex(curref, BNPARAM=BNPARAM)
            scores[curgroup,i] <- .find_nearest_quantile(ranked=curtest, index=curdex, quantile=quantile)
        }
    }

    if (any(!above.min)) {
        affected <- sort(unlist(by.group[!above.min]))
        markers.m1 <- relist(unlist(markers) - 1, markers)

        stuff <- recompute_scores(
            Exprs=exprs[,affected,drop=FALSE],
            Labels=labels[,affected,drop=FALSE] - 1L,
            References=all.ref,
            Genes=markers.m1,
            quantile=quantile
        )

        scores[affected,] <- t(stuff)
    }

    scores
}
