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
#' @param warn.lost Logical scalar indicating whether to emit a warning if markers from one reference in \code{trained} are \dQuote{lost} in other references.
#' @param allow.lost Logical scalar indicating whether to use lost markers in references where they are available. 
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
#' Here, the strategy is to perform classification separately within each reference, 
#' then collate the results to choose the label with the highest score across references.
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
#' @section Dealing with mismatching gene availabilities:
#' It is strongly recommended that the universe of genes be the same across all references in \code{trained}.
#' If this is not the case, the intersection of genes across all \code{trained} will be used in the recomputation.
#' This at least provides a common feature space for comparing correlations, 
#' though differences in the availability of markers between references may have unpredictable effects on the results
#' (and so a warning will be emitted by default, when when \code{warn.lost=TRUE}).
#'
#' That said, the intersection may be too string when dealing with many references with diverse feature annotations. 
#' In such cases, we can set \code{allow.lost=TRUE} so that the recomputation for each reference will use all available markers in that reference.
#' The idea here is to avoid penalizing all references by removing an informative marker when it is only absent in a single reference.
#' We hope that the recomputed scores are still roughly comparable if the number of lost markers is relatively low,
#' coupled with the use of ranks in the calculation of the Spearman-based scores to reduce the influence of individual markers.
#' This is perhaps as reliable as one might imagine, so setting \code{allow.lost=TRUE} should be considered a last resort.
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
combineRecomputedResults <- function(
    results, 
    test, 
    trained, 
    quantile=0.8, 
    assay.type.test="logcounts", 
    check.missing=TRUE, 
    allow.lost=FALSE, 
    warn.lost=TRUE,
    num.threads=1,
    BPPARAM=SerialParam())
{
    all.names <- c(list(colnames(test)), lapply(results, rownames))
    if (length(unique(all.names)) != 1) {
        stop("cell/cluster names in 'results' are not identical")
    }
    all.nrow <- c(ncol(test), vapply(results, nrow, 0L))
    if (length(unique(all.nrow)) != 1) {
        stop("numbers of cells/clusters in 'results' are not identical")
    }

    # Applying the integration.
    universe <- Reduce(union, c(list(rownames(test)), lapply(trained, function(x) rownames(x$ref))))
    ibuilt <- integrate_build(
        match(rownames(test), universe) - 1L,
        lapply(trained, function(x) x$ref),
        lapply(trained, function(x) match(rownames(x$ref), universe) - 1L), 
        lapply(trained, function(x) match(x$labels$full, x$labels$unique) - 1L),
        lapply(trained, function(x) x$built),
        nthreads = num.threads
    )

    test <- .to_clean_matrix(test, assay.type=assay.type.test, check.missing=check.missing, msg="test", BPPARAM=BPPARAM)
    collated <- vector("list", length(trained))
    for (i in seq_along(collated)) {
        collated[[i]] <- match(results[[i]]$labels, trained[[i]]$labels$unique) - 1L
    }

    irun <- integrate_run(test, collated, ibuilt, nthreads = num.threads) 
    scores <- irun$scores

    # Organizing the outputs.
    base.scores <- vector("list", length(results))
    for (r in seq_along(base.scores)) {
        mat <- results[[r]]$scores   
        mat[] <- NA_real_
        idx <- cbind(seq_len(nrow(mat)), collated[[i]] + 1L)
        mat[idx] <- scores[,r]
        base.scores[[r]] <- mat
    }

    all.scores <- do.call(cbind, base.scores)
    output <- DataFrame(scores = I(all.scores), row.names=rownames(results[[1]]))
    metadata(output)$label.origin <- .create_label_origin(base.scores)

    chosen <- irun$best + 1L
    cbind(output, .combine_result_frames(chosen, results))
}

#' @importFrom S4Vectors DataFrame
.combine_result_frames <- function(chosen, results) {
    has.first <- !is.null(results[[1]]$first.labels)
    has.pruned <- !is.null(results[[1]]$pruned.labels)

    # Organizing the statistics based on the chosen results.
    chosen.label <- chosen.first <- chosen.pruned <- rep(NA_character_, nrow(results[[1]]))

    for (u in unique(chosen)) {
        current <- chosen==u
        res <- results[[u]]
        chosen.label[current] <- res$labels[current]

        if (has.first) { # assume that either everyone has 'first', or no-one does.
            chosen.first[current] <- res$first.labels[current]
        }

        if (has.pruned) { # same for pruned.
            chosen.pruned[current] <- res$pruned.labels[current]
        }
    }

    output <- DataFrame(labels=chosen.label, row.names=rownames(results[[1]]))

    if (has.first) {
        output$first.labels <- chosen.first
        output <- output[,c("first.labels", "labels"),drop=FALSE]
    }

    if (has.pruned) {
        output$pruned.labels <- chosen.pruned
    }

    output$reference <- chosen

    if (is.null(names(results))) {
        names(results) <- sprintf("ref%i", seq_along(results))
    }
    output$orig.results <- do.call(DataFrame, lapply(results, I))

    output
}

#' @importFrom S4Vectors DataFrame
.create_label_origin <- function(collected.scores) {
    DataFrame(
        label=unlist(lapply(collected.scores, colnames)),
        reference=rep(seq_along(collected.scores), vapply(collected.scores, ncol, 0L))
    )
}
