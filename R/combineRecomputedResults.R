#' Combine SingleR results with recomputation
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} (usually against different references) into a single \linkS4class{DataFrame}.
#' The label from the results with the highest score for each cell is retained.
#' Unlike \code{\link{combineCommonResults}}, this does not assume that each run of \code{\link{classifySingleR}} was performed using the same set of common genes, instead recomputing the scores for comparison across references.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} when run on each reference separately.
#' @inheritParams SingleR
#' @param check.missing Deprecated and ignored, as any row filtering will cause mismatches with the \code{test.genes=} used in \code{\link{trainSingleR}}.
#' @param trained A list of \linkS4class{List}s containing the trained outputs of multiple references,
#' equivalent to either (i) the output of \code{\link{trainSingleR}} on multiple references with \code{recompute=TRUE},
#' or (ii) running \code{trainSingleR} on each reference separately and manually making a list of the trained outputs.
#' @param warn.lost Logical scalar indicating whether to emit a warning if markers from one reference in \code{trained} are absent in other references.
#' @param quantile Numeric scalar specifying the quantile of the correlation distribution to use for computing the score, see \code{\link{classifySingleR}}.
#' @param allow.lost Deprecated.
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
#' It may also contain \code{pruned.labels} if these were also present in \code{results}.
#'
#' The \code{\link{metadata}} contains \code{label.origin}, 
#' a DataFrame specifying the reference of origin for each label in \code{scores}.
#'
#' @details
#' Here, the strategy is to perform classification separately within each reference, 
#' then collate the results to choose the label with the highest score across references.
#' For a given cell in \code{test}, we extract its assigned label from \code{results} for each reference.
#' We also retrieve the marker genes associated with that label and take the union of markers across all references.
#' This defines a common feature space in which the score for each reference's assigned label is recomputed using \code{ref};
#' the label from the reference with the top recomputed score is then reported as the combined annotation for that cell.
#' 
#' A key aspect of this approach is that each entry of \code{results} is generated with reference-specific markers.
#' This avoids the inclusion of noise from irrelevant genes in the within-reference assignments.
#' Similarly, the common feature space for each cell is defined from the most relevant markers across all references,
#' analogous to one iteration of fine-tuning using only the best labels from each reference.
#' Compare this to the alternative approach of creating a common feature space, where we force all per-reference classifications to use the same set of markers;
#' this would slow down each individual classification as many more genes are involved.
#'
#' @section Dealing with mismatching gene availabilities:
#' It is recommended that the universe of genes be the same across all references in \code{trained}.
#' (Or, at the very least, markers used in one reference are available in the others.)
#' This ensures that a common feature space can be generated when comparing correlations across references.
#' Differences in the availability of markers between references will have unpredictable effects on the comparability of correlation scores,
#' so a warning will be emitted by default when \code{warn.lost=TRUE}.
#' Callers can protect against this by subsetting each reference to the intersection of features present across all references - this is done by default in \code{\link{SingleR}}.
#'
#' That said, this requirement may be too strict when dealing with many references with diverse feature annotations. 
#' In such cases, the recomputation for each reference will automatically use all available markers in that reference.
#' The idea here is to avoid penalizing all references by removing an informative marker when it is only absent in a single reference.
#' We hope that the recomputed scores are still roughly comparable if the number of lost markers is relatively low,
#' coupled with the use of ranks in the calculation of the Spearman-based scores to reduce the influence of individual markers.
#' This is perhaps as reliable as one might imagine.
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
#' @importFrom beachmat initializeCpp
combineRecomputedResults <- function(
    results, 
    test, 
    trained, 
    quantile=0.8, 
    assay.type.test="logcounts", 
    check.missing=FALSE,
    warn.lost=TRUE,
    allow.lost=FALSE, 
    num.threads = bpnworkers(BPPARAM),
    BPPARAM=SerialParam())
{
    test <- .to_clean_matrix(test, assay.type=assay.type.test, check.missing=FALSE, msg="test", BPPARAM=BPPARAM)

    # Applying the sanity checks.
    stopifnot(length(results) == length(trained))
    for (i in seq_along(results)) {
        curres <- results[[i]]
        if (ncol(test) != nrow(curres)) {
            stop("numbers of cells/clusters in 'results' are not identical")
        }
        if (!identical(rownames(curres), colnames(test))) {
            stop("cell/cluster names in 'results' are not identical")
        }

        curtrain <- trained[[i]]
        if (!all(curres$labels %in% curtrain$labels$unique)) {
            stop("not all labels in 'results' are present in 'trained'")
        }
        .check_test_genes(test, curtrain)
    }

    # Checking the genes.
    all.refnames <- lapply(trained, function(x) rownames(x$ref))
    if (warn.lost) {
        intersected <- Reduce(intersect, all.refnames)
        for (i in seq_along(trained)) {
            if (!all(trained[[i]]$markers$unique %in% intersected)) {
                warning("not all markers in 'trained' are available in each reference")
            }
        }
    }

    # Applying the integration.
    universe <- Reduce(union, c(list(rownames(test)), all.refnames))
    ibuilt <- train_integrated(
        test_features=match(rownames(test), universe) - 1L,
        references=lapply(trained, function(x) initializeCpp(x$ref)),
        ref_ids=lapply(all.refnames, function(x) match(x, universe) - 1L), 
        labels=lapply(trained, function(x) match(x$labels$full, x$labels$unique) - 1L),
        prebuilt=lapply(trained, function(x) rebuildIndex(x)$built),
        nthreads = num.threads
    )

    collated <- vector("list", length(trained))
    for (i in seq_along(collated)) {
        collated[[i]] <- match(results[[i]]$labels, trained[[i]]$labels$unique) - 1L
    }

    parsed <- initializeCpp(test)
    irun <- classify_integrated(
        test=parsed,
        results=collated,
        integrated_build=ibuilt,
        quantile=quantile,
        nthreads=num.threads
    ) 

    # Organizing the outputs.
    base.scores <- vector("list", length(results))
    for (r in seq_along(base.scores)) {
        mat <- results[[r]]$scores   
        mat[] <- NA_real_
        idx <- cbind(seq_len(nrow(mat)), collated[[r]] + 1L)
        mat[idx] <- irun$scores[,r]
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
    has.pruned <- !is.null(results[[1]]$pruned.labels)

    # Organizing the statistics based on the chosen results.
    chosen.label <- chosen.first <- chosen.pruned <- rep(NA_character_, nrow(results[[1]]))

    for (u in unique(chosen)) {
        current <- chosen==u
        res <- results[[u]]
        chosen.label[current] <- res$labels[current]

        if (has.pruned) { # same for pruned.
            chosen.pruned[current] <- res$pruned.labels[current]
        }
    }

    output <- DataFrame(labels=chosen.label, row.names=rownames(results[[1]]))

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
