#' Combine SingleR results with recomputation
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} (usually against different references) into a single \link[S4Vectors]{DataFrame}.
#' This involves recomputing the scores so that they are comparable across references.
#'
#' @param results A list of \link[S4Vectors]{DataFrame} prediction results as returned by \code{\link{classifySingleR}} when run on each reference separately.
#' This should have the same names as \code{trained}.
#' @inheritParams SingleR
#' @param check.missing Deprecated and ignored, as any row filtering will cause mismatches with the \code{test.genes=} used in \code{\link{trainSingleR}}.
#' @param trained A list of \link[S4Vectors]{List}s containing the trained outputs of multiple references,
#' equivalent to either (i) the output of \code{\link{trainSingleR}} on multiple references with \code{recompute=TRUE},
#' or (ii) running \code{trainSingleR} on each reference separately and manually making a list of the trained outputs.
#' The list should have the same names as \code{results}.
#' @param warn.lost Logical scalar indicating whether to emit a warning if markers from one reference in \code{trained} are absent in other references.
#' @param quantile Numeric scalar specifying the quantile of the correlation distribution to use for computing the score, see \code{\link{classifySingleR}}.
#' @param fine.tune A logical scalar indicating whether fine-tuning should be performed. 
#' @param tune.thresh A numeric scalar specifying the maximum difference from the maximum correlation to use in fine-tuning.
#' @param allow.lost Deprecated.
#'
#' @return A \link[S4Vectors]{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This mimics the output of \code{\link{classifySingleR}} and contains the following fields:
#' \itemize{
#' \item \code{scores}, a DataFrame of DataFrames containing the \emph{recomputed} scores for the best label in each reference.
#' Each nested DataFrame corresponds to a reference and contains \code{labels} (the best label for that cell in this reference) and \code{scores} (the recomputed score).
#' \item \code{labels}, a character vector containing the per-cell combined label across references.
#' \item \code{reference}, an integer vector specifying the reference from which the combined label was derived.
#' \item \code{delta.next}, a numeric vector containing the difference between the best and next-best score.
#' If \code{fine.tune=TRUE}, this is reported for scores after fine-tuning.
#' \item \code{orig.results}, a DataFrame containing \code{results}.
#' }
#' It may also contain \code{pruned.labels} if these were also present in \code{results}.
#'
#' @details
#' Here, the strategy is to perform classification separately within each reference, 
#' then collate the results to choose the label with the highest score across references.
#' For a given cell in \code{test}, we extract its assigned label from each reference in \code{results}, along with the marker genes associated with that label.
#' We take the union of the markers for the assigned labels across all references.
#' This defines a common feature space in which the score for each reference's assigned label is recomputed using \code{ref};
#' the label from the reference with the top recomputed score is then reported as the combined annotation for that cell.
#' 
#' A key aspect of this approach is that each entry of \code{results} is generated separately for each reference.
#' This avoids problems with unintersting technical or biological differences between references that could otherwise introduce noise by forcing irrelevant genes into the marker list. 
#' Similarly, the common feature space for each cell is defined from the most relevant markers across all references,
#' analogous to one iteration of fine-tuning using only the best labels from each reference.
#' Indeed, if fine-tuning is enabled, the common feature space is iteratively refined from the labels with the highest scores, using the same process described in \code{\link{classifySingleR}}.
#' This allows us to distinguish between closely-related labels from different references.
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
#' @references
#' Lun A, Bunis D, Andrews J (2020).
#' Thoughts on a more scalable algorithm for multiple references.
#' \url{https://github.com/SingleR-inc/SingleR/issues/94}
#'
#' @examples
#' # Making up data.
#' ref <- .mockRefData(nreps=8)
#' partition <- seq_len(ncol(ref)) %% 2
#' ref1 <- ref[,partition == 0]
#' ref2 <- ref[,partition == 1]
#' ref2$label <- tolower(ref2$label) # make lower-case for some more variety.
#'
#' test <- .mockTestData(ref)
#'
#' # Performing classification within each reference.
#' test <- scrapper::normalizeRnaCounts.se(test)
#'
#' ref1 <- scrapper::normalizeRnaCounts.se(ref1)
#' train1 <- trainSingleR(ref1, labels=ref1$label)
#' pred1 <- classifySingleR(test, train1)
#'
#' ref2 <- scrapper::normalizeRnaCounts.se(ref2)
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
#' @aliases
#' combineCommonResults
#' @export
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom beachmat initializeCpp
combineRecomputedResults <- function(
    results, 
    test, 
    trained, 
    quantile=0.8, 
    fine.tune=TRUE, 
    tune.thresh=0.05, 
    assay.type.test="logcounts", 
    check.missing=FALSE,
    warn.lost=TRUE,
    allow.lost=FALSE, 
    num.threads = 1,
    BPPARAM= NULL
) {
    num.threads <- .get_num_threads(num.threads, BPPARAM)
    test <- .to_clean_matrix(test, assay.type=assay.type.test, check.missing=FALSE, msg="test", num.threads=num.threads)

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

    if (!identical(names(results), names(trained))) {
        warning("'results' and 'trained' have different names")
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

    all.inter.test <- all.inter.ref <- vector("list", length(trained))
    test.genes <- rownames(test)
    for (i in seq_along(all.refnames)) {
        inter <- .create_intersection(test.genes, all.refnames[[i]])
        all.inter.test[[i]] <- inter$test - 1L
        all.inter.ref[[i]] <- inter$reference - 1L
    }

    # Applying the integration.
    ibuilt <- train_integrated(
        test_nrow=length(test.genes),
        test_features=all.inter.test,
        references=lapply(trained, function(x) initializeCpp(x$ref, .check.na=FALSE)),
        ref_features=all.inter.ref,
        labels=lapply(trained, function(x) match(x$labels$full, x$labels$unique) - 1L),
        prebuilt=lapply(trained, function(x) rebuildIndex(x)$built),
        nthreads = num.threads
    )

    collated <- vector("list", length(trained))
    for (i in seq_along(collated)) {
        collated[[i]] <- match(results[[i]]$labels, trained[[i]]$labels$unique) - 1L
    }

    parsed <- initializeCpp(test, .check.na=FALSE)
    irun <- classify_integrated(
        test=parsed,
        results=collated,
        integrated_build=ibuilt,
        quantile=quantile,
        use_fine_tune = fine.tune, 
        fine_tune_threshold = tune.thresh, 
        nthreads=num.threads
    ) 

    # Organizing the outputs.
    if (is.null(names(results))) {
        names(results) <- sprintf("ref%i", seq_along(results))
    }

    base.scores <- vector("list", length(results))
    names(base.scores) <- names(results)
    for (i in seq_along(base.scores)) {
        base.scores[[i]] <- DataFrame(labels=results[[i]]$labels, scores=irun$scores[,i])
    }

    all.scores <- DataFrame(lapply(base.scores, I))
    output <- DataFrame(scores = I(all.scores), row.names=rownames(results[[1]]))
    cbind(output, .combine_result_frames(irun$best + 1L, irun$delta, results))
}

#' @importFrom S4Vectors DataFrame
.combine_result_frames <- function(chosen, delta, results) {
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
    output$delta.next <- delta
    output$orig.results <- do.call(DataFrame, lapply(results, I))

    output
}
