#' Classify cells with SingleR
#' 
#' Assign labels to each cell in a test dataset, using a pre-trained classifier combined with an iterative fine-tuning approach.
#' 
#' @param test A numeric matrix of single-cell expression values where rows are genes and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param trained A \linkS4class{List} containing the output of the \code{\link{trainSingleR}} function.
#' Alternatively, a List of Lists produced by \code{\link{trainSingleR}} for multiple references.
#' @param quantile A numeric scalar specifying the quantile of the correlation distribution to use to compute the score for each label.
#' @param fine.tune A logical scalar indicating whether fine-tuning should be performed. 
#' @param tune.thresh A numeric scalar specifying the maximum difference from the maximum correlation to use in fine-tuning.
#' @param sd.thresh A numeric scalar specifying the threshold on the standard deviation, for use in gene selection during fine-tuning.
#' This is only used if \code{genes="sd"} when constructing \code{trained} and defaults to the value used in \code{\link{trainSingleR}}.
#' @param assay.type Integer scalar or string specifying the matrix of expression values to use if \code{test} is a \linkS4class{SummarizedExperiment}.
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param prune A logical scalar indicating whether label pruning should be performed.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifyign the parallelization scheme to use.
#' 
#' @return A \linkS4class{DataFrame} where each row corresponds to a cell in \code{test}.
#' In the case of a single reference, this contains:
#' \itemize{
#' \item \code{scores}, a numeric matrix of correlations at the specified \code{quantile} 
#' for each label (column) in each cell (row).
#' This will contain \code{NA}s if multiple references were supplied to \code{\link{trainSingleR}} with \code{recompute=TRUE}.
#' \item \code{first.labels}, a character vector containing the predicted label \emph{before} fine-tuning.
#' Only added if \code{fine.tune=TRUE}.
#' \item \code{tuned.scores}, a DataFrame containing \code{first} and \code{second}.
#' These are numeric vectors containing the best and next-best scores at the final round of fine-tuning for each cell.
#' Only added if \code{fine.tune=TRUE}.
#' \item \code{labels}, a character vector containing the predicted label based on the maximum entry in \code{scores}.
#' \item \code{pruned.labels}, a character vector containing the pruned labels where \dQuote{low-quality}.
#' els are replaced with \code{NA}s.
#' Only added if \code{prune=TRUE}.
#' }
#'
#' The \code{\link{metadata}} of the DataFrame contains:
#' \itemize{
#' \item \code{common.genes}, a character vector of genes used to compute the correlations prior to fine-tuning.
#' \item \code{de.genes}, a list of list of genes used to distinguish between each pair of labels.
#' Only returned if \code{genes="de"} when constructing \code{trained}, see \code{?\link{trainSingleR}} for more details. 
#' }
#'
#' In the case of multiple references, the output of \code{\link{combineCommonResults}} 
#' or \code{\link{combineRecomputedResults}} is returned,
#' depending on whether \code{recompute=TRUE} when constructing \code{trained}.
#' This is a \linkS4class{DataFrame} containing:
#' \itemize{
#' \item \code{scores}, a numeric matrix of scores for each cell (row) across all labels in all references (columns).
#' This will contain \code{NA}s if recomputation is performed.
#' \item \code{labels}, \code{first.labels} (if \code{fine.tune=TRUE}) and \code{pruned.labels} (if \code{prune=TRUE}),
#' containing the consolidated labels of varying flavors as described above.
#' \item \code{orig.results}, a \linkS4class{DataFrame} of DataFrames containing 
#' the results of running \code{\link{classifySingleR}} against each individual reference.
#' Each nested DataFrame has the same format as described above.
#' }
#' See \code{?\link{combineCommonResults}} and \code{?\link{combineRecomputedResults}} for more details.
#' 
#' @author Aaron Lun, based on the original \code{SingleR} code by Dvir Aran.
#'
#' @details
#' Consider each cell in the test set \code{test} and each label in the training set.
#' We compute Spearman's rank correlations between the test cell and all cells in the training set with the given label,
#' based on the expression profiles of the genes selected by \code{trained}.
#' The score is defined as the quantile of the distribution of correlations, as specified by \code{quantile}.
#' (Technically, we avoid explicitly computing all correlations by using a nearest neighbor search, but the resulting score is the same.)
#' After repeating this across all labels, the label with the highest score is used as the prediction for that cell.
#'
#' If \code{fine.tune=TRUE}, an additional fine-tuning step is performed for each cell to improve resolution.
#' We identify all labels with scores that are no more than \code{tune.thresh} below the maximum score.
#' These labels are used to identify a fresh set of marker genes, and the calculation of the score is repeated using only these genes.
#' The aim is to refine the choice of markers and reduce noise when distinguishing between closely related labels.
#' The best and next-best scores are reported in the output for use in diagnostics, e.g., \code{\link{pruneScores}}.
#'
#' The default \code{assay.type} is set to \code{"logcounts"} simply for consistency with \code{\link{trainSingleR}}.
#' In practice, the raw counts (for UMI data) or the transcript counts (for read count data) can also be used without normalization and log-transformation.
#' Any monotonic transformation will have no effect the calculation of the correlation values other than for some minor differences due to numerical precision.
#'
#' If \code{prune=TRUE}, label pruning is performed as described in \code{\link{pruneScores}} with default arguments.
#' This aims to remove low-quality labels that are ambiguous or correspond to misassigned cells.
#' However, the default settings can be somewhat aggressive and discard otherwise useful labels in some cases - see \code{?\link{pruneScores}} for details.
#'
#' If \code{trained} was generated from multiple references, the per-reference statistics are combined into a single DataFrame of results.
#' This is done using \code{\link{combineRecomputedResults}} if \code{recompute=TRUE} in \code{\link{trainSingleR}},
#' otherwise it is done using \code{\link{combineCommonResults}}.
#'
#' @examples
#' # Mocking up data with log-normalized expression values:
#' ref <- .mockRefData()
#' test <- .mockTestData(ref)
#'
#' ref <- scuttle::logNormCounts(ref)
#' test <- scuttle::logNormCounts(test)
#'
#' # Setting up the training:
#' trained <- trainSingleR(ref, label=ref$label)
#'
#' # Performing the classification:
#' pred <- classifySingleR(test, trained)
#' table(predicted=pred$labels, truth=test$label)
#' 
#' @seealso
#' \code{\link{trainSingleR}}, to prepare the training set for classification.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores.
#'
#' \code{\link{combineCommonResults}}, to combine results from multiple references.
#'
#' @export
#' @useDynLib SingleR
classifySingleR <- function(
    test, 
    trained, 
    quantile=0.8, 
    fine.tune=TRUE, 
    tune.thresh=0.05, 
    sd.thresh=NULL, 
    prune=TRUE, 
    assay.type="logcounts", 
    check.missing=TRUE, 
    BPPARAM=SerialParam()) 
{
    test <- .to_clean_matrix(test, assay.type, check.missing, msg="test", BPPARAM=BPPARAM)

    # Unfortunately, we can't test for List, because each trained structure is
    # also a list; so we just check whether the 'built' field exists.
    if (solo <- !is.null(trained$built)) { 
        trained <- list(trained)
    }

    results <- lapply(trained, FUN=.classify_internals, test=test, quantile=quantile, 
        fine.tune=fine.tune, tune.thresh=tune.thresh, prune=prune)

    if (solo) {
        results[[1]]
    } else {
        combineRecomputedResults(results, test=test, trained=trained, 
            check.missing=FALSE, quantile=quantile)
    } 
}

#' @importFrom S4Vectors DataFrame metadata metadata<- I
.classify_internals <- function(test, trained, quantile, fine.tune, tune.thresh=0.05, prune=TRUE) {
    m <- match(trained$markers$unique, rownames(test))
    if (anyNA(m)) {
        stop("'rownames(test)' does not contain all genes used in 'trained'")
    }

    out <- run(test, m - 1L, trained$built, quantile, fine.tune, tune.thresh)
    colnames(out$scores) <- trained$labels$unique
    output <- DataFrame(
        scores = I(out$scores), 
        labels = trained$labels$unique[out$best + 1L],
        gaps = out$delta
    )

    if (prune) {
        output$pruned.labels <- output$labels
        output$pruned.labels[pruneScores(output)] <- NA_character_
    }

    rownames(output) <- colnames(test)
    metadata(output)$common.genes <- trained$markers$unique
    metadata(output)$de.genes <- trained$markers$full

    output
}
