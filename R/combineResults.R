#' Combining results from different references
#'
#' It is often desirable to combine information from separate references,
#' thus improving the quality and breadth of the cell type annotation.
#' However, it is not trivial due to the presence of batch effects across references 
#' (from differences in technology, experimental protocol or the biological system)
#' as well as differences in the annotation vocabulary between investigators.
#' This page describes some of the considerations with choosing a strategy
#' to combine information from multiple reference datasets.
#'
#' @section Option 1 - using reference-specific labels:
#' This option nests each label within each reference data, (e.g., \dQuote{Ref1-Bcell} vs \dQuote{Ref2-Bcell}).
#' It is most applicable if there are relevant biological differences between the references,
#' e.g., one reference is concerned with healthy tissue while the other reference considers diseased tissue.
#'
#' In practical terms, this option is easily implemented by just \code{cbind}ing the expression matrices together 
#' and \code{paste}ing the reference name onto the corresponding character vector of labels. 
#' There is no need for time-consuming label harmonization between references.
#'
#' However, the fact that we are comparing across references means that the marker set is likely to contain genes responsible for uninteresting batch effects. 
#' This will increase noise during the calculation of the score in each reference, possibly leading to a loss of precision and a greater risk of technical variation dominating the classification results.
#'
#' @section Option 2 - using harmonized labels:
#' This option also involves combining the reference datasets into a single matrix but with harmonization of the labels so that the same cell type is given the same label across references. 
#' This would allow feature selection methods to identify robust sets of label-specific markers that are more likely to generalize to other datasets. 
#' It would also simplify interpretation, as there is no need to worry about the reference from which the labels came.
#'
#' The most obvious problem with this approach is that it assumes that harmonized labels are available.
#' This is not always the case due to differences in naming schemes (e.g. \code{"B cell"} vs \code{"B"}) between references.
#' Another problem is that of differences in label resolution across references (e.g., how to harmonize \code{"B cell"} to another reference that splits to \code{"naive B cell"} and \code{"mature B cell"}).
#'
#' To mitigate this, \pkg{SingleR} datasets (e.g., \code{\link{ImmGenData}}) have all their labels mapped to the Cell Ontology,
#' allowing the use of standard terms to refer to the same cell type across references.
#' Users can then traverse the ontology graph to achieve a consistent label resolution across references.
#'
#' @section Option 3 - comparing scores across the union of markers:
#' This option involves performing classification separately within each reference, then collating the results to choose the label with the highest score across references. 
#' This is a relatively expedient approach that avoids the need for explicit harmonization while also reduces the potential for reference-specific markers.
#'
#' It leaves a mixture of labels in the final results that is up to the user to resolve, though perhaps this may be considered a feature as it smoothly handles differences in resolution between references, e.g., a cell that cannot be resolved as a CD4+ or CD8+ T cell may simply fall back to \code{"T cell"}.
#' It will also be somewhat suboptimal if there are many reference-specific labels, as markers are not identified with the aim of distinguishing a label in one reference from another label in another reference.
#'
#' @author Aaron Lun
#' @name combine-predictions
#' @seealso
#' \code{\link{combineUnifiedResults}} and \code{\link{combineSeparateResults}},
#' for the functions that implement variants of Option 3.
#'
#' \code{\link{matchReferences}}, to harmonize labels between reference datasets.
NULL

#' Combine SingleR results
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} (usually against different references) into a single \linkS4class{DataFrame}.
#' The label from the results with the highest score for each cell is retained.
#' This assumes that each run of \code{\link{classifySingleR}} was performed using the same set of common genes,
#' hence the \code{Common} in the function name.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} when run on each reference separately.
#'
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This has the same fields as the output of \code{\link{classifySingleR}}, where the \code{scores} are combined across all \code{results}.
#' The set of labels for each cell are those from the DataFrame with the largest maximum score.
#' The original results are available in the \code{orig.results} field.
#' 
#' @details
#' Labels are combined across \code{results} based on the highest score in each reference (see comments in \code{?\link{combine-predictions}}).
#' Each result should be generated from training sets that use a common set of genes during classification, i.e., \code{common.genes} should be the same in the \code{trained} argument to each \code{\link{classifySingleR}} call.
#' This is because the scores are not comparable across results if they were generated from different sets of genes.
#'
#' It is unlikely that this method will be called directly by the end-user.
#' Users are advised to use the multi-reference mode of \code{\link{SingleR}}, \code{\link{trainSingleR}} and/or \code{\link{classifySingleR}}, which will take care of the use of a common set of genes before calling this function to combine results across references.
#'
#' If this function must be called manually, users should ensure that \code{common.genes} is the same for all calls used to generate \code{results}.
#' This is most easily achieved by calling \code{\link{trainSingleR}} on each reference; replacing each \code{common.genes} with the union of all \code{common.genes}; and then calling \code{\link{classifySingleR}} on the test with the modified training objects.
#' The resulting DataFrames can then be passed as \code{results} above.
#'
#' @author Jared Andrews
#'
#' @examples
#' ##############################
#' ## Mocking up training data ##
#' ##############################
#'
#' Ngroups <- 5
#' Ngenes <- 1000
#' means <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
#' means[1:900,] <- 0
#' colnames(means) <- LETTERS[1:5]
#'
#' g <- rep(LETTERS[1:5], each=4)
#' g2 <- rep(LETTERS[6:10], each=4)
#' ref1 <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*length(g), 
#'     lambda=10*2^means[,g]), ncol=length(g))),
#'     colData=DataFrame(label=g)
#' )
#' ref2 <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*length(g2), 
#'     lambda=10*2^means[,g]), ncol=length(g2))),
#'     colData=DataFrame(label=g2)
#' )
#' rownames(ref1) <- sprintf("GENE_%s", seq_len(nrow(ref1)))
#' rownames(ref2) <- sprintf("GENE_%s", seq_len(nrow(ref2)))
#'
#' ref1 <- scater::logNormCounts(ref1)
#' ref2 <- scater::logNormCounts(ref2)
#'
#' ###############################
#' ## Mocking up some test data ##
#' ###############################
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' means <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
#' means[1:900] <- 0
#' colnames(means) <- LETTERS[1:5]
#' test <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' 
#' rownames(test) <- sprintf("GENE_%s", seq_len(nrow(test)))
#' test <- scater::logNormCounts(test)
#' 
#' ###############################
#' ## Performing classification ##
#' ###############################
#' 
#' pred1 <- SingleR(test, ref1, labels=ref1$label)
#' pred2 <- SingleR(test, ref2, labels=ref2$label)
#'
#' pred3 <- SingleR(test, ref1, labels=ref1$label, 
#'     method="cluster", clusters=test$label) 
#' pred4 <- SingleR(test, ref2, labels=ref2$label, 
#'     method="cluster", clusters=test$label) 
#'
#' ###############################
#' ##     Combining results     ##
#' ###############################
#'
#' pred.single <- combineResults(list("pred1" = pred1, "pred2" = pred2))
#' pred.clust <- combineResults(list("pred3" = pred3, "pred4" = pred4))
#'
#' @seealso
#' \code{\link{SingleR}} and \code{\link{classifySingleR}}, for generating predictions to use in \code{results}.
#'
#' @export
#' @importFrom S4Vectors DataFrame metadata metadata<-
combineCommonResults <- function(results) {
    if (length(unique(lapply(results, rownames))) != 1) {
        stop("cell/cluster names are not identical")
    }

    if (length(unique(lapply(results, function (x) sort(metadata(x)$common.genes)))) != 1) {
        # This should be changed to 'stop' before release/after merge with PR #60.
        warning("common genes are not identical")
    }

    has.first <- FALSE
    has.pruned <- FALSE
    has.de <- FALSE

    if (!is.null(results[[1]]$first.labels)) {
        has.first <- TRUE
    }

    if (!is.null(results[[1]]$pruned.labels)) {
        has.pruned <- TRUE
    }

    if (!is.null(metadata(results[[1]])$de.genes)) {
        has.de <- TRUE
        de.names <- sapply(results, function(x) names(metadata(x)$de.genes))
    }

    ncells <- nrow(results[[1]])

    last.best <- rep(-Inf, ncells) 
    chosen.label <- chosen.first <- chosen.pruned <- rep(NA_character_, ncells)
    collected.scores <- list()
    collected.de <- list()

    for (i in seq_along(results)) {
        res <- results[[i]]
        scores <- res$scores
        curbest <- scores[cbind(seq_len(ncells), max.col(scores))]
        collected.scores[[i]] <- scores
        better <- curbest > last.best
         
        chosen.label[better] <- res$labels[better]

        if (has.first) { # either everyone has 'first', or no-one does.
            chosen.first[better] <- res$first.labels[better]
        }

        if (has.pruned) { # same for pruned.
            chosen.pruned[better] <- res$pruned.labels[better]
        }

        if (has.de) {
            collected.de[[i]] <- metadata(res)$de.genes
        }

        last.best <- curbest
    }

    all.scores <- do.call(cbind, collected.scores)

    output <- DataFrame(scores = I(all.scores), row.names=rownames(results[[1]]))
    metadata(output)$common.genes <- metadata(results[[1]])$common.genes

    if (has.de) {
        metadata(output)$de.genes <- do.call(c, collected.de)
        names(metadata(output)$de.genes) <- de.names
    }

    if (has.first) {
        output$first.labels <- chosen.first
    }

    output$labels <- chosen.label

    if (has.pruned) {
        output$pruned.labels <- chosen.pruned
    }

    output$orig.results <- do.call(DataFrame, lapply(results, I))

    output
}
