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
#' It is also logistically simpler as it allows each reference to be processed separately (more or less, depending on the exact algorithm) for embarrassing parallelization.
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

#' Combine SingleR results with common genes
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} (usually against different references) into a single \linkS4class{DataFrame}.
#' The label from the results with the highest score for each cell is retained.
#' This assumes that each run of \code{\link{classifySingleR}} was performed using a common set of marker genes,
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
#' pred.single <- combineCommonResults(list("pred1" = pred1, "pred2" = pred2))
#' pred.clust <- combineCommonResults(list("pred3" = pred3, "pred4" = pred4))
#'
#' @seealso
#' \code{\link{SingleR}} and \code{\link{classifySingleR}}, for generating predictions to use in \code{results}.
#'
#' \code{\link{combineRecomputedResults}}, for another approach to combining predictions.
#'
#' @export
#' @importFrom S4Vectors DataFrame metadata metadata<-
combineCommonResults <- function(results) {
    if (length(unique(lapply(results, rownames))) != 1) {
        stop("cell/cluster names are not identical")
    }

    all.common <- lapply(results, function (x) sort(metadata(x)$common.genes))
    if (length(unique(all.common)) != 1) {
        # This should be changed to 'stop' before release/after merge with PR #60.
        warning("common genes are not identical")
    }

    ncells <- nrow(results[[1]])
    collected.scores <- collected.best <- vector("list", length(results))
    for (i in seq_along(results)) {
        scores <- results[[i]]$scores
        collected.best[[i]] <- scores[cbind(seq_len(ncells), max.col(scores))]
        collected.scores[[i]] <- scores
    }

    all.scores <- do.call(cbind, collected.scores)
    output <- DataFrame(scores = I(all.scores), row.names=rownames(results[[1]]))
    metadata(output)$common.genes <- all.common[[1]]

    chosen <- max.col(do.call(cbind, collected.best))
    cbind(output, .combine_result_frames(chosen, results))
}

#' @importFrom S4Vectors DataFrame metadata metadata<-
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
    }

    if (has.pruned) {
        output$pruned.labels <- chosen.pruned
    }

    # Collating some DE statistics.
    if (has.de <- !is.null(metadata(results[[1]])$de.genes)) {
        collected.de <- vector("list", length(results))
        for (i in seq_along(results)) {
            collected.de[[i]] <- metadata(results[[i]])$de.genes
        }
        metadata(output)$de.genes <- do.call(c, collected.de)
    }

    output$reference <- chosen
    output$orig.results <- do.call(DataFrame, lapply(results, I))

    output
}

#' Combine SingleR results with recomputation
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} (usually against different references) into a single \linkS4class{DataFrame}.
#' The label from the results with the highest score for each cell is retained.
#' Unlike \code{\link{combineCommonResults}}, this does not assume that each run of \code{\link{classifySingleR}} was performed using the same set of common genes, instead recomputing the scores for comparison across references.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} when run on each reference separately.
#' @inheritParams SingleR
#' @param ref A list or \linkS4class{List} of the same length as \code{results}.
#' Each entry can be a SummarizedExperiment object or numeric matrix containing the reference used to obtain the annotations in the corresponding entry of \code{results}.
#' See the equivalent argument in \code{\link{trainSingleR}} for details.
#' @param labels A list of the same length as \code{ref},
#' where each element should contain a character vector or factor specifying the label for the corresponding entry of \code{ref}.
#'
#' @return A \linkS4class{DataFrame} is returned containing combined annotations for each cell or cluster (row).
#' This has the following fields:
#' \itemize{
#' \item \code{labels}, a character vector containing the best label chosen across all references.
#' \item \code{reference}, an integer vector specifying the reference result from which the best label originates.
#' \item \code{inter.delta.med}, a numeric vector containing the delta value across references.
#' \item \code{intra.delta.med}, a numeric vector containing the delta value within the originating reference.
#' }
#'
#' @details
#' This function implements a variant of Option 3 described in \code{?\link{combine-predictions}}).
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
#' No modification of the output of \code{\link{trainSingleR}} is required as there is no need for a set of common marker genes during the within-reference classifications.
#' However, it is strongly recommended that the universe of genes be the same across all references.
#' Differences in the availability of genes between references may lead to unpredictable combining results.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{SingleR}} and \code{\link{classifySingleR}}, for generating predictions to use in \code{results}.
#'
#' \code{\link{combineCommonResults}}, for another approach to combining predictions.
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
#' combined <- combineRecomputedResults(list(pred1, pred2), test=test,
#'     ref=list(ref1, ref2), labels=list(ref1$label, ref2$label))
#' head(combined)
#'
#' @export
#' @importFrom S4Vectors DataFrame selfmatch metadata
#' @importFrom BiocParallel bpiterate SerialParam
#' @importFrom BiocNeighbors KmknnParam
combineRecomputedResults <- function(results, test, ref, labels, quantile=0.8, 
    assay.type.test="logcounts", assay.type.ref="logcounts", check.missing=TRUE,
    BNPARAM=KmknnParam(), BPPARAM=SerialParam())
{
    all.names <- c(list(colnames(test)), lapply(results, rownames))
    if (length(unique(all.names)) != 1) {
        stop("cell/cluster names are not identical")
    }

    collated <- lapply(results, function(x) x$labels) 
    names(collated) <- make.names(seq_along(collated))
    collated <- DataFrame(collated)
    groups <- selfmatch(collated)
    by.group <- split(seq_along(groups), groups)

    test <- .to_clean_matrix(test, assay.type=assay.type.test, check.missing=check.missing, msg="test")
    ref <- lapply(ref, FUN=.to_clean_matrix, assay.type=assay.type.ref, check.missing=check.missing, msg="ref")

    available <- rownames(test)
    for (r in ref) { 
        available <- intersect(available, rownames(r)) 
    }

    ##############################################
    # Using iterator to avoid having to serialize the entire object to the workers.
    # Rather, only the data for the relevant markers for the current labels is serialized.
    envir <- new.env()
    envir$counter <- 1L
    ITER <- function() {
        if (envir$counter > length(by.group)) {
            return(NULL)
        }
        curgroup <- by.group[[envir$counter]]
        curlabels <- collated[curgroup[1],,drop=FALSE]
        envir$counter <- envir$counter + 1L

        curmarkers <- character(0)
        for (i in seq_along(ref)) {
            M <- metadata(results[[i]])$de.genes
            if (is.null(M)) {
                to.add <- metadata(results[[i]])$common.genes
            } else {
                curlab <- curlabels[[i]]
                to.add <- M[[curlab]]
                to.add <- unlist(to.add, use.names=FALSE)
            }
            curmarkers <- union(curmarkers, to.add)
        }

        curmarkers <- intersect(curmarkers, available)
        all.ref <- vector("list", length(ref))
        for (i in seq_along(ref)) {
            keep <- curlabels[[i]]==labels[[i]]
            all.ref[[i]] <- ref[[i]][curmarkers,keep,drop=FALSE]
        }

        list(
            test=test[curmarkers,curgroup,drop=FALSE],
            ref=do.call(cbind, all.ref),
            label=rep(seq_along(all.ref), vapply(all.ref, ncol, 0L)),
            names=unlist(curlabels)
        )
    }

    combined <- bpiterate(
        ITER=ITER,
        FUN=function(data, BNPARAM, quantile, ...) {
            test <- data$test
            ref <- data$ref
            label <- data$label
            names <- data$names

            # We have no interest looking for marker genes between references,
            # as this would leave us susceptible to batch effects. So, don't
            # bother to detect genes or to fine-tune, we just want the scores
            # and the identity of the label that delivers the maximum score.
            trained <- trainSingleR(ref, label, genes = "all", 
                check.missing = FALSE, BNPARAM = BNPARAM)
            output <- classifySingleR(test, trained, quantile = quantile, fine.tune = FALSE, 
                prune = FALSE, check.missing = FALSE)
   
            output$labels <- as.integer(output$labels)
            metadata(output)$names <- names
            output
        },
        BPPARAM=BPPARAM,
        quantile=quantile,
        BNPARAM=BNPARAM
    )

    ##############################################
    # Organizing scores for output.

    base.scores <- lapply(results, function(x) {
        mat <- x$scores   
        mat[] <- NA_real_
        mat
    })

    for (i in seq_along(combined)) {
        chosen <- by.group[[i]]
        stats <- combined[[i]]
        names <- metadata(stats)$names
        for (r in seq_along(names)) {
            base.scores[[r]][chosen,names[r]] <- stats$scores[,r]
        }
    }

    all.scores <- do.call(cbind, base.scores)
    output <- DataFrame(scores = I(all.scores), row.names=rownames(results[[1]]))

    o <- order(unlist(by.group))
    chosen <- unlist(lapply(combined, function(x) x$labels))[o]
    cbind(output, .combine_result_frames(chosen, results))
}
