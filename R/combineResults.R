#' Combine SingleR results
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} into a single \linkS4class{DataFrame}.
#' The label from the results with the highest score for each cell is retained.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} or \code{\link{SingleR}}.
#'
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This is identical to the output of \code{\link{classifySingleR}}, albeit with \code{scores} combined and data from the results DataFrame
#' with the highest score carried through for each cell. 
#' The original results are available in the \code{orig.results} field.
#' 
#' @details
#' Given the importance of closely-related reference profiles for proper labeling, use of multiple reference sets can be helpful.
#' The results provided as input for this function should be generated from training sets that \strong{use a common set of genes}.
#'
#' @section Method rationale:
#' There are three obvious options for combining reference datasets or classification results stemming from disparate references:
#'
#' \strong{Option 1} would be to combine the reference datasets into a single matrix and treat each label as though it's specific to the reference
#' from which it originated (e.g. Ref1-Bcell vs Ref2-Bcell), which is easily accomplished by \code{paste}ing the reference name onto the
#' corresponding set of labels. 
#' This option could be useful if the difference between the reference sets were important, 
#' and it also avoids the need for label harmonization between references.
#'
#' However, this method will be prone to enrichment for genes responsible for uninteresting batch effects between the references. 
#' This method would likely lead to a loss of precision and additional noise and risks the potential for technical variation to drive classification.
#'
#' \strong{Option 2} would also include combining the reference datasets into a single matrix but would utilize label harmonization so that
#' the same cell type is given the same label across references. 
#' This would allow feature selection methods to identify robust sets of label-specific markers that are more likely to generalize to other datasets. 
#' It would also simplify interpretation, as there is no need to worry about the reference from which the labels came.
#'
#' The main obstacle to this method is the diffculty and annoyance of harmonization. 
#' Putting aside trivial differences in naming schemes (e.g. "B cell" vs "B"), additional challenges like differences in label resolution 
#' across references (e.g., how to harmonize "B cell" to another reference that splits to "naive B cell" and "mature B cell"), 
#' different sorting strategies, or subtle biological difference that require domain expertise.
#'
#' \strong{Option 3} (the method that this function uses) would be to perform classification separately within each reference, then collate
#' the results to choose the label with the highest score across references. 
#' This avoids the potential for reference-specific markers (like option 1) and the need for explicit harmonization (like option 2).
#' This option leaves a mixture of labels in the final results that is up to the user to resolve.
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
#' \code{\link{matchReferences}}, to harmonize labels between reference datasets.
#' \code{\link{SingleR}}, for generating predictions.
#'
#' @importFrom S4Vectors DataFrame metadata metadata<-
#'
#' @export
combineResults <- function(results) {
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

    output <- DataFrame(scores = I(all.scores), labels = chosen.label)
    metadata(output)$common.genes <- metadata(results[[1]])$common.genes

    if (has.de) {
        metadata(output)$de.genes <- do.call(c, collected.de)
        names(metadata(output)$de.genes) <- de.names
    }

    if (has.first) {
        output$first.labels <- chosen.first
    }

    if (has.pruned) {
        output$pruned.labels <- chosen.pruned
    }

    output$orig.results <- do.call(DataFrame, lapply(results, I))

    output
}
