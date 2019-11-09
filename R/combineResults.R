#' Combine SingleR results
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} into a single \linkS4class{DataFrame}.
#' The label with the highest score across predictions for each cell is retained.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} or \code{\link{SingleR}}.
#' 
#' @details
#' Given the importance of closely-related reference profiles for proper labeling, use of multiple reference sets can be helpful.
#' The results provided as input for this function should be generated from training sets that \strong{use a common set of genes}.
#' 
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#'
#' Fields are:
#' \itemize{
#' \item \code{scores}, a numeric matrix of correlations at the specified \code{quantile} for each label (column) in each cell (row), 
#' combined from each object in \code{results}.
#' \item \code{labels}, a character vector containing the predicted label based on the maximum entry in \code{scores} 
#' after combining scores from each results object in \code{results}.
#' }
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
#' trained1 <- trainSingleR(ref1, ref1$label)
#' trained2 <- trainSingleR(ref2, ref2$label)
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
#' pred.single <- combineResults(list(pred1, pred2))
#' pred.clust <- combineResults(list(pred3, pred4))
#'
#' @seealso
#' \code{\link{matchReferences}}, to harmonize labels between reference datasets.
#' \code{\link{SingleR}}, for generating predictions.
#'
#' @importFrom S4Vectors DataFrame
#'
#' @export
combineResults <- function(results) {
    num.features <- sapply(results, nrow)
    if (abs(max(num.features) - min(num.features)) != 0) {
        stop("Results objects contain different numbers of cells or clusters.")
    }

    pred.scores <- list()
    feature.names <- list()

    for (i in seq_along(results)) {
        res <- results[[i]]
        pred.scores[[i]] <- res$scores

        feature.names[[i]] <- rownames(res)
    }

    all.scores <- do.call(cbind, pred.scores)
    best.labels <- colnames(all.scores)[max.col(all.scores)]

    output <- DataFrame(scores=I(all.scores), labels=best.labels)

    if (length(feature.names) > 0) {
        if (length(unique(feature.names)) == 1) {
            rownames(output) <- feature.names[[1]]
        } else {
            warning("Results objects row names do not match, skipping row name assignment.")
        }
    }

    return(output)
}
