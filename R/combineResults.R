#' Combine SingleR results
#'
#' Combine results from multiple runs of \code{\link{classifySingleR}} into a single \linkS4class{DataFrame}.
#' The label with the highest score across predictions for each cell is retained.
#'
#' @param results A list of \linkS4class{DataFrame} prediction results as returned by \code{\link{classifySingleR}} or \code{\link{SingleR}}.
#' 
#' @details
#' Given the importance of closely-related reference profiles for proper labeling, use of multiple reference sets can be helpful.
#' 
#' @return A \linkS4class{DataFrame} is returned containing the annotation statistics for each cell or cluster (row).
#' This is identical to the output of \code{\link{classifySingleR}}.
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
#' means1 <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
#' means2 <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
#' means1[1:900,] <- 0
#' means2[1:900,] <- 0
#' colnames(means1) <- LETTERS[1:5]
#' colnames(means2) <- LETTERS[6:10]
#'
#' g1 <- rep(LETTERS[1:5], each=4)
#' g2 <- rep(LETTERS[6:10], each=4)
#' ref1 <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*length(g1), 
#'         lambda=10*2^means1[,g1]), ncol=length(g1))),
#'     colData=DataFrame(label=g1)
#' )
#' ref2 <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*length(g2), 
#'         lambda=10*2^means2[,g2]), ncol=length(g2))),
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
#' test <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means1[,g]), ncol=N)),
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
#' ###############################
#' ##     Combining results     ##
#' ###############################
#'
#' pred <- combineResults(list(pred1, pred2))
#'
#' @seealso
#' \code{\link{matchReference}}, to harmonize labels between reference datasets.
#' \code{\link{SingleR}}, for generating predictions.
#'
#'
#' @export
combineResults <- function(results) {
    pred.scores <- list()

    for (p in results) {
        pred.scores[[p]] <- p$scores
    }

    all.scores <- do.call(cbind, pred.scores)
    best.labels <- colnames(all.scores)[max.col(all.scores)]

    output <- DataFrame(scores=I(all.scores), first.labels=best.labels, labels=best.labels)

    output
}

