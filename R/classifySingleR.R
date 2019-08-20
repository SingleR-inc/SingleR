#' Classify cells with SingleR
#' 
#' Assign labels to each cell in a test dataset, using a pre-trained classifier combined with an iterative fine-tuning approach.
#' 
#' @param test A numeric matrix of single-cell expression values where rows are genes and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param trained A \linkS4class{List} containing the output of the \code{\link{trainSingleR}} function.
#' @param quantile A numeric scalar specifying the quantile of the correlation distribution to use to compute the score for each label.
#' @param fine.tune A logical scalar indicating whether fine-tuning should be performed. 
#' @param tune.thresh A numeric scalar specifying the maximum difference from the maximum correlation to use in fine-tuning.
#' @param sd.thresh A numeric scalar specifying the threshold on the standard deviation, for use in gene selection during fine-tuning.
#' This is only used if \code{genes="sd"} when constructing \code{trained} and defaults to the value used in \code{\link{trainSingleR}}.
#' @param assay.type Integer scalar or string specifying the matrix of expression values to use if \code{test} is a \linkS4class{SummarizedExperiment}.
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifyign the parallelization scheme to use.
#' 
#' @return A \linkS4class{DataFrame} where each row corresponds to a cell in \code{test}.
#' If \code{fine.tune=FALSE}, fields are:
#' \itemize{
#' \item \code{scores}, a matrix of correlations at the specified \code{quantile} for each label (column) in each cell (row).
#' \item \code{labels}, the predicted label based on the maximum entry in \code{scores}.
#' }
#'
#' If \code{fine.tune=TRUE}, fields are:
#' \itemize{
#' \item \code{scores}, a matrix of correlations as above.
#' \item \code{first.labels}, the predicted label based on the maximum entry in \code{scores}.
#' \item \code{labels}, the predicted label after fine-tuning.
#' }
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
#'
#' The default \code{assay.type} is set to \code{"logcounts"} simply for consistency with \code{\link{trainSingleR}}.
#' In practice, the raw counts (for UMI data) or the transcript counts (for read count data) can also be used without normalization and log-transformation.
#' Any monotonic transformation will have no effect the calculation of the correlation values other than for some minor differences due to numerical precision.
#'
#' @examples
#' ###########################################
#' ## Mocking up some example training data ##
#' ###########################################
#'
#' Ngroups <- 5
#' Ngenes <- 1000
#' means <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
#' means[1:900,] <- 0
#' colnames(means) <- LETTERS[1:5]
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' sce <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' 
#' rownames(sce) <- sprintf("GENE_%s", seq_len(nrow(sce)))
#' sce <- scater::logNormCounts(sce)
#' 
#' trained <- trainSingleR(sce, sce$label)
#'
#' ##################################################
#' ## Mocking up some test data for classification ##
#' ##################################################
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' test <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' 
#' rownames(test) <- sprintf("GENE_%s", seq_len(nrow(test)))
#' test <- scater::logNormCounts(test)
#' 
#' pred <- classifySingleR(test, trained)
#' table(predicted=pred$labels, truth=g)
#' 
#' @export
#' @importFrom BiocNeighbors KmknnParam bndistance 
#' @importFrom S4Vectors List DataFrame
#' @importFrom SummarizedExperiment colData<- colData 
#' @importFrom BiocParallel SerialParam bpstart bpisup bpstop bplapply
classifySingleR <- function(test, trained, quantile=0.8, 
    fine.tune=TRUE, tune.thresh=0.05, sd.thresh=NULL,
    assay.type="logcounts", check.missing=TRUE, BPPARAM=SerialParam()) 
{
    test <- .to_clean_matrix(test, assay.type, check.missing, msg="test")

    # Don't globally subset 'x', as fine-tuning requires all genes
    # when search.mode='sd'.
    ref.genes <- trained$common.genes
    if (!all(ref.genes %in% rownames(test))) {
        stop("'rownames(test)' does not contain all genes used in 'trained'")
    }

    # Initial search in rank space.
    ranked <- .scaled_colranks_safe(test[ref.genes,,drop=FALSE])
    all.indices <- trained$nn.indices

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Parallelizing across labels rather than cells, as we often have few cells but many labels.
    scores <- bplapply(all.indices, FUN=.find_nearest_quantile, ranked=ranked, quantile=quantile, BPPARAM=BPPARAM)
    scores <- do.call(cbind, scores)
    rownames(scores) <- rownames(ranked)

    # Fine-tuning with an iterative search in lower dimensions.
    labels <- colnames(scores)[max.col(scores)]
    if (fine.tune) {
        search.mode <- trained$search$mode
        if (search.mode=="de") {
            tuned <- .fine_tune_de(exprs=test, scores=scores, references=trained$original.exprs, 
                quantile=quantile, tune.thresh=tune.thresh, de.info=trained$search$extra,
                BPPARAM=BPPARAM)
        } else if (search.mode=="sd") {
            if (is.null(sd.thresh)) {
                sd.thresh <- trained$search$args$sd.thresh
            }
            tuned <- .fine_tune_sd(exprs=test, scores=scores, references=trained$original.exprs, 
                quantile=quantile, tune.thresh=tune.thresh, median.mat=trained$search$extra,
                sd.thresh=sd.thresh, BPPARAM=BPPARAM)
        } else {
            stop(sprintf("unrecognised search mode '%s' when fine-tuning", search.mode))
        }

        new.labels <- colnames(scores)[tuned[[1]]+1L]
        output <- DataFrame(scores=I(scores), first.labels=labels, 
            tuning.scores=I(DataFrame(first=tuned[[2]], second=tuned[[3]])),
            labels=new.labels)
    } else {
        output <- DataFrame(scores=I(scores), labels=labels)
    }

    rownames(output) <- colnames(test)
    output
}

#' @importFrom BiocNeighbors queryKNN
.find_nearest_quantile <- function(ranked, index, quantile) {
    k <- max(1, round(nrow(index) * (1-quantile)))
    nn.d <- queryKNN(query=ranked, k=k, BNINDEX=index, get.index=FALSE, get.distance=FALSE)
    1 - 2*nn.d^2 # see https://arxiv.org/abs/1208.3145
}
