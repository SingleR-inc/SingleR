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
#' @param prune A logical scalar indicating whether label pruning should be performed.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifyign the parallelization scheme to use.
#' 
#' @return A \linkS4class{DataFrame} where each row corresponds to a cell in \code{test}.
#' If \code{fine.tune=FALSE}, fields are:
#' \itemize{
#' \item \code{scores}, a numeric matrix of correlations at the specified \code{quantile} for each label (column) in each cell (row).
#' \item \code{labels}, a character vector containing the predicted label based on the maximum entry in \code{scores}.
#' \item \code{prune.scores}, a character vector containing the pruned labels where \dQuote{low-quality} labels are replaced with \code{NA}s.
#' Only added if \code{prune=TRUE}.
#' }
#'
#' If \code{fine.tune=TRUE}, fields are:
#' \itemize{
#' \item \code{scores}, a numeric matrix of correlations as above.
#' \item \code{first.labels}, a character vector containing the predicted label \emph{before} fine-tuning.
#' \item \code{tuned.scores}, a DataFrame containing \code{first} and \code{second}.
#' These are numeric vectors containing the best and next-best scores at the final round of fine-tuning for each cell.
#' \item \code{labels}, a character vector containing the predicted label after fine-tuning.
#' \item \code{prune.scores}, a character vector of pruned labels as above.
#' Only added if \code{prune=TRUE}.
#' }
#'
#' In both cases, the \code{\link{metadata}} of the DataFrame contains:
#' \itemize{
#' \item \code{common.genes}, a character vector of genes used to compute the correlations prior to fine-tuning.
#' \item \code{de.genes}, a list of list of genes used to distinguish between each pair of labels.
#' Only returned if \code{genes="de"} when constructing \code{trained}, see \code{?\link{trainSingleR}} for more details. 
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
#' ref <- SummarizedExperiment(
#'     list(counts=matrix(rpois(1000*length(g), 
#'         lambda=10*2^means[,g]), ncol=length(g))),
#'     colData=DataFrame(label=g)
#' )
#' rownames(ref) <- sprintf("GENE_%s", seq_len(nrow(ref)))
#' 
#' ref <- scater::logNormCounts(ref)
#' trained <- trainSingleR(ref, ref$label)
#'
#' ###############################
#' ## Mocking up some test data ##
#' ###############################
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
#' ###############################
#' ## Performing classification ##
#' ###############################
#' 
#' pred <- classifySingleR(test, trained)
#' table(predicted=pred$labels, truth=g)
#' 
#' @seealso
#' \code{\link{trainSingleR}}, to prepare the training set for classification.
#'
#' \code{\link{pruneScores}}, to remove low-quality labels based on the scores.
#'
#' @export
#' @importFrom BiocNeighbors KmknnParam bndistance 
#' @importFrom S4Vectors List DataFrame metadata metadata<-
#' @importFrom SummarizedExperiment colData<- colData 
#' @importFrom BiocParallel SerialParam bpstart bpisup bpstop bplapply
classifySingleR <- function(test, trained, quantile=0.8, fine.tune=TRUE, 
    tune.thresh=0.05, sd.thresh=NULL, prune=TRUE, 
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
    if (length(scores)) { 
        scores <- do.call(cbind, scores)
        labels <- colnames(scores)[max.col(scores)]
    } else {
        scores <- matrix(0, ncol(test), 0) 
        labels <- rep(NA_character_, ncol(test))
    }
    rownames(scores) <- rownames(ranked)

    # Fine-tuning with an iterative search in lower dimensions.
    search.mode <- trained$search$mode
    if (fine.tune) {
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

        if (ncol(scores)) {
            new.labels <- colnames(scores)[tuned[[1]]+1L]
        } else {
            new.labels <- rep(NA_character_, nrow(scores))
        }
        output <- DataFrame(scores=I(scores), first.labels=labels, 
            tuning.scores=I(DataFrame(first=tuned[[2]], second=tuned[[3]])),
            labels=new.labels)
    } else {
        output <- DataFrame(scores=I(scores), labels=labels)
    }

    if (prune) {
        output$pruned.labels <- output$labels
        output$pruned.labels[pruneScores(output)] <- NA_character_
    }

    rownames(output) <- colnames(test)
    metadata(output)$common.genes <- ref.genes
    if(search.mode=="de") {
        metadata(output)$de.genes <- trained$search$extra
    }

    output
}

#' @importFrom BiocNeighbors queryKNN
.find_nearest_quantile <- function(ranked, index, quantile) {
    # We want to find the cells with correlations just before and after 'quantile'.
    # Given correlations 'rho', the quantile value for each observation is:
    # 
    #     (seq_along(rho)-1)/(length(rho)-1)
    #
    # This means that we can just find the floor and ceiling of:
    #
    #     quantile * (length(rho) -1) + 1
    #
    # to obtain the relevant values. Specifically, we consider the floor 
    # as this represents the furthest neighbor in terms of distance. We
    # then convert this to the 'k' in a k-nearest neighbor search.
    denom <- nrow(index) - 1L
    qn <- as.integer(denom * quantile) + 1L
    k <- max(1L, nrow(index) - qn + 1L)

    nn.d <- queryKNN(query=ranked, k=k, last=2, BNINDEX=index, get.index=FALSE, warn.ties=FALSE)$distance
    rho <- 1 - 2*nn.d^2 # see https://arxiv.org/abs/1208.3145

    if (k==1) {
        drop(rho)
    } else {
        # Linear interpolation between the floor/ceiling elements to obtain the
        # quantile. The right weight is that of the higher correlation, and the
        # left weight is that of the lower correlation. 
        rightweight <- quantile - (qn-1)/denom
        furtherweight <- qn/denom - quantile
        (rho[,1] * rightweight + rho[,2] * furtherweight)/(rightweight + furtherweight)
    }
}
