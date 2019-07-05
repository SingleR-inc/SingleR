#' Classify cells with SingleR
#' 
#' @param x A numeric matrix of single-cell expression values (usually log-transformed or otherwise variance-stabilized),
#' where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param trained A \linkS4class{List} containing the output of the \code{\link{trainSingleR}} function.
#' @param quantile A numeric scalar specifying the quantile of the correlation distribution to use to compute the score for each label.
#' @param fine.tune A logical scalar indicating whether fine-tuning should be performed. 
#' @param tune.thresh A numeric scalar specifying the maximum difference from the maximum correlation to use in fine-tuning.
#' @param sd.thresh A numeric scalar specifying the threshold on the standard deviation, for use in gene selection during fine-tuning.
#' This is only used if \code{genes="sd"} when constructing \code{trained} and defaults to the value used in \code{\link{trainSingleR}}.
#' @param assay.type Integer scalar or string specifying the matrix of expression values to use if \code{x} is a \linkS4class{SingleCellExperiment}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifyign the parallelization scheme to use.
#' 
#' @return A \linkS4class{DataFrame} where each row corresponds to a cell in \code{x}.
#' Fields are:
#' \itemize{
#' \item \code{scores}, a matrix of correlations at the specified \code{quantile} for each label (column) in each cell (row).
#' \item \code{labels}, the predicted label based on the maximum entry in \code{scores}.
#' }
#' 
#' @author Aaron Lun, based on the original \code{SingleR} code by Dvir Aran.
#'
#' @details
#' Consider each cell in the test set \code{x} and each label in the training set.
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
#' sce <- SingleCellExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' rownames(sce) <- sprintf("GENE_%s", seq_len(nrow(sce)))
#' 
#' trained <- trainSingleR(sce, sce$label)
#'
#' ##################################################
#' ## Mocking up some test data for classification ##
#' ##################################################
#'
#' N <- 100
#' g <- sample(LETTERS[1:5], N, replace=TRUE)
#' test <- SingleCellExperiment(
#'     list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
#'     colData=DataFrame(label=g)
#' )
#' rownames(test) <- sprintf("GENE_%s", seq_len(nrow(test)))
#' 
#' pred <- classifySingleR(test, trained)
#' table(predicted=pred$labels, truth=g)
#' 
#' @export
#' @importFrom BiocNeighbors KmknnParam bndistance queryKNN
#' @importFrom S4Vectors List DataFrame
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- colData assay
#' @importFrom BiocParallel SerialParam
classifySingleR <- function(x, trained, quantile=0.8, 
    fine.tune=TRUE, tune.thresh=0.05, sd.thresh=NULL,
    assay.type=1, BPPARAM=SerialParam()) 
{
    if (is.null(rownames(x))) {
        stop("'x' must have row names")
    }
    if (is(x, "SingleCellExperiment")) {
        original <- x
        x <- assay(x, i=assay.type)
    }

    # Don't globally subset 'x', as fine-tuning requires all genes
    # when search.mode='sd'.
    ref.genes <- trained$common.genes
    if (!all(ref.genes %in% rownames(x))) {
        stop("'rownames(x)' does not contain all genes used in 'trained'")
    }

    # Initial search in rank space.
    ranked <- .scaled_colranks_safe(x[ref.genes,,drop=FALSE])
    all.indices <- trained$nn.indices

    scores <- matrix(0, nrow(ranked), length(all.indices), dimnames=list(rownames(ranked), names(all.indices)))
    for (u in names(all.indices)) {
        curdex <- all.indices[[u]]
        k <- max(1, round(nrow(curdex) * (1-quantile)))
        nn.out <- queryKNN(query=ranked, k=k, BNINDEX=curdex, get.index=FALSE, BPPARAM=BPPARAM)
        scores[,u] <- 1 - 2*nn.out$distance[,ncol(nn.out$distance)]^2 # see https://arxiv.org/abs/1208.3145
    }

    # Fine-tuning with an iterative search in lower dimensions.
    labels <- colnames(scores)[max.col(scores)]
    if (fine.tune) {
        search.mode <- trained$search$mode
        if (search.mode=="de") {
            new.labels <- .fine_tune_de(exprs=x, scores=scores, references=trained$original.exprs, 
                quantile=quantile, tune.thresh=tune.thresh, de.info=trained$search$extra,
                BPPARAM=BPPARAM)
        } else if (search.mode=="sd") {
            if (is.null(sd.thresh)) {
                sd.thresh <- trained$search$args$sd.thresh
            }
            new.labels <- .fine_tune_sd(exprs=x, scores=scores, references=trained$original.exprs, 
                quantile=quantile, tune.thresh=tune.thresh, median.mat=trained$search$extra,
                sd.thresh=sd.thresh, BPPARAM=BPPARAM)
        } else {
            stop(sprintf("unrecognised search mode '%s' when fine-tuning", search.mode))
        }

        output <- DataFrame(scores=I(scores), labels=new.labels, first.labels=labels)
    } else {
        output <- DataFrame(scores=I(scores), labels=labels)
    }

    # Restoring all failed cells with "NA" values.
    expander <- rep(NA_integer_, ncol(x))
    expander[!sr.out$failed] <- seq_len(nrow(output))
    output <- output[expander,]
    output
}
