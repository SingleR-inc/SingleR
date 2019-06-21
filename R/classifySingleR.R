#' Classify cells with SingleR
#'
#' Predict cell type annotations using a pre-trained SingleR classifier.
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
#' test <- classifySingleR(test, trained)
#' table(predicted=test$SingleR$labels, truth=g)
#' 
#' @export
#' @importFrom BiocNeighbors KmknnParam bndistance queryKNN
#' @importFrom S4Vectors List DataFrame
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- colData assay
#' @importFrom BiocParallel SerialParam
classifySingleR <- function(x, trained, quant.thresh=0.2, 
    fine.tune=TRUE, tune.thresh=0.05, sd.thresh=NULL,
    assay.type=1, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    if (is.null(rownames(x))) {
        stop("'x' must have row names")
    }
    if (was.sce <- is(x, "SingleCellExperiment")) {
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
    sr.out <- .scaled_colranks_safe(x[ref.genes,,drop=FALSE])
    ranked <- sr.out$mat[!sr.out$failed,,drop=FALSE]
    all.indices <- trained$nn.indices

    scores <- matrix(0, nrow(ranked), length(all.indices), dimnames=list(rownames(ranked), names(all.indices)))
    for (u in names(all.indices)) {
        curdex <- all.indices[[u]]
        k <- max(1, round(nrow(curdex) * quant.thresh))
        nn.out <- queryKNN(query=ranked, k=k, BNINDEX=curdex, get.index=FALSE, BPPARAM=BPPARAM)
        scores[,u] <- 1 - 2*nn.out$distance[,ncol(nn.out$distance)]^2 # see https://arxiv.org/abs/1208.3145
    }

    # Fine-tuning with an iterative search in lower dimensions.
    labels <- colnames(scores)[max.col(scores)]
    if (fine.tune) {
        to.use <- x[,!sr.out$failed,drop=FALSE]
        search.mode <- trained$search$mode

        if (search.mode=="de") {
            new.labels <- .fine_tune_de(exprs=to.use, scores=scores, references=trained$original.exprs, 
                quant.thresh=quant.thresh, tune.thresh=tune.thresh, de.info=trained$search$extra,
                BPPARAM=BPPARAM)
        } else if (search.mode=="sd") {
            if (is.null(sd.thresh)) {
                sd.thresh <- trained$search$args$sd.thresh
            }
            new.labels <- .fine_tune_sd(exprs=to.use, scores=scores, references=trained$original.exprs, 
                quant.thresh=quant.thresh, tune.thresh=tune.thresh, median.mat=trained$search$extra,
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

    if (was.sce) {
        colData(original)$SingleR <- output
        return(original)
    } else {
        return(output)
    }
}
