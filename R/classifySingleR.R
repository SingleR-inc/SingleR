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
#' ########################
#' ## Doing the training ##
#' ########################
#'
#' trained <- trainSingleR(sce, sce$label)
#'
#' sce <- classifySingleR(sce, trained)
#' table(predicted=sce$SingleR$labels, truth=g)
#' 
#' @export
#' @importFrom scran scaledColRanks
#' @importFrom BiocNeighbors KmknnParam bndistance queryKNN
#' @importFrom S4Vectors List DataFrame
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData<- colData assay
#' @importFrom BiocParallel SerialParam
classifySingleR <- function(x, trained, quant.thresh=0.2, fine.tune=TRUE, tune.thresh=0.05, assay.type=1, 
    BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    if (is.null(rownames(x))) {
        stop("'x' must have row names")
    }
    if (was.sce <- is(x, "SingleCellExperiment")) {
        original <- x
        x <- assay(x, i=assay.type)
    }

    ref.genes <- trained$common.genes
    if (all(ref.genes %in% rownames(x))) {
        x <- x[ref.genes,,drop=FALSE]
    } else {
        stop("'rownames(x)' does not contain all genes used in 'trained'")
    }

    # Initial search in rank space.
    ranked <- scaledColRanks(x, transposed=TRUE)
    all.indices <- trained$nn.indices
    scores <- matrix(0, nrow(ranked), length(all.indices), dimnames=list(rownames(ranked), names(all.indices)))
    for (u in names(all.indices)) {
        curdex <- all.indices[[u]]
        k <- max(1, round(nrow(curdex) * quant.thresh))
        nn.out <- queryKNN(query=ranked, k=k, BNINDEX=curdex, get.index=FALSE)
        scores[,u] <- 1 - 2*nn.out$distance[,ncol(nn.out$distance)]^2 # see https://arxiv.org/abs/1208.3145
    }

    # Fine-tuning with an iterative search in lower dimensions.
    labels <- colnames(scores)[max.col(scores)]
    if (fine.tune) {
        new.labels <- .fine_tune_data(exprs=x, scores=scores, genes=trained$search.mode, 
            extras=trained$extra, references=trained$original.exprs, 
            quant.thresh=quant.thresh, tune.thresh=tune.thresh, 
            sd.thresh=trained$sd.thresh, BPPARAM=BPPARAM)
        output <- DataFrame(scores=I(scores), labels=new.labels, first.labels=labels)
    } else {
        output <- DataFrame(scores=I(scores), labels=labels)
    }

    if (was.sce) {
        colData(original)$SingleR <- output
        return(original)
    } else {
        return(output)
    }
}
