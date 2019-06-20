#' Train the SingleR classifier 
#'
#' Train the SingleR classifier on a reference dataset with known labels.
#' 
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
#' trained
#' trained$indices
#' 
#' @export
#' @importFrom scran scaledColRanks
#' @importFrom BiocNeighbors KmknnParam bndistance buildIndex
#' @importFrom S4Vectors List
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
trainSingleR <- function(x, labels, genes="de", sd.thresh=1, assay.type=1, BNPARAM=KmknnParam()) {
    if (is.null(rownames(x))) {
        stop("'x' must have row names")
    }
    if (is(x, "SingleCellExperiment")) {
        x <- assay(x, i=assay.type)
    }

    # Choosing the gene sets of interest. 
    if (is.list(genes)) {
        genes <- "de"
        extra <- genes
        x <- x[unique(unlist(extra)),,drop=FALSE]
    } else if (is.character(genes)) {
        genes <- match.arg(genes, c("de", "sd", "all"))
        if (genes=="de") {
            extra <- getGenesByDE(x, labels)
            x <- x[unique(unlist(extra)),,drop=FALSE]
        } else if (genes=="sd") {
            sd.out <- getGenesBySD(x, labels, sd.thresh=sd.thresh)
            x <- x[sd.out$genes,,drop=FALSE]
            extra <- sd.out$mat
        }
    }

    # Converting to scaled ranks and building an index for each label.
    ranked <- scaledColRanks(x, transposed=TRUE)
    if (bndistance(BNPARAM)!="Euclidean") {
        stop("'bndistance(BNPARAM)' must be 'Euclidean'") # for distances to be convertible to Spearman rank correlations.
    }
    indices <- List()
    for (u in unique(labels)) {
        indices[[u]] <- buildIndex(ranked[labels==u,,drop=FALSE], BNPARAM=BNPARAM)
    }

    List(
        original=x,
        genes=genes,
        indices=indices,
        extra=extra
    )
}

#' @importFrom utils head
getGenesByDE <- function(x, labels) {
    ulabels <- unique(labels)
    n <- round(500*(2/3)^log2(ncol(x)))
    mat <- medianMatrix(x, labels)

    collected <- list()
    for (i in ulabels) {
        subcollected <- list()
        for (j in setdiff(ulabels, i)) {
            s <- sort(mat[,i] - mat[,j], decreasing=TRUE)
            s <- s[s>0]
            subcollected[[j]] <- head(names(s), n)
        }
        collected[[i]] <- subcollected
    }
    collected
}

#' @importFrom DelayedMatrixStats rowSds
getGenesBySE <- function(x, labels, sd.thresh=1) {
    mat <- medianMatrix(x, labels)
    sd <- rowSds(mat)
    list(mat=mat, genes=rownames(mat)[sd > sd.thresh])
}
