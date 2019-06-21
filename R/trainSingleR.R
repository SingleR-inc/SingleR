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
#' @importFrom BiocNeighbors KmknnParam bndistance buildIndex KmknnParam
#' @importFrom S4Vectors List
#' @importFrom SummarizedExperiment assay
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
        extra <- genes
        common <- unique(unlist(extra))
        genes <- "de"
    } else if (is.character(genes)) {
        genes <- match.arg(genes, c("de", "sd", "all"))
        if (genes=="de") {
            extra <- .get_genes_by_de(x, labels)
            common <- unique(unlist(extra))
        } else if (genes=="sd") {
            sd.out <- .get_genes_by_sd(x, labels, sd.thresh=sd.thresh)
            common <- sd.out$genes
            extra <- sd.out$mat
        } else {
            common <- rownames(x)
            extra <- NULL
        }
    }

    if (bndistance(BNPARAM)!="Euclidean") {
        stop("'bndistance(BNPARAM)' must be 'Euclidean'") # for distances to be convertible to Spearman rank correlations.
    }
    indices <- original <- List()
    for (u in unique(labels)) {
        current <- x[,labels==u,drop=FALSE] # don't subset by 'common' here, as this loses genes for fine-tuning when genes='sd'.
        original[[u]] <- current
        sr.out <- .scaled_colranks_safe(current[common,,drop=FALSE])
        indices[[u]] <- buildIndex(sr.out$mat[!sr.out$failed,,drop=FALSE], BNPARAM=BNPARAM)
    }

    List(
        common.genes=as.character(common),
        original.exprs=original,
        search.mode=genes,
        nn.indices=indices,
        extra=extra,
        BNPARAM=BNPARAM
    )
}

#' @importFrom utils head
.get_genes_by_de <- function(x, labels) {
    ulabels <- unique(labels)
    n <- round(500*(2/3)^log2(ncol(x)))
    mat <- .median_by_label(x, labels)

    collected <- list()
    for (i in ulabels) {
        subcollected <- list()
        for (j in ulabels) {
            s <- sort(mat[,i] - mat[,j], decreasing=TRUE)
            s <- s[s>0]
            subcollected[[j]] <- as.character(head(names(s), n))
        }
        collected[[i]] <- subcollected
    }
    collected
}

#' @importFrom DelayedMatrixStats rowSds
.get_genes_by_sd <- function(x, labels, sd.thresh=1) {
    mat <- .median_by_label(x, labels)
    sd <- rowSds(mat)
    list(mat=mat, genes=as.character(rownames(mat)[sd > sd.thresh]))
}

#' @importFrom DelayedMatrixStats rowMedians
.median_by_label <- function(mat, labels) {
    ulabels <- unique(labels)
    output <- matrix(0, nrow(mat), length(ulabels))
    rownames(output) <- rownames(mat)
    colnames(output) <- ulabels

    for (u in ulabels) {
        output[,u] <- rowMedians(mat, cols=u==labels)
    }
    output
}

#' @importFrom DelayedMatrixStats colRanks rowSds
.scaled_colranks_safe <- function(x) {
    out <- colRanks(x, ties.method="average")
    center <- (nrow(x) + 1)/2
    sum.sq <- rowVars(out, center=center) * (nrow(x)-1)
    out <- (out - center)/(sqrt(sum.sq) * 2)
    list(mat=out, failed=(sum.sq < 1e-8))
}
