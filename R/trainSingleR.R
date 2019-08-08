#' Train the SingleR classifier 
#'
#' Train the SingleR classifier on a reference dataset with known labels.
#' 
#' @param x A numeric matrix of single-cell expression values (usually log-transformed or otherwise variance-stabilized),
#' where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param labels A character vector or factor of known labels for all cells in \code{x}.
#' @param genes A string specifying the feature selection method to be used, see Details.
#' 
#' Alternatively, a list of lists of character vectors containing DE genes between pairs of labels.
#' 
#' Alternatively, a list of character vectors containing marker genes for each label.
#' @param sd.thresh A numeric scalar specifying the minimum threshold on the standard deviation per gene.
#' Only used when \code{genes="sd"}.
#' @param de.n An integer scalar specifying the number of DE genes to use when \code{genes="de"}.
#' Defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
#' @param assay.type An integer scalar or string specifying the assay of \code{x} containing the relevant expression matrix,
#' if \code{x} is a \linkS4class{SingleCellExperiment} object.
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for building nearest neighbor indices.
#'
#' @return A \linkS4class{List} containing:
#' \describe{
#' \item{\code{common.genes}:}{A character vector of all genes that were chosen by the designated feature selection method.}
#' \item{\code{nn.indices}:}{A List of \linkS4class{BiocNeighborIndex} objects containing pre-constructed neighbor search indices.}
#' \item{\code{original.exprs}:}{A List of numeric matrices where each matrix contains all cells for a particular label.}
#' \item{\code{extra}:}{A List of additional information on the feature selection, for use by \code{\link{classifySingleR}}.
#' This includes the selection method in \code{genes} and method-specific structures that can be re-used during classification.}
#' }
#'
#' @details
#' This function uses a training data set to select interesting features and construct nearest neighbor indices in rank space.
#' The resulting objects can be re-used across multiple classification steps with different test data sets via \code{\link{classifySingleR}}.
#' This improves efficiency by avoiding unnecessary repetition of steps during the downstream analysis.
#' 
#' Several options are available for feature selection, usually assuming that \code{x} has been log-transformed or equivalent.
#' \itemize{
#' \item \code{genes="de"} identifies genes that are differentially expressed between labels.
#' This is done by identifying the median expression within each label, and computing differences between medians for each pair of labels.
#' For each label, the top \code{de.n} genes with the largest differences compared to another label are chosen as markers to distinguish the two labels.
#' The set of all features is defined as the union of markers from all pairwise comparisons.
#' \item \code{genes="sd"} identifies genes that are highly variable across labels.
#' This is done by identifying the median expression within each label, and computing the standard deviation in the medians across all labels.
#' The set of all features is defined as those genes with standard deviations above \code{sd.thresh}.
#' \item \code{genes="all"} will not perform any feature selection.
#' }
#'
#' Alternatively, \code{genes} can be a named list of named lists of character vectors:
#' \preformatted{genes <- list(
#'    A = list(A = character(0), B = "GENE_1", C = c("GENE_2", "GENE_3")),
#'    B = list(A = "GENE_100", B = character(0), C = "GENE_200"),
#'    C = list(A = c("GENE_4", "GENE_5"), B = "GENE_5", C = character(0))
#' )
#' }
#' If we consider the entry \code{genes$A$B}, this contains marker genes for label \code{"A"} against label \code{"B"}.
#' That is, these genes are upregulated in \code{"A"} compared to \code{"B"}.
#' The outer list should have one list per label, and each inner list should have one character vector per label.
#' (Obviously, a label cannot have markers against itself, so this is just set to \code{character(0)}.)
#'
#' Alternatively, \code{genes} can be a named list of character vectors:
#' \preformatted{genes <- list(
#'      A = c("GENE_1", "GENE_2", "GENE_3"),
#'      B = c("GENE_100", "GENE_200"),
#'      C = c("GENE_4", "GENE_5")
#' )
#' }
#' The entry \code{genes$A} represent the genes that are upregulated in \code{A} compared to some or all other labels. 
#' This allows the function to handle pre-defined marker lists for specific cell populations.
#' However, it obviously captures less information than marker sets for the pairwise comparisons.
#'
#' If \code{genes} explicitly contains gene identities (as character vectors), \code{x} can be the raw counts or any monotonic transformation thereof.
#' Use of log-transformed counts is only necessary for automated marker detection when \code{genes="de"} or \code{"sd"}.
#'
#' @author Aaron Lun, based on the original \code{SingleR} code by Dvir Aran.
#' 
#' @seealso
#' \code{\link{classifySingleR}}, where the output of this function gets used.
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
#' # Normalizing and log-transforming for automated marker detection.
#' sce <- scater::logNormCounts(sce)
#'
#' trained <- trainSingleR(sce, sce$label)
#' trained
#' trained$nn.indices
#' length(trained$common.genes)
#'
#' # Alternatively, supplying a set of label-specific markers
#' # (in this case, the last 100 genes are DE for each cluster).
#' last.100 <- tail(rownames(sce), 100)
#' markers <- rep(list(last.100), 5)
#' names(markers) <- LETTERS[1:5]
#' 
#' trained <- trainSingleR(sce, sce$label, genes=markers)
#' trained$common.genes
#' 
#' @export
#' @importFrom BiocNeighbors KmknnParam bndistance buildIndex KmknnParam
#' @importFrom S4Vectors List
trainSingleR <- function(x, labels, genes="de", sd.thresh=1, de.n=NULL, 
    assay.type="logcounts", check.missing=TRUE, BNPARAM=KmknnParam()) 
{
    x <- .to_clean_matrix(x, assay.type, check.missing, msg="x")

    # Choosing the gene sets of interest. 
    args <- list()
    if (is.list(genes)) {
        is.char <- vapply(genes, is.character, TRUE)
        if (all(is.char)) {
            genes <- .convert_per_label_set(genes)
        } else if (any(is.char)) {
            stop("'genes' must be a list of character vectors or a list of list of vectors")
        }
        .validate_de_gene_set(genes, labels)
        extra <- genes
        common <- unique(unlist(extra))
        genes <- "de"
    } else if (is.character(genes)) {
        genes <- match.arg(genes, c("de", "sd", "all"))
        if (genes=="de") {
            extra <- .get_genes_by_de(x, labels, de.n=de.n)
            common <- unique(unlist(extra))
        } else if (genes=="sd") {
            sd.out <- .get_genes_by_sd(x, labels, sd.thresh=sd.thresh)
            common <- sd.out$genes
            extra <- sd.out$mat
            args$sd.thresh <- sd.thresh
        } else {
            common <- rownames(x)
            extra <- .median_by_label(x, labels)
            genes <- "sd"
            args$sd.thresh <- sd.thresh
        }
    }

    if (bndistance(BNPARAM)!="Euclidean") {
        stop("'bndistance(BNPARAM)' must be 'Euclidean'") # for distances to be convertible to Spearman rank correlations.
    }
    indices <- original <- List()
    for (u in unique(labels)) {
        # Don't subset by 'common' here, as this loses genes for fine-tuning when genes='sd'.
        current <- x[,labels==u,drop=FALSE] 

        # Coerce to double to make life easier for the C++ code later.
        if (!is.double(current[0,])) {
            current <- current + 0
        }

        original[[u]] <- current
        sr.out <- .scaled_colranks_safe(current[common,,drop=FALSE])
        indices[[u]] <- buildIndex(sr.out, BNPARAM=BNPARAM)
    }

    List(
        common.genes=as.character(common),
        original.exprs=original,
        nn.indices=indices,
        search=List(mode=genes, args=args, extra=extra)
    )
}

#' @importFrom utils head
.get_genes_by_de <- function(x, labels, de.n=NULL) {
    ulabels <- unique(labels)
    mat <- .median_by_label(x, labels)
    if (is.null(de.n)) {
        de.n <- round(500*(2/3)^log2(ncol(mat)))
    }

    collected <- list()
    for (i in ulabels) {
        subcollected <- list()
        for (j in ulabels) {
            s <- sort(mat[,i] - mat[,j], decreasing=TRUE)
            s <- s[s>0]
            subcollected[[j]] <- as.character(head(names(s), de.n))
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

.convert_per_label_set <- function(genes) {
    # Converting into a list-of-lists format so that it plays nice with downstream methods.
    # This is done by saving the markers on the diagonal so that each label's markers are
    # included in the set of genes to use during fine-tuning. Don't exclude the diagonal!
    all.labs <- names(genes)
    for (i in all.labs) {
        empty <- rep(list(character(0)), length(all.labs))
        names(empty) <- all.labs
        empty[[i]] <- genes[[i]]
        genes[[i]] <- empty
    }
    genes
}

.validate_de_gene_set <- function(genes, labels) {
    ulabels <- unique(labels)
    if (!all(ulabels %in% names(genes))) {
        stop("need marker gene information for each label")
    }
    for (u in ulabels) {
        if (!all(ulabels %in% names(genes[[u]]))) {
            stop("need marker genes between each pair of labels")
        }
    }
}

#' @importFrom DelayedMatrixStats rowMedians
#' @importFrom DelayedArray DelayedArray
.median_by_label <- function(mat, labels) {
    ulabels <- unique(labels)
    output <- matrix(0, nrow(mat), length(ulabels))
    rownames(output) <- rownames(mat)
    colnames(output) <- ulabels

    for (u in ulabels) {
        output[,u] <- rowMedians(DelayedArray(mat), cols=u==labels)
    }
    output
}

#' @importFrom DelayedMatrixStats colRanks rowVars
#' @importFrom DelayedArray DelayedArray
.scaled_colranks_safe <- function(x) {
    out <- colRanks(DelayedArray(x), ties.method="average")
    center <- (nrow(x) + 1)/2
    sum.sq <- rowVars(out, center=center) * (nrow(x)-1)
    sum.sq <- pmax(1e-8, sum.sq)
    (out - center)/(sqrt(sum.sq) * 2)
}
