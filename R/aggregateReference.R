#' Aggregate reference samples
#'
#' Aggregate reference samples for a given label by averaging their count profiles.
#' This can be done with varying degrees of resolution to preserve the within-label heterogeneity.
#'
#' @param ref A numeric matrix of reference expression values, usually containing log-expression values.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param labels A character vector or factor of known labels for all cells in \code{ref}.
#' @param ncenters Integer scalar specifying the maximum number of aggregated profiles to produce for each label.
#' @param power Numeric scalar between 0 and 1 indicating how much aggregation should be performed, see Details.
#' @param rank Integer scalar specfiying the number of principal components to use during clustering.
#' @param assay.type An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object.
#' @param ntop Integer scalar specifying the number of highly variable genes to use for the PCA step.
#' @param subset.row Integer, character or logical vector indicating the rows of \code{ref} to use for k-means clustering. 
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how parallelization should be performed.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object indicating which SVD algorithm should be used in \code{\link{runPCA}}.
#' 
#' @details
#' With single-cell reference datasets, it is often useful to aggregate individual cells into pseudo-bulk samples to serve as a reference.
#' This improves speed in downstream assignment with \code{\link{classifySingleR}} or \code{\link{SingleR}}.
#' The most obvious aggregation is to simply average all counts for all cells in a label to obtain a single pseudo-bulk profile.
#' However, this discards information about the within-label heterogeneity (e.g., the \dQuote{shape} and spread of the population in expression space)
#' that may be informative for assignment, especially for closely related labels.
#'
#' The default approach in this function is to create a series of pseudo-bulk samples to represent each label.
#' This is achieved by performing vector quantization via k-means clustering on all cells in a particular label.
#' Cells in each cluster are subsequently averaged to create one pseudo-bulk sample that serves as a representative for that location in the expression space.
#' This reduces the number of separate observations (for speed) while preserving some level of population heterogeneity (for fidelity).
#' 
#' The number of pseudo-bulk samples per label is controlled by \code{ncenters}.
#' By default, we set the number of clusters to \code{X^power} where \code{X} is the number of cells for that label.
#' This ensures that labels with more cells have more resolved representatives.
#' If \code{ncenters} is greater than the number of samples for a label and/or \code{power=1}, no aggregation is performed.
#' Setting \code{power=0} will aggregate all cells of a label into a single pseudo-bulk profile.
#'
#' In practice, k-means clustering is actually performed on the first \code{rank} principal components as computed using \code{\link{runPCA}}.
#' The use of PCs compacts the data for more efficient operation of \code{\link{kmeans}};
#' it also removes some of the high-dimensional noise to highlight major factors of within-label heterogenity.
#' Note that the PCs are only used for clustering and the full expression profiles are still used for the final averaging.
#' Users can disable the PCA step by setting \code{rank=Inf}.
#'
#' By default, we speed things up by only using the top \code{ntop} genes with the largest variances in the PCA.
#' More subsetting of the matrix prior to the PCA can be achieved by setting \code{subset.row} to an appropriate indexing vector.
#' This option may be useful for clustering based on known genes of interest but retaining all genes in the aggregated results.
#' (If both options are set, subsetting by \code{subset.row} is done first, and then the top \code{ntop} genes are selected.)
#' In both cases, though, the aggregation is performed on the full expression profiles.
#'
#' We use the average rather than the sum in order to be compatible with \code{\link{trainSingleR}}'s internal marker detection.
#' Moreover, unlike counts, the sum of transformed and normalized expression values generally has little meaning.
#' We do not use the median to avoid consistently obtaining zeros for lowly expressed genes.
#' 
#' @return
#' A \linkS4class{SummarizedExperiment} object with a \code{"logcounts"} assay containing a matrix of aggregated expression values,
#' and a \code{label} column metadata field specifying the label corresponding to each column.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' 
#' # Making up some labels for demonstration purposes:
#' labels <- sample(LETTERS, ncol(sce), replace=TRUE)
#'
#' # Aggregation at different resolutions:
#' (aggr <- aggregateReference(sce, labels, power=0.5))
#'
#' (aggr <- aggregateReference(sce, labels, power=0))
#'
#' # No aggregation:
#' (aggr <- aggregateReference(sce, labels, power=1))
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam bpnworkers bpmapply
#' @importFrom BiocSingular bsparam
aggregateReference <- function(ref, labels, ncenters=NULL, power=0.5, ntop=1000, assay.type="logcounts", 
    rank=20, subset.row=NULL, check.missing=TRUE, BPPARAM=SerialParam(), BSPARAM=bsparam())
{
    by.label <- split(seq_along(labels), labels)
    if (is(ref, "SummarizedExperiment")) {
        ref <- assay(ref, i=assay.type)
    }

    FRAGMENTER <- function(x, i, j) {
        x <- x[,j,drop=FALSE] # this is usually smaller, so do this first.
        if (!is.null(i)) {
            x <- x[i,,drop=FALSE]
        }
        x
    }

    all.seeds <- .define_seeds(length(by.label))

    if (is.null(BPPARAM) || is(BPPARAM, "SerialParam") || is(BPPARAM, "MulticoreParam")) {
        output.vals <- bpmapply(chosen=by.label, .seed=all.seeds,
            FUN=function(chosen, x, ...) .aggregate_internal(FRAGMENTER(x, i=subset.row, j=chosen), ...), 
            MoreArgs=list(x=ref, ncenters=ncenters, power=power, rank=rank, check.missing=check.missing, ntop=ntop, BSPARAM=BSPARAM),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    } else {
        by.mat <- lapply(by.label, FRAGMENTER, x=ref, i=subset.row)
        output.vals <- bpmapply(by.mat, .seed=all.seeds, FUN=.aggregate_internal,
            MoreArgs=list(ncenters=ncenters, power=power, rank=rank, check.missing=check.missing, ntop=ntop, BSPARAM=BSPARAM),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    }

    if (length(output.vals)==0L) {
        output.vals[[1]] <- matrix(0, nrow(ref), 0, dimnames=list(rownames(ref), NULL))
    }

    first <- labels[vapply(by.label, function(i) i[1], 0L)]
    num <- vapply(output.vals, ncol, 0L)
    output.labels <- rep(first, num)
    
    output <- SummarizedExperiment(list(logcounts=do.call(cbind, output.vals)), colData=DataFrame(label=output.labels))
    colnames(output) <- sprintf("%s.%s", output.labels, sequence(num))
    output
}

#' @importFrom stats kmeans
#' @importFrom utils head
#' @importFrom Matrix t rowMeans
#' @importFrom BiocSingular runPCA
#' @importFrom DelayedArray sweep colsum DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom beachmat realizeFileBackedMatrix
.aggregate_internal <- function(current, ncenters, power, rank, ntop, check.missing, BSPARAM, .seed) {
    oldseed <- .get_seed()
    on.exit(.set_seed(oldseed))

    old <- RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old[1]), add=TRUE, after=FALSE)
    assign(".Random.seed", .seed, envir = .GlobalEnv)

    old.bp <- getAutoBPPARAM()
    setAutoBPPARAM(SerialParam())
    on.exit(setAutoBPPARAM(old.bp), add=TRUE, after=FALSE)

    cur.ncenters <- ncenters
    if (is.null(cur.ncenters)) {
        cur.ncenters <- floor(ncol(current)^power)
    }

    if (cur.ncenters <= 1) {
        val <- matrix(rowMeans(current), dimnames=list(rownames(current), NULL))
    } else if (cur.ncenters >= ncol(current)) {
        val <- current
    } else {
        to.use <- current 

        if (ntop <= nrow(current)) {
            o <- order(rowVars(to.use), decreasing=TRUE)
            to.use <- to.use[head(o, ntop),,drop=FALSE]
        }

        to.use <- t(to.use)
        to.use <- realizeFileBackedMatrix(to.use)
        
        # Identifying the top PCs to avoid realizing the entire matrix.
        if (rank <= min(dim(to.use))-1L) {
            pcs <- runPCA(to.use, rank=rank, get.rotation=FALSE, BSPARAM=BSPARAM)$x
        } else {
            pcs <- as.matrix(to.use)
        }

        clustered <- kmeans(pcs, centers=cur.ncenters)
        val <- colsum(DelayedArray(current), clustered$cluster)
        tab <- table(clustered$cluster)[colnames(val)]
        val <- sweep(val, 2, tab, "/")
    }

    val
}

#' @importFrom stats runif
#' @importFrom parallel nextRNGStream
.define_seeds <- function(n) {
    if (!n) {
        return(list())
    }

    runif(1) # bumping the RNG so that repeated calls to this function generate different results.

    oldseed <- .get_seed()
    on.exit(.set_seed(oldseed))

    old <- RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old[1]), add=TRUE, after=FALSE)

    seeds <- vector("list", n)
    seeds[[1L]] <- .Random.seed
    for (i in seq_len(n - 1L)) {
        seeds[[i + 1L]] <- nextRNGStream(seeds[[i]])
    }

    seeds
}

.get_seed <- function( ) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
}

.set_seed <- function(oldseed) {
    if (!is.null(oldseed)) {
        assign(".Random.seed", oldseed, envir = .GlobalEnv)
    } else {
        rm(.Random.seed, envir = .GlobalEnv)
    }
}
