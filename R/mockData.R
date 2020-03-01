#' Mock data for examples
#'
#' Make up some test and reference data for the various examples in the \pkg{SingleR} package.
#'
#' @param ngroups Integer scalar specifying the number of groups.
#' @param nreps Integer scalar specifying the number of replicates per group.
#' @param ngenes Integer scalar specifying the number of genes in the dataset.
#' @param prop Numeric scalar specifying the proportion of genes that are DE between groups.
#' @param mock.ref A \linkS4class{SummarizedExperiment} object produced by \code{.mockRefData}.
#' @param ncells Integer scalar specifying the number of cells to simulate.
#'
#' @details
#' This functions are simply provided to simulate some data in the Examples of the documentation.
#' The simulations are very simple and should not be used for performance comparisons.
#'
#' @return 
#' Both functions return a \linkS4class{SummarizedExperiment} object containing simulated counts in the \code{counts} assay,
#' with the group assignment of each sample in the \code{"label"} field of the \code{\link{colData}}.
#' 
#' @author Aaron Lun
#' @examples
#' ref <- .mockRefData()
#' test <- .mockTestData(ref)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom stats rpois rnorm
#' @rdname mockData
.mockRefData <- function(ngroups=5, nreps=4, ngenes=1000, prop=0.5) {
    nmarkers <- ngenes*prop
    nmarkers.per.group <- ceiling(nmarkers/ngroups)

    means <- matrix(0, ncol=ngroups, nrow=ngenes,
        dimnames=list(
            sprintf("GENE_%i", seq_len(ngenes)),
            LETTERS[seq_len(ngroups)]
        )
    )

    counter <- 0L
    for (i in seq_len(ngroups)) {
        means[counter + seq_len(nmarkers.per.group),i] <- rnorm(nmarkers.per.group)
        counter <- counter + nmarkers.per.group
    }

    g <- rep(colnames(means), each=nreps)
    mat <- matrix(rpois(1000*length(g), lambda=10*2^means[,g]), ncol=length(g))
    rownames(mat) <- rownames(means)

    SummarizedExperiment(
        list(counts=mat),
        colData=DataFrame(label=g),
        metadata=list(means=means)
    )
}

#' @export
#' @importFrom S4Vectors metadata
#' @rdname mockData
.mockTestData <- function(mock.ref, ncells=100) {
    means <- metadata(mock.ref)$means

    g <- sample(colnames(means), ncells, replace=TRUE)
    mat <- matrix(rpois(nrow(means)*ncells, lambda=2^means[,g]), ncol=ncells)
    rownames(mat) <- rownames(means)

    SummarizedExperiment(
        list(counts=mat),
        colData=DataFrame(label=g)
    )
}
