#' Train the SingleR classifier 
#'
#' Train the SingleR classifier on one or more reference datasets with known labels.
#' 
#' @param ref A numeric matrix of expression values where rows are genes and columns are reference samples (individual cells or bulk samples).
#' Each row should be named with the gene name.
#' In general, the expression values are expected to be normalized and log-transformed, see Details.
#' 
#' Alternatively, a \link[SummarizedExperiment]{SummarizedExperiment} object containing such a matrix.
#'
#' Alternatively, a list or \link[S4Vectors]{List} of SummarizedExperiment objects or numeric matrices containing multiple references.
#' @param labels A character vector or factor of known labels for all samples in \code{ref}.
#' 
#' Alternatively, if \code{ref} is a list, \code{labels} should be a list of the same length.
#' Each element should contain a character vector or factor specifying the labels for the columns of the corresponding element of \code{ref}.
#' @param genes A string containing \code{"de"}, indicating that markers should be calculated from \code{ref}.
#' For back compatibility, other string values are allowed but will be ignored with a deprecation warning.
#' 
#' Alternatively, if \code{ref} is \emph{not} a list, \code{genes} can be either:
#' \itemize{
#' \item A list of lists of character vectors containing DE genes between pairs of labels.
#' \item A list of character vectors containing marker genes for each label.
#' }
#' 
#' If \code{ref} \emph{is} a list, \code{genes} can be a list of length equal to \code{ref}.
#' Each element of the list should be one of the two above choices described for non-list \code{ref},
#' containing markers for labels in the corresponding entry of \code{ref}.
#' @param sd.thresh Deprecated and ignored.
#' @param de.method String specifying how DE genes should be detected between pairs of labels.
#' Defaults to \code{"classic"}, which sorts genes by the log-fold changes and takes the top \code{de.n}.
#' Other options are \code{"wilcox"} and \code{"t"}, see Details.
#' Ignored if \code{genes} is a list of markers/DE genes.
#' @param de.n An integer scalar specifying the number of DE genes to use when \code{genes="de"}.
#' If \code{de.method="classic"}, defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
#' Otherwise, defaults to 10.
#' Ignored if \code{genes} is a list of markers/DE genes.
#' @param de.args Named list of additional arguments to pass to \code{\link[scrapper]{scoreMarkers}} when \code{de.method="wilcox"} or \code{"t"}.
#' Ignored if \code{genes} is a list of markers/DE genes.
#' @param aggr.ref Logical scalar indicating whether references should be aggregated to pseudo-bulk samples for speed, see \code{\link{aggregateReference}}.
#' @param aggr.args Further arguments to pass to \code{\link{aggregateReference}} when \code{aggr.ref=TRUE}.
#' @param recompute Deprecated and ignored.
#' @param assay.type An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \link[SummarizedExperiment]{SummarizedExperiment} object (or is a list that contains one or more such objects).
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values.
#' If true and any missing values are found, the rows containing these values are silently removed.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying how the neighbor search index should be constructed.
#' @param approximate Deprecated, use \code{BNPARAM} instead.
#' @param num.threads Integer scalar specifying the number of threads to use for index building.
#' @param hint.sce Boolean indicating whether to print a hint to change \code{de.method=} when any entry of \code{ref} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' It will also suggest setting \code{aggr.ref=TRUE} for greater efficiency when \code{ref} contains 1000 cells or more.
#' @param BPPARAM A \link[BiocParallel]{BiocParallelParam} object specifying how parallelization should be performed when \code{check.missing = TRUE}.
#' @param restrict A character vector of gene names to use for marker selection.
#' By default, all genes in \code{ref} are used.
#' @param test.genes Character vector of the names of the genes in the test dataset, i.e., the row names of \code{test} in \code{\link{classifySingleR}}.
#' If \code{NULL}, it is assumed that the test dataset and \code{ref} have the same genes in the same row order.
#'
#' @return 
#' For a single reference, a \link[S4Vectors]{List} is returned containing:
#' \describe{
#' \item{\code{built}:}{An external pointer to various indices in C++ space.
#' Note that this cannot be serialized and should be removed prior to any \code{\link{saveRDS}} step.}
#' \item{\code{ref}:}{The reference expression matrix.
#' This may have fewer columns than the input \code{ref} if \code{aggr.ref = TRUE}.}
#' \item{\code{markers}:}{A list containing \code{unique}, a character vector of all marker genes used in training;
#' and \code{full}, a list of list of character vectors containing the markers for each pairwise comparison between labels.}
#' \item{\code{labels}:}{A list containing \code{unique}, a character vector of all unique reference labels;
#' and \code{full}, a character vector containing the assigned label for each column in \code{ref}.}
#' }
#'
#' For multiple references, a List of Lists is returned where each internal List corresponds to a reference in \code{ref} and has the same structure as described above.
#'
#' @details
#' This function uses a training data set to select interesting features and construct nearest neighbor indices in rank space.
#' The resulting objects can be re-used across multiple classification steps with different test data sets via \code{\link{classifySingleR}}.
#' This improves efficiency by avoiding unnecessary repetition of steps during the downstream analysis.
#' 
#' The automatic marker detection (\code{genes="de"}) identifies genes that are differentially expressed between pairs of labels in the reference dataset.
#' The expression values are expected to be log-transformed and normalized.
#' For each pair of labels, the top \code{de.n} genes with strongest upregulation in one label are chosen as markers to distinguish it from the other label.
#' The exact ranking depends on the \code{de.method=} argument:
#' \itemize{
#' \item The default \code{de.method="classic"} will use \code{\link{getClassicMarkers}} to compute the median expression for each label and each gene.
#' Then, for each pair of labels, the top \code{de.n} genes with the largest positive differences are chosen as markers to distinguish the first label from the second.
#' This is intended for reference datasets derived from bulk transcriptomic data (e.g., microarrays) with a high density of non-zero values. 
#' It is less effective for single-cell data, where it is not uncommon to have more than 50\% zero counts for a given gene such that the median is also zero for each group.
#' \item \code{de.method="wilcox"} will rank genes based on the area under the curve (AUC) in each pairwise comparison between labels.
#' The top \code{de.n} genes with the largest AUCs above 0.5 are chosen as markers for the first label compared to the second.
#' This is analogous to ranking on significance in the Wilcoxon ranked sum test and is intended for use with single-cell data.
#' The exact calculaton is performed using the \code{\link[scrapper]{scoreMarkers}} function.
#' \item \code{de.method="t"} will rank genes on the Cohen's d in each pairiwse comparison.
#' The top \code{de.n} genes with the largest positive Cohen's d are chosen as markers for the first label compared to the second.
#' This is roughly analogous to ranking on significance in the t-test and is faster than the AUC.
#' The exact calculaton is performed using the \code{\link[scrapper]{scoreMarkers}} function.
#' }
#' Alternatively, users can detect markers externally and pass a list of markers to \code{genes} (see \dQuote{Custom gene specification}).
#'
#' Classification with \code{classifySingleR} assumes that the test dataset contains all marker genes that were detected from the reference.
#' If the test and reference datasets do not have the same genes in the same order, we can set \code{test.genes} to the row names of the test dataset.
#' This will instruct \code{trainSingleR} to only consider markers that are present in the test dataset.
#' Any subsequent call to \code{classifySingleR} will also check that \code{test.genes} is consistent with \code{rownames(test)}.
#'
#' On a similar note, if \code{restrict} is specified, marker selection will only be performed using the specified subset of genes.
#' This can be convenient for ignoring inappropriate genes like pseudogenes or predicted genes.
#' It has the same effect as filtering out undesirable rows from \code{ref} prior to calling \code{trainSingleR}.
#' Unlike \code{test.genes}, setting \code{restrict} does not introduce further checks on \code{rownames(test)} in \code{classifySingleR}.
#'
#' @section Custom gene specification:
#' Rather than relying on the in-built feature selection, users can pass in their own features of interest to \code{genes}.
#' The function expects a named list of named lists of character vectors, with each vector containing the DE genes between a pair of labels.
#' For example:
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
#' Alternatively, \code{genes} can be a named list of character vectors containing per-label markers.
#' For example:
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
#' If \code{genes} is manually passed, \code{ref} can contain the raw counts or any monotonic transformation thereof.
#' There is no need to supply (log-)normalized expression values for the benefit of the automatic marker detection.
#' Similarly, for manual \code{genes}, the values of \code{de.method}, \code{de.n} and \code{sd.thresh} have no effect.
#'
#' Check out the Examples to see how manual \code{genes} can be passed to \code{trainSingleR}.
#'
#' @section Dealing with multiple references:
#' The default \pkg{SingleR} policy for dealing with multiple references is to perform the classification for each reference separately and combine the results 
#' (see \code{?\link{combineRecomputedResults}} for an explanation).
#' To this end, if \code{ref} is a list with multiple references, marker genes are identified separately within each reference if \code{genes = NULL}.
#' Rank calculation and index construction is then performed within each reference separately.
#' The result is identical to \code{lapply}ing over a list of references and runing \code{trainSingleR} on each reference.
#'
#' Alternatively, \code{genes} can still be used to explicitly specify marker genes for each label in each of multiple references.
#' This is achieved by passing a list of lists to \code{genes},
#' where each inner list corresponds to a reference in \code{ref} and can be of any format described in \dQuote{Custom feature specification}.
#' Thus, it is possible for \code{genes} to be - wait for it - a list (per reference) of lists (per label) of lists (per label) of character vectors.
#'
#' @section Aggregating single-cell references:
#' It is generally unnecessary to have single-cell resolution on the reference profiles.
#' We can instead set \code{aggr.ref=TRUE} to aggregate per-cell references into a set of pseudo-bulk profiles using \code{\link{aggregateReference}}.
#' This improves classification speed while using vector quantization to preserve within-label heterogeneity and mitigate the loss of information.
#' Note that any aggregation is done \emph{after} marker gene detection; this ensures that the relevant tests can appropriately penalize within-label variation.
#' Users should also be sure to set the seed as the aggregation involves randomization.
#'
#' @author Aaron Lun, based on the original \code{SingleR} code by Dvir Aran.
#' 
#' @seealso
#' \code{\link{classifySingleR}}, where the output of this function gets used.
#'
#' \code{\link{combineRecomputedResults}}, to combine results from multiple references.
#'
#' \code{\link{rebuildIndex}}, to rebuild the index after external memory is invalidated.
#'
#' @examples
#' # Making up some data for a quick demonstration.
#' ref <- .mockRefData()
#'
#' # Normalizing and log-transforming for automated marker detection.
#' ref <- scrapper::normalizeRnaCounts.se(ref)
#'
#' trained <- trainSingleR(ref, ref$label)
#' trained
#' length(trained$markers$unique)
#'
#' # Alternatively, supplying a custom set of markers from pairwise comparisons.
#' all.labels <- unique(ref$label)
#' custom.markers <- list()
#' for (x in all.labels) {
#'     current.markers <- lapply(all.labels, function(x) sample(rownames(ref), 20))
#'     names(current.markers) <- all.labels
#'     current.markers[[x]] <- character(0)
#'     custom.markers[[x]] <- current.markers
#' }
#' custom.trained <- trainSingleR(ref, ref$label, genes=custom.markers)
#'
#' # Alternatively, supplying a custom set of markers for each label.
#' custom.markers <- list()
#' for (x in all.labels) {
#'     custom.markers[[x]] <- sample(rownames(ref), 20)
#' }
#' custom.trained <- trainSingleR(ref, ref$label, genes=custom.markers)
#' 
#' @export
#' @importFrom S4Vectors List isSingleString metadata metadata<-
#' @importFrom BiocNeighbors defineBuilder AnnoyParam KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors List
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedArray DelayedArray
trainSingleR <- function(
    ref, 
    labels, 
    test.genes=NULL,
    genes="de", 
    sd.thresh=NULL, 
    de.method=c("classic", "wilcox", "t"), 
    de.n=NULL, 
    de.args=list(),
    aggr.ref=FALSE, 
    aggr.args=list(), 
    recompute=TRUE, 
    restrict=NULL,
    assay.type="logcounts", 
    check.missing=TRUE,
    hint.sce=TRUE,
    approximate = FALSE,
    num.threads = bpnworkers(BPPARAM),
    BNPARAM = NULL,
    BPPARAM = SerialParam()
) {
    de.method <- match.arg(de.method)

    if (solo <- !.is_list(ref)) {
        ref <- list(ref)
        labels <- list(labels)
        if (!is.character(genes)) {
            genes <- list(genes)
        }
    }

    if (isSingleString(genes)) {
        genes <- rep(genes, length(ref))
    } else if (length(genes)!=length(ref)) {
        stop("list-like 'genes' should be the same length as 'ref'")
    }

    if (is.null(BNPARAM)) {
        if (approximate) {
            BNPARAM <- AnnoyParam()
        } else {
            BNPARAM <- KmknnParam()
        }
    }

    output <- vector("list", length(ref))
    names(output) <- names(ref)
    for (l in seq_along(ref)) {
        curref <- ref[[l]]
        if (is(curref, "SingleCellExperiment") && hint.sce) {
            hints <- character()
            what <- "SingleCellExperiment"
            if (de.method == "classic") {
                hints <- c(hints, "'de.method = \"t\"' or \"wilcox\"")
            }
            if (ncol(curref) >= 1000) {
                what <- paste("large", what)
                hints <- c(hints, "'aggr.ref = TRUE' for speed")
            }
            if (length(hints)) {
                msg <- sprintf("Detected a %s as the reference dataset, consider setting %s in trainSingleR().", what, paste(hints, collapse=" and "))
                msg <- paste(msg, "If you know better, this hint can be disabled with 'hint.sce=FALSE'.")
                message(paste(strwrap(msg, 80), collapse="\n"))
                hint.sce <- FALSE # only need to print this message once.
            }
        }

        curref <- .to_clean_matrix(curref, assay.type, check.missing, msg="ref", num.threads=num.threads)

        # Removing duplicated names and missing labels.
        if (anyDuplicated(rownames(curref))) {
            keep <- !duplicated(rownames(curref))
            curref <- DelayedArray(curref)[keep,,drop=FALSE]
        }

        curlabels <- as.character(labels[[l]])
        stopifnot(length(curlabels) == ncol(curref))
        keep <- !is.na(curlabels)
        if (!all(keep)) {
            curref <- DelayedArray(curref)[,keep,drop=FALSE]
            curlabels <- curlabels[keep]
        }

        markers <- .identify_genes(
            ref=curref,
            labels=curlabels,
            genes=genes[[l]],
            de.method=de.method,
            de.n=de.n,
            de.args=de.args,
            restrict=restrict,
            test.genes=test.genes,
            num.threads=num.threads
        )

        if (aggr.ref) {
            aggr <- do.call(aggregateReference, c(list(ref=quote(curref), label=curlabels, check.missing=FALSE, num.threads=num.threads), aggr.args))
            curref <- assay(aggr)
            curlabels <- aggr$label
        }

        ulabels <- .get_levels(curlabels)
        built <- .build_index(
            ref=curref,
            labels=curlabels,
            ulabels=ulabels,
            test.genes=test.genes,
            markers=markers,
            BNPARAM=BNPARAM,
            num.threads=num.threads
        )

        output[[l]] <- List(
            built = built,
            ref = curref,
            labels = list(full = curlabels, unique = ulabels),
            markers = list(full = markers, unique = rownames(curref)[get_ref_subset(built) + 1]),
            options = list(BNPARAM = BNPARAM, test.genes = test.genes)
        )
    }

    if (solo) {
        output[[1]]
    } else {
        List(output)
    }
}

.identify_genes <- function(ref, labels, genes, de.method, de.n, test.genes, restrict, de.args, num.threads) {
    if (length(labels)!=ncol(ref)) {
        stop("number of labels must be equal to number of cells")
    }

    # Note that the genes are reported as names rather than indexing, so these
    # row-subsetting operations don't require re-indexing of the output. We use
    # DelayedArrays to avoid making copies of the data.
    if (!is.null(restrict)) {
        ref <- DelayedArray(ref)[rownames(ref) %in% restrict,,drop=FALSE]
    }
    if (!is.null(test.genes)) {
        ref <- DelayedArray(ref)[rownames(ref) %in% test.genes,,drop=FALSE]
    }

    if (.is_list(genes)) {
        is.char <- vapply(genes, is.character, TRUE)
        if (all(is.char)) {
            genes <- .convert_per_label_set(genes)
        } else if (any(is.char)) {
            stop("'genes' must be a list of character vectors or a list of list of vectors")
        }

        genes <- lapply(genes, as.list) # to convert from List of Lists.
        genes <- .validate_de_gene_set(genes, labels)

        # Ensure that the user hasn't supplied genes that aren't available.
        rn <- rownames(ref)
        genes <- lapply(genes, function(l) lapply(l, intersect, rn))
    } else { 
        genes <- match.arg(genes, c("de", "sd", "all"))
        if (genes != "de") {
            .Deprecated(old="genes = \"", genes, "\"")
        } 
        genes <- .get_genes_by_de(ref, labels, de.n=de.n, de.method=de.method, de.args=de.args, num.threads=num.threads)
    }

    genes
}

#' @importFrom beachmat initializeCpp
.build_index <- function(ref, labels, ulabels, markers, test.genes, BNPARAM, num.threads) {
    for (m in seq_along(markers)) {
        current <- markers[[m]]
        for (n in seq_along(current)) {
            idx <- match(current[[n]], rownames(ref))
            if (anyNA(idx)) {
                stop("could not find '", current[[n]][which(is.na(idx))[1]], "' in 'rownames(ref)'")
            }
            current[[n]] <- idx - 1L
        }
        markers[[m]] <- current
    }

    if (is.null(test.genes)) {
        test.genes <- ref.genes <- seq_len(nrow(ref))
    } else {
        intersection <- .create_intersection(test.genes, rownames(ref))
        test.genes <- intersection$test
        ref.genes <- intersection$reference
    }

    builder <- defineBuilder(BNPARAM)
    parsed <- initializeCpp(ref, .check.na=FALSE)
    train_single(
        test_features=test.genes - 1L, 
        ref=parsed,
        ref_features=ref.genes - 1L,
        labels=match(labels, ulabels) - 1L,
        markers=markers,
        builder=builder$builder,
        nthreads=num.threads
    )
}

.get_levels <- function(labels) sort(unique(labels))

# Unfortunately, we can't test for List, because each trained structure is
# also a list; so we just check whether the 'ref' field exists.
.is_solo <- function(trained) !is.null(trained$ref)

#' @importFrom utils head
.get_genes_by_de <- function(ref, labels, de.method="classic", de.n=NULL, de.args=list(), num.threads=1) {
    if (de.method=="classic") {
        return(getClassicMarkers(ref=ref, labels=labels, de.n=de.n, check.missing=FALSE, num.threads=num.threads))
    }

    if (de.method=="t") {
        compute.auc <- FALSE
        compute.cohens.d <- TRUE
        effect.size <- "cohens.d"
    } else {
        compute.auc <- TRUE
        compute.cohens.d <- FALSE 
        effect.size <- "auc"
    }

    if (is.null(de.n)) {
        de.n <- 10
    }

    pairwise <- do.call(
        scrapper::scoreMarkers,
        c(
            list(
                ref,
                groups=labels,
                num.threads=num.threads,
                all.pairwise=de.n,
                compute.group.mean=FALSE,
                compute.group.detected=FALSE,
                compute.delta.detected=FALSE,
                compute.delta.mean=FALSE,
                compute.auc=compute.auc,
                compute.cohens.d=compute.cohens.d
            ),
            de.args
        )
    )[[effect.size]]

    all.labels <- names(pairwise)
    all.genes <- rownames(ref)
    for (g1 in all.labels) {
        for (g2 in all.labels) {
            if (g1 == g2) {
                pairwise[[g1]][[g2]] <- character(0)
            } else {
                pairwise[[g1]][[g2]] <- all.genes[pairwise[[g1]][[g2]]$index]
            }
        }
    }

    pairwise
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
    ulabels <- .get_levels(labels)
    if (!all(ulabels %in% names(genes))) {
        stop("need marker gene information for each label")
    }

    genes <- genes[ulabels]
    for (u in ulabels) {
        if (!all(ulabels %in% names(genes[[u]]))) {
            stop("need marker genes between each pair of labels")
        }
        genes[[u]] <- genes[[u]][ulabels]
    }

    genes
}
