#' Train the SingleR classifier 
#'
#' Train the SingleR classifier on one or more reference datasets with known labels.
#' 
#' @param ref A numeric matrix of expression values where rows are genes and columns are reference samples (individual cells or bulk samples).
#' Each row should be named with the gene name.
#' In general, the expression values are expected to be log-transformed, see Details.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#'
#' Alternatively, a list or \linkS4class{List} of SummarizedExperiment objects or numeric matrices containing multiple references.
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
#' Setting to \code{"wilcox"} or \code{"t"} will use Wilcoxon ranked sum test or Welch t-test between labels, respectively,
#' and take the top \code{de.n} upregulated genes per comparison.
#' Ignored if \code{genes} is a list of markers/DE genes.
#' @param de.n An integer scalar specifying the number of DE genes to use when \code{genes="de"}.
#' If \code{de.method="classic"}, defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
#' Otherwise, defaults to 10.
#' Ignored if \code{genes} is a list of markers/DE genes.
#' @param de.args Named list of additional arguments to pass to \code{\link[scran]{pairwiseTTests}} or \code{\link[scran]{pairwiseWilcox}} when \code{de.method="wilcox"} or \code{"t"}.
#' Ignored if \code{genes} is a list of markers/DE genes.
#' @param aggr.ref Logical scalar indicating whether references should be aggregated to pseudo-bulk samples for speed, see \code{\link{aggregateReference}}.
#' @param aggr.args Further arguments to pass to \code{\link{aggregateReference}} when \code{aggr.ref=TRUE}.
#' @param recompute Deprecated and ignored.
#' @param assay.type An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object (or is a list that contains one or more such objects).
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values.
#' If true and any missing values are found, the rows containing these values are silently removed.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying how the neighbor search index should be constructed.
#' @param approximate Deprecated, use \code{BNPARAM} instead.
#' @param num.threads Integer scalar specifying the number of threads to use for index building.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' Relevant for marker detection if \code{genes = NULL}, aggregation if \code{aggr.ref = TRUE}, and \code{NA} checking if \code{check.missing = TRUE}.
#' @param restrict A character vector of gene names to use for marker selection.
#' By default, all genes in \code{ref} are used.
#' @param test.genes Character vector of the names of the genes in the test dataset, i.e., the row names of \code{test} in \code{\link{classifySingleR}}.
#' If \code{NULL}, it is assumed that the test dataset and \code{ref} have the same genes in the same row order.
#'
#' @return 
#' For a single reference, a \linkS4class{List} is returned containing:
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
#' The automatic marker detection (\code{genes="de"}) identifies genes that are differentially expressed between labels.
#' This is done by identifying the median expression within each label, and computing differences between medians for each pair of labels.
#' For each label, the top \code{de.n} genes with the largest differences compared to another label are chosen as markers to distinguish the two labels.
#' The expression values are expected to be log-transformed and normalized.
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
#' @section Custom feature specification:
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
#' If \code{genes} is manually passed, \code{ref} can be the raw counts or any monotonic transformation thereof.
#' There is no need to supply (log-)normalized expression values for the benefit of the automatic marker detection.
#' Similarly, for manual \code{genes}, the values of \code{de.method}, \code{de.n} and \code{sd.thresh} have no effect.
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
#' @section Note on single-cell references:
#' The default marker selection is based on log-fold changes between the per-label medians and is very much designed with bulk references in mind.
#' It may not be effective for single-cell reference data where it is not uncommon to have more than 50\% zero counts for a given gene such that the median is also zero for each group.
#' Users are recommended to either set \code{de.method} to another DE ranking method, or detect markers externally and pass a list of markers to \code{genes} (see Examples).
#'
#' In addition, it is generally unnecessary to have single-cell resolution on the reference profiles.
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
#' ref <- scuttle::logNormCounts(ref)
#'
#' trained <- trainSingleR(ref, ref$label)
#' trained
#' length(trained$markers$unique)
#'
#' # Alternatively, computing and supplying a set of label-specific markers.
#' by.t <- scran::pairwiseTTests(assay(ref, 2), ref$label, direction="up")
#' markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)
#' trained <- trainSingleR(ref, ref$label, genes=markers)
#' length(trained$markers$unique)
#' 
#' @export
#' @importFrom S4Vectors List isSingleString metadata metadata<-
#' @importFrom BiocNeighbors defineBuilder AnnoyParam KmknnParam
#' @importFrom BiocParallel SerialParam bpisup bpstart bpstop
#' @importFrom beachmat initializeCpp
#' @importFrom S4Vectors List
#' @importFrom SummarizedExperiment assay
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
    approximate = FALSE,
    num.threads = bpnworkers(BPPARAM),
    BNPARAM = NULL,
    BPPARAM = SerialParam()) 
{
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

    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    output <- vector("list", length(ref))
    names(output) <- names(ref)
    for (l in seq_along(ref)) {
        curref <- .to_clean_matrix(ref[[l]], assay.type, check.missing, msg="ref", BPPARAM=BPPARAM)

        curlabels <- as.character(labels[[l]])
        stopifnot(length(curlabels) == ncol(curref))
        keep <- !is.na(curlabels)
        if (!all(keep)) {
            curref <- curref[,keep,drop=FALSE]
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
            BPPARAM=BPPARAM
        )

        if (aggr.ref) {
            aggr <- do.call(aggregateReference, c(list(ref=quote(curref), label=curlabels, check.missing=FALSE, BPPARAM=BPPARAM), aggr.args))
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

.identify_genes <- function(ref, labels, genes, de.method, de.n, test.genes, restrict, de.args, BPPARAM) {
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
        genes <- .get_genes_by_de(ref, labels, de.n=de.n, de.method=de.method, de.args=de.args, BPPARAM=BPPARAM)
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
        stopifnot(nrow(test) == nrow(ref))
        test.genes <- ref.genes <- seq_len(nrow(ref))
    } else {
        intersection <- .create_intersection(test.genes, rownames(ref))
        test.genes <- intersection$test
        ref.genes <- intersection$reference
    }

    builder <- defineBuilder(BNPARAM)
    parsed <- initializeCpp(ref)
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
.get_genes_by_de <- function(ref, labels, de.method="classic", de.n=NULL, de.args=list(), BPPARAM=SerialParam()) {
    if (de.method=="classic") {
        getClassicMarkers(ref=ref, labels=labels, de.n=de.n, check.missing=FALSE, BPPARAM=BPPARAM)
    } else {
        if (de.method=="t") {
            FUN <- scran::pairwiseTTests
        } else {
            FUN <- scran::pairwiseWilcox
        }

        pairwise <- do.call(FUN, c(list(x=ref, groups=labels, direction="up", log.p=TRUE, BPPARAM=BPPARAM), de.args))
        if (is.null(de.n)) {
            de.n <- 10
        }

        collected <- scran::getTopMarkers(pairwise$statistics, pairwise$pairs, n=de.n, pval.field="log.p.value", fdr.field="log.FDR", fdr.threshold=log(0.05))
        lapply(collected, as.list)
    }
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
