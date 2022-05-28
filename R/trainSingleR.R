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
#' Alternatively, a list or \linkS4class{List} of SummarizedExperiment objects or numeric matrices containing multiple references,
#' in which case the row names are expected to be the same across all objects.
#' @param labels A character vector or factor of known labels for all samples in \code{ref}.
#' 
#' Alternatively, if \code{ref} is a list, \code{labels} should be a list of the same length.
#' Each element should contain a character vector or factor specifying the label for the corresponding entry of \code{ref}.
#' @param genes A string specifying the feature selection method to be used, see Details.
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
#' @param sd.thresh A numeric scalar specifying the minimum threshold on the standard deviation per gene.
#' Only used when \code{genes="sd"}.
#' @param de.method String specifying how DE genes should be detected between pairs of labels.
#' Defaults to \code{"classic"}, which sorts genes by the log-fold changes and takes the top \code{de.n}.
#' Setting to \code{"wilcox"} or \code{"t"} will use Wilcoxon ranked sum test or Welch t-test between labels, respectively,
#' and take the top \code{de.n} upregulated genes per comparison.
#' @param de.n An integer scalar specifying the number of DE genes to use when \code{genes="de"}.
#' If \code{de.method="classic"}, defaults to \code{500 * (2/3) ^ log2(N)} where \code{N} is the number of unique labels.
#' Otherwise, defaults to 10.
#' @param de.args Named list of additional arguments to pass to \code{\link[scran]{pairwiseTTests}} or \code{\link[scran]{pairwiseWilcox}} when \code{de.method="wilcox"} or \code{"t"}.
#' @param aggr.ref Logical scalar indicating whether references should be aggregated to pseudo-bulk samples for speed, see \code{\link{aggregateReference}}.
#' @param aggr.args Further arguments to pass to \code{\link{aggregateReference}} when \code{aggr.ref=TRUE}.
#' @param recompute Logical scalar indicating whether to set up indices for later recomputation of scores,
#' when \code{ref} contains multiple references from which the individual results are to be combined.
#' (See the difference between \code{\link{combineCommonResults}} and \code{\link{combineRecomputedResults}}.)
#' @param assay.type An integer scalar or string specifying the assay of \code{ref} containing the relevant expression matrix,
#' if \code{ref} is a \linkS4class{SummarizedExperiment} object (or is a list that contains one or more such objects).
#' @param check.missing Logical scalar indicating whether rows should be checked for missing values (and if found, removed).
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for building nearest neighbor indices.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' @param restrict A character vector of gene names to use for marker selection.
#' By default, all genes in \code{ref} are used.
#'
#' @return 
#' For a single reference, a \linkS4class{List} is returned containing:
#' \describe{
#' \item{\code{common.genes}:}{A character vector of all genes that were chosen by the designated feature selection method.}
#' \item{\code{nn.indices}:}{A List of \linkS4class{BiocNeighborIndex} objects containing pre-constructed neighbor search indices.}
#' \item{\code{original.exprs}:}{A List of numeric matrices where each matrix contains all cells for a particular label.}
#' \item{\code{search}:}{A List of additional information on the feature selection, for use by \code{\link{classifySingleR}}.
#' This includes \code{mode}, a string containing the selection method;
#' \code{args}, method-specific arguments that can be re-used during classification;
#' and \code{extras}, method-specific structures that can be re-used during classification.}
#' }
#'
#' For multiple references, a List of Lists is returned where each internal List corresponds to a reference in \code{ref} and has the same structure as described above.
#'
#' @details
#' This function uses a training data set to select interesting features and construct nearest neighbor indices in rank space.
#' The resulting objects can be re-used across multiple classification steps with different test data sets via \code{\link{classifySingleR}}.
#' This improves efficiency by avoiding unnecessary repetition of steps during the downstream analysis.
#' 
#' Several options are available for feature selection:
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
#' If \code{genes="de"} or \code{"sd"}, the expression values are expected to be log-transformed and normalized.
#'
#' If \code{restrict} is specified, \code{ref} is subsetted to only include the rows with names that are in \code{restrict}.
#' Marker selection and all subsequent classification will be performed using this restrictive subset of genes.
#' This can be convenient for ensuring that only appropriate genes are used (e.g., not pseudogenes or predicted genes).
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
#' Similarly, for manual \code{genes}, \code{de.n} and \code{sd.thresh} have no effect.
#'
#' @section Dealing with multiple references:
#' The default \pkg{SingleR} policy for dealing with multiple references is to perform the classification for each reference separately and combine the results (see \code{?\link{combineRecomputedResults}} for an explanation).
#' To this end, if \code{ref} is a list with multiple references, marker genes are identified separately within each reference when \code{genes="de"} or \code{"sd"}.
#' Rank calculation and index construction is then performed within each reference separately.
#'
#' Alternatively, \code{genes} can still be used to explicitly specify marker genes for each label in each of multiple references.
#' This is achieved by passing a list of lists to \code{genes},
#' where each inner list corresponds to a reference in \code{ref} and can be of any format described in \dQuote{Custom feature specification}.
#' Thus, it is possible for \code{genes} to be - wait for it - a list (per reference) of lists (per label) of lists (per label) of character vectors.
#'
#' If \code{recompute=TRUE}, the output is exactly equivalent to running \code{trainSingleR} on each reference separately.
#' If \code{recompute=FALSE}, \code{trainSingleR} is also run each reference but the difference is that the final \code{common} set of genes consists of the union of common genes across all references.
#' This is necessary to ensure that correlations are computed from the same set of genes across reference and are thus reasonably comparable in \code{\link{combineCommonResults}}.
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
#' \code{\link{combineCommonResults}} and \code{\link{combineRecomputedResults}}, to combine results from multiple references.
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
#' trained$nn.indices
#' length(trained$common.genes)
#'
#' # Alternatively, computing and supplying a set of label-specific markers.
#' by.t <- scran::pairwiseTTests(assay(ref, 2), ref$label, direction="up")
#' markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)
#' trained <- trainSingleR(ref, ref$label, genes=markers)
#' length(trained$common.genes)
#' 
#' @export
#' @importFrom BiocNeighbors KmknnParam 
#' @importFrom S4Vectors List isSingleString metadata metadata<-
#' @importFrom BiocParallel SerialParam bpisup bpstart bpstop
trainSingleR <- function(
    ref, 
    labels, 
    genes="de", 
    sd.thresh=1, 
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
    BNPARAM=KmknnParam(), 
    BPPARAM=SerialParam()) 
{
    de.method <- match.arg(de.method)

    if (solo <- !.is_list(ref)) {
        ref <- list(ref)
        labels <- list(labels)
        if (!is.character(genes)) {
            genes <- list(genes)
        }
    }

    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    ref <- lapply(ref, FUN=.to_clean_matrix, assay.type=assay.type, 
        check.missing=check.missing, msg="ref", BPPARAM=BPPARAM)

    # Cleaning the genes.
    gns <- lapply(ref, rownames)
    if (length(unique(gns))!=1L) {
        stop("row names are not identical across references")
    }

    if (!is.null(restrict)) {
        keep <- gns[[1]] %in% restrict
        ref <- lapply(ref, FUN="[", i=keep, , drop=FALSE)
    }

    if (isSingleString(genes)) {
        genes <- rep(genes, length(ref))
    } else if (length(genes)!=length(ref)) {
        stop("list-like 'genes' should be the same length as 'ref'")
    }

    # Cleaning the labels.
    labels <- lapply(labels, as.character)
    if (length(labels)!=length(ref)) {
        stop("lists in 'labels' and 'ref' should be of the same length")
    }

    for (l in seq_along(labels)) {
        keep <- !is.na(labels[[l]])
        if (!all(keep)) {
            labels[[l]] <- labels[[l]][keep]
            ref[[l]] <- ref[[l]][,keep,drop=FALSE]
        }
    }

    gene.info <- mapply(FUN=.identify_genes, ref=ref, labels=labels, genes=genes,
        MoreArgs=list(sd.thresh=sd.thresh, de.method=de.method, de.n=de.n, de.args=de.args, BPPARAM=BPPARAM),
        SIMPLIFY=FALSE)

    if (!solo && !recompute) {
        .Deprecated("'recompute = FALSE'")
    }
    output <- mapply(FUN=.build_trained_index, ref=ref, labels=labels, markers=gene.info,
        MoreArgs = list(aggr.ref=aggr.ref, aggr.args=aggr.args, BPPARAM=BPPARAM, approximate = approximate), 
        SIMPLIFY=FALSE)

    if (solo) {
        output[[1]]
    } else {
        final <- List(output)
        metadata(final)$recompute <- recompute
        final
    }
}

.identify_genes <- function(ref, labels, genes="de", sd.thresh=1, de.method="classic", de.n=NULL, de.args=list(), BPPARAM=BPPARAM) {
    if (length(labels)!=ncol(ref)) {
        stop("number of labels must be equal to number of cells")
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

#' @importFrom S4Vectors List
#' @importFrom SummarizedExperiment assay
.build_trained_index <- function(ref, labels, markers, aggr.ref, aggr.args, search.info, approximate = FALSE, BPPARAM = SerialParam()) {
    if (aggr.ref) {
        aggr <- do.call(aggregateReference, c(list(ref=quote(ref), label=labels, check.missing=FALSE, BPPARAM=BPPARAM), aggr.args))
        ref <- assay(aggr)
        labels <- aggr$label
    }

    if (anyNA(labels)) {
        keep <- !is.na(labels)
        ref <- ref[,keep,drop=FALSE]
        labels <- labels[keep]
    }
    ulabels <- .get_levels(labels)

    original.markers <- markers
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

    built <- prebuild(ref, match(labels, ulabels) - 1L, markers, approximate = approximate)

    List(
        built = built,
        ref = ref,
        labels = list(full = labels, unique = ulabels),
        markers = list(full = original.markers, unique = rownames(ref)[get_subset(built) + 1])
    )
}

.get_levels <- function(labels) sort(unique(labels))

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

        pairwise <- do.call(FUN, c(list(x=ref, groups=labels, direction="up", BPPARAM=BPPARAM), de.args))
        if (is.null(de.n)) {
            de.n <- 10
        }

        collected <- scran::getTopMarkers(pairwise$statistics, pairwise$pairs, n=de.n)
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
