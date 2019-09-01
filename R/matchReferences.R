#' Match labels from two references
#'
#' Match labels from a pair of references, corresponding to the same underlying cell type or state
#' but with differences in nomenclature.
#' 
#' @param ref1,ref2 Numeric matrices of single-cell (usually log-transformed) expression values where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix
#' @param labels1,labels2 A character vector or factor of known labels for all cells in \code{ref1} and \code{ref2}, respectively.
#' @param ... Further arguments to pass to \code{\link{SingleR}}.
#'
#' @return A numeric matrix containing a probability table of mutual assignment.
#' Values close to 1 represent a 1:1 mapping between labels across the two references.
#'
#' @details
#' It is often the case that two references contain the same cell types for the same biological system,
#' but the two sets of labels differ in their nomenclature.
#' This makes it difficult to compare results from different references.
#' It also interferes with attempts to combine multiple datasets to create a larger, more comprehensive reference.
#' 
#' The \code{matchReferences} function attempts to facilitate matching of labels across two reference datasets.
#' It does so by using one of the references (say, \code{ref1}) to assign its labels to the other (\code{ref2}).
#' For each label X in \code{labels2}, we compute the probability of assigning a sample of X to each label Y in \code{labels1}.
#' We also use \code{ref2} to assign labels to \code{ref1}, to obtain the probability of assigning a sample of Y to label X. 
#'
#' We then consider the probability of mutual assignment, i.e., assigning a sample of X to Y \emph{and} a sample of Y to X.
#' This is computed by simply taking the product of the two probabilities mentioned earlier.
#' The output matrix contains mutual assignment probabilities for all pairs of X (rows) and Y (columns).
#'
#' The mutual assignment probabilities are only high if there is a 1:1 mapping between labels.
#' A perfect mapping manifests as probabilities of 1 in the relevant entries of the output matrix.
#' Lower values are expected for ambiguous mappings and near-zero values for labels that are specific to one reference.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{SingleR}}, to do the actual cross-assignment.
#' 
#' @examples
#' example(SingleR, echo=FALSE)
#' test$label <- paste0(test$label, "_X") # modifying the labels.
#' matchReferences(test, sce, labels1=test$label, labels2=sce$label)
#' @export
matchReferences <- function(ref1, ref2, labels1, labels2, ...) {
    first <- SingleR(test=ref1, ref=ref2, labels=labels2, ...)
    second <- SingleR(test=ref2, ref=ref1, labels=labels1, ...)

    f1 <- factor(labels1)
    f2 <- factor(labels2)
    tab1 <- table(f1, factor(first$labels, levels(f2)))
    tab2 <- table(f2, factor(second$labels, levels(f1)))

    tab1 <- tab1/rowSums(tab1)
    tab2 <- tab2/rowSums(tab2)
    output <- tab1 * t(tab2)

    output <- unclass(output)
    dimnames(output) <- unname(dimnames(output))
    output
}
