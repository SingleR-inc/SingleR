#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedMatrixStats rowAnyNAs
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom DelayedArray DelayedArray
.to_clean_matrix <- function(x, assay.type, check.missing, msg="x") {
    if (is.null(rownames(x)) && nrow(x)) { # zero-length matrices have NULL dimnames.
        stop(sprintf("'%s' must have row names", msg))
    }
    if (is(x, "SummarizedExperiment")) {
        x <- assay(x, i=assay.type)
    }

    # Stripping out genes with NA's from 'x'.
    if (check.missing) {
        discard <- rowAnyNAs(DelayedArray(x))
        if (any(discard)) {
            warning(sprintf("'%s' contains rows with missing values", msg))
            x <- x[!discard,,drop=FALSE]
        }
    }

    x
}

.ensure_named <- function(results) {
    if (is.null(rownames(results))) {
        rownames(results) <- seq_len(nrow(results))
    }
    results
}

.name_unless_NULL <- function(target, names) {
    if (!is.null(target)) {
        names(target) <- names
    }
    target
}

.values_title <- function(is.combined, ref.use, ref.names, value.name){
    if (ref.use==0) {
        if (is.combined) {
            front <- "Combined "
        } else {
            front <- ""
        }
    } else {
        front <- paste0(ref.names[ref.use], " ")
    }
    paste0(front, value.name)
}

#' @importFrom DelayedArray is_sparse
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.realize_references <- function(ref) {
    lapply(ref, function(x) {
        if (is_sparse(x)) {
            as(x, "dgCMatrix")
        } else {
            as.matrix(x)
        }
    })
}
