#' @importFrom beachmat initializeCpp tatami.row.nan.counts
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedMatrixStats rowAnyNAs
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom DelayedArray DelayedArray
.to_clean_matrix <- function(x, assay.type, check.missing, msg="x", num.threads=1) {
    if (is.null(rownames(x)) && nrow(x)) { # zero-length matrices have NULL dimnames.
        stop(sprintf("'%s' must have row names", msg))
    }

    if (is(x, "SummarizedExperiment")) {
        x <- assay(x, i=assay.type)
    }

    if (.is_data_frame(x)) {
        x <- as.matrix(x)
        if (!is.numeric(x)) {
            stop("failed to convert data.frame into a numeric matrix")
        }
    }

    # Stripping out genes with NA's from 'x'.
    if (check.missing) {
        ptr <- initializeCpp(x, .check.na=TRUE)
        keep <- tatami.row.nan.counts(ptr, num.threads=num.threads) == 0
        if (!all(keep)) {
            # Returning a DelayedArray to avoid making an actual subset.
            x <- DelayedArray(x)[keep,,drop=FALSE]
        }
    }

    x
}

.create_intersection <- function(test, reference) {
    # Effectively an NA-safe intersect() that preserves ordering in 'test'.
    common <- test[test %in% reference]
    common <- common[!is.na(common)]
    common <- common[!duplicated(common)]

    # match() takes the first occurrence, consistent with internal behavior in singlepp.
    list(
        test = match(common, test),
        reference = match(common, reference)
    )
}

#' @importFrom methods is
#' @importClassesFrom S4Vectors List
.is_list <- function(val) {
    (is.list(val) || is(val, "List")) && !.is_data_frame(val)
}

#' @importFrom methods is
#' @importClassesFrom S4Vectors DataFrame
.is_data_frame <- function(val) {
    is.data.frame(val) || is(val, "DataFrame")
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
.realize_reference <- function(x) {
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x)
    }
}

.get_num_threads <- function(num.threads, BPPARAM) {
    if (!is.null(BPPARAM)) {
        num.threads <- BiocParallel::bpnworkers(BPPARAM)
    }
    num.threads
}

utils::globalVariables(".data")
