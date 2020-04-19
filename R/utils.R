#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedMatrixStats rowAnyNAs
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom DelayedArray DelayedArray
.to_clean_matrix <- function(x, assay.type, check.missing, msg="x") {
    if (is.null(rownames(x))) {
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

.grab_results <- function(results, index) {
    if (index == 0 || is.null(results$orig.results)) {
        return(results)
    } else {
        return(results$orig.results[[index]])
    }
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

.values_title <- function(results, scores.use, val.name){
    target_bit <- switch(as.character(scores.use==0),
        "TRUE" = ifelse(is.null(results$orig.results), "", "Final"),
        "FALSE" = paste0("Ref #", scores.use))
    paste(target_bit, val.name)
}

.calls_title <- function(results, calls.use, val.name, show = "blank", scores.use = 0){
    if (show == "delta.next") {
        calls.use <- scores.use
    }

    if (!is.null(results$orig.results)) {
        if (calls.use == 0) {
            return(paste0("Final ", val.name))
        } else {
            return(paste0("Ref #", calls.use, " ", val.name))
        }
    } else {
        return(val.name)
    }
}
