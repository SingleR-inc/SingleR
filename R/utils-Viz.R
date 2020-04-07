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
