#' Rebuild the index
#'
#' Rebuild the index (or indices), typically after restarting the R session.
#' This is because the indices are held in external memory and are not serialized correctly by R.
#'
#' @param trained List containing the output of \code{\link{trainSingleR}},
#' possibly after some operations that invalidate the indices.
#' @param num.threads Integer specifying the number of threads to use for training.
#'
#' @return \code{trained} is returned with valid indices.
#' If it already had valid indices, this function is a no-op.
#' 
#' @author Aaron Lun
#' @examples
#' # Making up the training set.
#' ref <- .mockRefData()
#' ref <- scrapper::normalizeRnaCounts.se(ref)
#' trained <- trainSingleR(ref, ref$label)
#' trained$built # a valid address
#' 
#' # Saving and reloading the index.
#' tmp <- tempfile(fileext=".rds")
#' saveRDS(trained, file=tmp)
#' reloaded <- readRDS(tmp)
#' reloaded$built # not valid anymore
#'
#' rebuilt <- rebuildIndex(reloaded)
#' rebuilt$built # back to validity
#'
#' @export
rebuildIndex <- function(trained, num.threads=1) {
    if (.is_solo(trained)) {
        trained <- .rebuild_index(trained, num.threads = num.threads)
    } else {
        for (i in seq_along(trained)) {
            trained[[i]] <- .rebuild_index(trained[[i]], num.threads = num.threads)
        }
    }
    trained
}

.rebuild_index <- function(trained, num.threads) {
    if (!is(trained$built, "externalptr") || !is_valid_built(trained$built)) {
        trained$built <- .build_index(
            ref=trained$ref, 
            markers=trained$markers$full, 
            labels=trained$labels$full, 
            ulabels=trained$labels$unique, 
            test.genes=trained$options$test.genes,
            num.threads=num.threads) 
    }
    trained
}
