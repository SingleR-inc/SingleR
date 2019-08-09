#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SummarizedExperiment SummarizedExperiment
.create_se <- function(dataset, hub = ExperimentHub(), assays="normcounts",
                       has.rowdata=TRUE, has.coldata=TRUE, suffix=NULL) {
    host <- file.path("SingleR", dataset)
    if (is.null(suffix)) {
        suffix <- ""
    } else {
        suffix <- paste0("-", suffix)
    }
    
    all.assays <- list()
    for (a in assays) {
        all.assays[[a]] <- hub[hub$rdatapath==file.path(host, sprintf("%s%s.rds", a, suffix))][[1]]
    }
    
    args <- list()
    if (has.coldata) {
        args$colData <- hub[hub$rdatapath==file.path(host, sprintf("coldata%s.rds", suffix))][[1]]
    }
    if (has.rowdata) {
        args$rowData <- hub[hub$rdatapath==file.path(host, sprintf("rowdata%s.rds", suffix))][[1]]
    }
    
    do.call(SummarizedExperiment, c(list(assays=all.assays), args))
}
