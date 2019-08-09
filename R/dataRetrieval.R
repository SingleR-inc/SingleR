########################################################################
### Functions for the retrieval of individual data sets from ExpHub
########################################################################

#' Obtain the HPCA data
#'
#' Download and cache the normalized expression values of the data stored in
#' the Human Primary Cell Atlas. The data will be downloaded from ExperimentHub,
#' returning a \linkS4class{SummarizedExperiment} object for further use.
#'
#' @details
#' This function provides normalized expression values 713 microarray sampes of 
#' the Human Primary Cell Atlas (HPCA) (Mabbott et al., 2013).
#' These 713 samples were processed and normalized as described in Aran, Looney &
#' Liu et al. (2019) and each sample has been assigned to one of 38 main cell types
#' and 169 subtypes.
#' The cell type labels are stored in the colData of the returned \linkS4class{SummarizedExperiment}.
#'
#'
#' @return A \linkS4class{SummarizedExperiment} object.
#'
#' @author Friederike Duendar
#'
#' @references
#' Mabbott et al. (2013).
#' An expression atlas of human primary cells: Inference of gene function from coexpression networks.
#' \emph{BMC Genomics}. doi: 10.1186/1471-2164-14-632
#' 
#' Processing described in Aran, Looney & Liu et al. (2019). 
#' \emph{Nature Immunology} 20, 163–172. doi: 10.1038/s41590-018-0276-y
#' 
#' @examples
#' hpca <- refData_hpca()
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
HumanPrimaryCellAtlasData <- function() {
    version <- "1.0.0"
    se <- .create_se(file.path("hpca", version), has.rowdata=FALSE)
}


#####################################################################
### Helper function
#####################################################################

#' Create a SummarizedExperiment object using data from ExperimentHub
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
