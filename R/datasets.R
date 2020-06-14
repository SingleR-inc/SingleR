#' Reference dataset extractors
#'
#' These dataset getter functions are deprecated as they have been migrated to the \pkg{celldex} package
#' for more general use throughout the Bioconductor package ecosystem.
#'
#' @param ... Further arguments to pass to the \pkg{celldex} function of the same name.
#'
#' @return A \linkS4class{SummarizedExperiment} object containing the reference dataset.
#'
#' @author Aaron Lun
#' 
#' @name datasets
NULL

#' @export
#' @rdname datasets
HumanPrimaryCellAtlasData <- function(...) {
    .Deprecated(new="celldex::HumanPrimaryCellAtlasData")
    celldex::HumanPrimaryCellAtlasData(...)
}

#' @export
#' @rdname datasets
BlueprintEncodeData <- function(...) {
    .Deprecated(new="celldex::BlueprintEncodeData")
    celldex::BlueprintEncodeData(...)
}

#' @export
#' @rdname datasets
ImmGenData <- function(...) {
    .Deprecated(new="celldex::ImmGenData")
    celldex::ImmGenData(...)
}

#' @export
#' @rdname datasets
MouseRNAseqData <- function(...) {
    .Deprecated(new="celldex::MouseRNAseqData")
    celldex::MouseRNAseqData(...)
}

#' @export
#' @rdname datasets
DatabaseImmuneCellExpressionData <- function(...) {
    .Deprecated(new="celldex::DatabaseImmuneCellExpressionData")
    celldex::DatabaseImmuneCellExpressionData(...)
}

#' @export
#' @rdname datasets
NovershternHematopoieticData <- function(...) {
    .Deprecated(new="celldex::NovershternHematopoieticData")
    celldex::NovershternHematopoieticData(...)
}

#' @export
#' @rdname datasets
MonacoImmuneData <- function(...) {
    .Deprecated(new="celldex::MonacoImmuneData")
    celldex::MonacoImmuneData(...)
}
