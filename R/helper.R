#' Retrieve SingleR reference dataset from Github
#'
#' \code{getReferenceDataset} retrieves \code{SingleR} reference datasets from
#' Github. 
#'
#' @details
#' Reference datasets available for retrieval by this function include:
#' \describe{
#'     \item{"hpca"}{Human Primary Cell Atlas (HPCA): a collection of Gene 
#'         Expression Omnibus (GEO datasets), which contains 713 microarray 
#'         samples classified to 38 main cell types and further annotated to 169
#'         subtypes.}
#'     \item{"blueprint_encode"}{Blueprint + ENCODE datasets: Blueprint 
#'         Epigenomics, 144 RNA-seq pure immune samples annotated to 28 cell 
#'         types. ENCODE, 115 RNA-seq pure stroma and immune samples annotated 
#'         to 17 cell types. Altogether, 259 samples with 43 cell types.}
#'     \item{"immgen"}{Immunological Genome Project (ImmGen): 830 microarray 
#'         samples, which we classified to 20 main cell types and further
#'         annotated to 253 subtypes.}
#'     \item{"mouse.rnaseq"}{A dataset of 358 mouse RNA-seq samples annotated to
#'         28 cell types. This dataset was collected, processed, and shared 
#'         courtesy of Bérénice Benayoun. This data set is especially useful for 
#'         brain-related samples.}
#' }
#'
#'
#' @param dataset String for reference dataset to retrieve. 
#' @return A \linkS4class{List} containing:
#' \describe{
#'     \item{\code{data}:}{Numeric matrix containing expression values for 
#'         all genes for each sample.}
#'     \item{\code{types}:}{A character vector containing the specific cell type
#'         for each sample.}
#'     \item{\code{main_types}:}{A character vector containing the broad cell 
#'         type for each sample.}
#'     \item{\code{name}:}{Name of the reference dataset.}
#'     \item{\code{sd.thres}:}{Currently ignored. Numeric scalar to use as the 
#'         standard deviation threshold for selecting variable genes across all 
#'         samples in the reference dataset. See \code{\link{trainSingleR}}.}
#'     \item{\code{de.genes}:}{A list of lists of character vectors containing 
#'         DE genes between pairs of \code{types} labels. See 
#'         \code{\link{trainSingleR}}.}
#'     \item{\code{de.genes.main}:}{A list of lists of character vectors 
#'         containing DE genes between pairs of \code{main_types} labels. 
#'         See \code{\link{trainSingleR}}.}
#' }
#'
#' @author Jared Andrews
#'
#' @examples
#'
#' immgen <- getReferenceDataset(dataset = "immgen")
#' mouse.rnaseq <- getReferenceDataset(dataset = "mouse.rnaseq")
#'
#' @export
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' 
getReferenceDataset <- function(dataset = c("hpca", "blueprint_encode", 
    "immgen", "mouse.rnaseq")) {
    dataset <- match.arg(dataset)
    if (length(dataset) > 1) {
        stop("No dataset specified.")
    }

    full.url <- sprintf("https://github.com/dviraran/SingleR/blob/master/data/%s.rda?raw=true", dataset)

    bfc <- BiocFileCache(ask=FALSE)
    ref <- bfcrpath(bfc, full.url)

    env <- new.env()
    load(ref, envir = env)
    ref.set <- get(dataset, envir = env)

    return(ref.set)
}