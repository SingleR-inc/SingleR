#' Retrieve SingleR reference datasets from Github
#'
#' \code{getSinglerRef} retrieves \code{SingleR} reference datasets from Github. 
#'
#' @param dataset String referring to SingleR reference dataset to retrieve. 
#'     Options include:
#'     \itemize{
#'       \item "hpca"
#'       \item "blueprint_encode"
#'       \item "immgen"
#'       \item "mouse.rnaseq"
#'     }
#' @return SingleR reference dataset.
#'
#' @author Jared Andrews
#'
#' @examples
#'
#' immgen <- getSinglerRef(dataset = "immgen")
#' mouse.rnaseq <- getSinglerRef(dataset = "mouse.rnaseq")
#'
#' @export
#' @importFrom curl curl_download
#' 
getSinglerRef <- function(dataset = "hpca") {
	  opts <- c("hpca", "blueprint_encode", "immgen", "mouse.rnaseq")
	  if (!(dataset %in% opts)) {
		    stop("Invalid SingleR reference set. Valid options are 'hpca',", 
			      "'blueprint_encode', 'immgen', or 'mouse.rnaseq'.")
	  }
	  base.url <- "https://github.com/dviraran/SingleR/blob/master/data/"
	  url.tail <- ".rda?raw=true"

	  full.url <- paste0(base.url, dataset, url.tail)

	  bfc <- BiocFileCache(ask=FALSE)
	  ref <- bfcrpath(bfc, full.url)

	  load(ref, envir = environment())
	  ref.set <- get(dataset, envir = environment())
	
	  return(ref.set)
}