#' @importFrom beachmat getExecutor
.onLoad <- function(libname, pkgname) {
    set_executor(getExecutor())
}
