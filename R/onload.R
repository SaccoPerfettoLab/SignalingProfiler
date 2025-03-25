
#' Package Load Hook
#'
#' This function is executed when the package is loaded. It sets up the reticulate
#' virtual environment `"r-signalingprofiler"` if available.
#'
#' @param libname The name of the library where the package is installed.
#' @param pkgname The name of the package being loaded.
#'
#' @details
#' This function is called internally when the package is loaded. It ensures that
#' the `"r-signalingprofiler"` virtual environment is used, but does not
#' enforce its availability (`required = FALSE`).
#'
#' If the virtual environment does not exist, it will not trigger an error,
#' allowing the package to function without interruption.
#'
#' @keywords internal
#' @noRd
#'
.onLoad <- function(libname, pkgname) {
  reticulate::use_virtualenv("r-signalingprofiler", required = FALSE)
}
