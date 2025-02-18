
.onLoad <- function(libname, pkgname) {
  reticulate::use_virtualenv("r-signalingprofiler", required = FALSE)
}
