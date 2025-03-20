
#' Install Python Dependencies for SignalingProfiler
#'
#' This function installs the required Python dependencies for
#' the SignalingProfiler package within a specified virtual environment.
#' If the default environment name ("r-signalingprofiler") is used and it already exists,
#' the function removes and recreates it before installing dependencies.
#'
#' @param envname "r-signalingprofiler"
#' @param new_env Logical. If `TRUE` (default when `envname` is `"r-signalingprofiler"`),
#'   the existing virtual environment will be removed before installation.
#'
#' @return This function does not return a value; it installs Python dependencies.
#'
#' @details
#' his function ensures that the required Python dependencies (`requests` and `pandas`)
#' are installed in the specified virtual environment. If `new_env` is `TRUE`,
#' the function first removes any existing virtual environment with the same name before proceeding.
#'
#' @export
#'
#' @examples
#' install_sp_py()
#' install_sp_py(envname = "custom-env", new_env = FALSE)
install_sp_py <- function(envname = "r-signalingprofiler",
                          new_env = identical(envname, "r-signalingprofiler")) {

  if(new_env && reticulate::virtualenv_exists(envname)){
    reticulate::virtualenv_remove(envname)
  }

  reticulate::py_install("requests", envname = envname, ...)
  reticulate::py_install("pandas", envname = envname, ...)
}
