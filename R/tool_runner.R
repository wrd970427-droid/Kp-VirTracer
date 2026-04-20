#' Run External Command
#'
#' @param command Executable path.
#' @param args Character vector of arguments.
#' @param wd Working directory.
#' @param error_on_status Whether to fail when non-zero.
#' @return Processx result object.
run_external_tool <- function(command, args = character(), wd = ".", error_on_status = FALSE) {
  processx::run(
    command = command,
    args = args,
    wd = wd,
    echo = FALSE,
    echo_cmd = FALSE,
    error_on_status = error_on_status
  )
}
