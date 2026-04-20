#' Check Required External Tools
#'
#' @param cfg Configuration list.
#' @return Named character vector of resolved executables.
resolve_tool_path <- function(path_or_bin) {
  if (is.null(path_or_bin) || is.na(path_or_bin) || path_or_bin == "") {
    return("")
  }
  if (file.exists(path_or_bin)) {
    return(path_or_bin)
  }
  Sys.which(path_or_bin)
}

check_required_tools <- function(cfg) {
  check_external_tools(cfg)
}

check_external_tools <- function(cfg) {
  req <- unlist(cfg$tools, use.names = TRUE)
  resolved <- vapply(req, resolve_tool_path, character(1))
  missing <- names(resolved)[resolved == ""]
  if (length(missing) > 0) {
    msg <- glue::glue("Missing required tools: {paste(missing, collapse = ', ')}")
    if (isTRUE(cfg$runtime$strict_tools)) {
      stop(msg, call. = FALSE)
    } else {
      log_warn(msg)
    }
  }
  resolved
}
