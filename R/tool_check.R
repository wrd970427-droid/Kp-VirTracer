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

#' Safe lookup in named tool vector/list (missing names -> NULL, no error).
get_tool <- function(tools, name) {
  if (is.null(tools) || is.null(name) || !nzchar(name)) {
    return(NULL)
  }
  nms <- names(tools)
  if (is.null(nms) || !(name %in% nms)) {
    return(NULL)
  }
  val <- unname(tools[[name]])
  if (length(val) == 0 || is.null(val) || is.na(val) || !nzchar(val)) {
    return(NULL)
  }
  as.character(val[[1]])
}

check_required_tools <- function(cfg) {
  check_external_tools(cfg)
}

check_external_tools <- function(cfg) {
  req <- unlist(cfg$tools, use.names = TRUE)
  resolved <- vapply(req, resolve_tool_path, character(1))
  optional_tools <- c("copla")
  missing <- names(resolved)[resolved == ""]
  missing_required <- setdiff(missing, optional_tools)
  missing_optional <- intersect(missing, optional_tools)

  if (length(missing_required) > 0) {
    msg <- glue::glue("Missing required tools: {paste(missing_required, collapse = ', ')}")
    if (isTRUE(cfg$runtime$strict_tools)) {
      stop(msg, call. = FALSE)
    } else {
      log_warn(msg)
    }
  }
  if (length(missing_optional) > 0) {
    log_warn(glue::glue("Optional tools not found: {paste(missing_optional, collapse = ', ')}"))
  }
  resolved
}
