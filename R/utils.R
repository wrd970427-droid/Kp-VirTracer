# Utility helpers for KpVirTracer.

assert_dir_exists <- function(path, label = "directory") {
  if (!dir.exists(path)) {
    stop(glue::glue("Missing {label}: {path}"), call. = FALSE)
  }
  invisible(TRUE)
}
