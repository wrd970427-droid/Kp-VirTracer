# Logging helpers for KpVirTracer.

log_info <- function(msg) cli::cli_alert_info(msg)
log_warn <- function(msg) cli::cli_alert_warning(msg)
log_error <- function(msg) cli::cli_alert_danger(msg)
