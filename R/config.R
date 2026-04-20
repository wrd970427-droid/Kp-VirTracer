#' Load KpVirTracer Configuration
#'
#' Read YAML config and merge with package defaults.
#'
#' @param config Optional user config path.
#' @param threads Thread override.
#' @return A named list.
load_kp_virtracer_config <- function(config = NULL, threads = 4L) {
  default_path <- system.file("extdata", "config.yaml.example", package = "KpVirTracer")
  if (identical(default_path, "")) {
    default_path <- file.path("inst", "extdata", "config.yaml.example")
  }
  default_cfg <- yaml::read_yaml(default_path)
  if (!is.null(config)) {
    user_cfg <- yaml::read_yaml(config)
    default_cfg <- utils::modifyList(default_cfg, user_cfg)
  }
  default_cfg$threads <- as.integer(threads)
  if (is.null(default_cfg$runtime)) {
    default_cfg$runtime <- list()
  }
  if (is.null(default_cfg$runtime$strict_tools)) {
    default_cfg$runtime$strict_tools <- FALSE
  }
  default_cfg
}
