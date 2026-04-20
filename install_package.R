#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  script_path <- normalizePath("install_package.R", winslash = "/", mustWork = FALSE)
}
project_root <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "DESCRIPTION"))) {
  stop("Please run this script from the Kp-VirTracer project root.", call. = FALSE)
}

if (!requireNamespace("remotes", quietly = TRUE) && !requireNamespace("devtools", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

if (requireNamespace("remotes", quietly = TRUE)) {
  remotes::install_local(project_root, upgrade = "never", force = TRUE)
} else {
  devtools::install(project_root, upgrade = "never", force = TRUE)
}

message("[INFO] KpVirTracer installed successfully.")
message("[INFO] Run examples:")
message("  Rscript run_kpvirtracer.R --help")
message("  R -q -e \"library(KpVirTracer); run_kp_virtracer(input='input_dir', output='output_dir')\"")
