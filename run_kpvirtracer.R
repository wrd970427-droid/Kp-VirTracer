#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(optparse)))

parser <- OptionParser(
  usage = "Rscript run_kpvirtracer.R --input <dir> --output <dir> [options]"
)
parser <- add_option(parser, "--input", type = "character", help = "Input FASTA directory")
parser <- add_option(parser, "--output", type = "character", help = "Output directory")
parser <- add_option(parser, "--threads", type = "integer", default = 4L, help = "Number of threads [default %default]")
parser <- add_option(parser, "--config", type = "character", default = NULL, help = "Optional config YAML path")
parser <- add_option(parser, "--virulence-db", type = "character", default = NULL, help = "BLAST DB prefix built by makeblastdb -out (required)")
parser <- add_option(parser, "--force", action = "store_true", default = FALSE, help = "Overwrite existing output")
parser <- add_option(parser, "--resume", action = "store_true", default = FALSE, help = "Resume from previous run")
parser <- add_option(parser, "--verbose", action = "store_true", default = FALSE, help = "Verbose mode")

opts <- parse_args(parser)
if (is.null(opts$input) || is.null(opts$output)) {
  print_help(parser)
  stop("Both --input and --output are required.", call. = FALSE)
}

if (!requireNamespace("KpVirTracer", quietly = TRUE)) {
  stop("Package 'KpVirTracer' is not installed. Run: Rscript install_package.R", call. = FALSE)
}

KpVirTracer::run_kp_virtracer(
  input = opts$input,
  output = opts$output,
  threads = opts$threads,
  config = opts$config,
  virulence_db = opts$`virulence-db`,
  force = opts$force,
  resume = opts$resume,
  verbose = opts$verbose
)

message("[INFO] KpVirTracer run completed.")
