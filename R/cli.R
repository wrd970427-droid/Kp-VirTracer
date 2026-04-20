#' KpVirTracer CLI Main Entry
#'
#' Parse command-line arguments and dispatch the pipeline runner.
#'
#' @param args Character vector of CLI arguments.
#' @return Invisibly returns `NULL`.
#' @export
kp_virtracer_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  parser <- optparse::OptionParser(
    usage = "KpVirTracer --input <dir> --output <dir> [options]"
  )
  parser <- optparse::add_option(parser, "--input", type = "character", help = "Input FASTA directory")
  parser <- optparse::add_option(parser, "--output", type = "character", help = "Output directory")
  parser <- optparse::add_option(parser, "--threads", type = "integer", default = 4L)
  parser <- optparse::add_option(parser, "--config", type = "character", default = NULL)
  parser <- optparse::add_option(
    parser, "--virulence-db", type = "character", default = NULL,
    help = "BLAST DB prefix built by makeblastdb -out (required)"
  )
  parser <- optparse::add_option(parser, "--force", action = "store_true", default = FALSE)
  parser <- optparse::add_option(parser, "--resume", action = "store_true", default = FALSE)
  parser <- optparse::add_option(parser, "--keep-temp", action = "store_true", default = FALSE)
  parser <- optparse::add_option(parser, "--verbose", action = "store_true", default = FALSE)

  opts <- optparse::parse_args(parser, args = args)
  if (is.null(opts$input) || is.null(opts$output)) {
    optparse::print_help(parser)
    stop("Both --input and --output are required.", call. = FALSE)
  }
  run_kp_virtracer(
    input = opts$input,
    output = opts$output,
    threads = opts$threads,
    config = opts$config,
    virulence_db = opts$`virulence-db`,
    force = opts$force,
    resume = opts$resume,
    keep_temp = opts$`keep-temp`,
    verbose = opts$verbose
  )
  invisible(NULL)
}
