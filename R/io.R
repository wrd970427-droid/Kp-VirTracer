#' Discover FASTA Files
#'
#' @param input_dir Directory with `.fa`, `.fasta`, `.fna` files.
#' @return Tibble manifest skeleton.
discover_fasta_files <- function(input_dir) {
  assert_dir_exists(input_dir, "input directory")
  pattern <- "\\.(fa|fasta|fna)$"
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop("No FASTA files found in input directory.", call. = FALSE)
  }
  sample_ids <- fs::path_ext_remove(basename(files))
  if (any(duplicated(sample_ids))) {
    dups <- unique(sample_ids[duplicated(sample_ids)])
    stop(glue::glue("Duplicate sample IDs detected: {paste(dups, collapse = ', ')}"), call. = FALSE)
  }
  tibble::tibble(
    Sample_ID = sample_ids,
    Input_FASTA = files,
    File_Size = file.info(files)$size,
    Status = "pending"
  )
}

init_output_structure <- function(output_dir, force = FALSE) {
  if (dir.exists(output_dir) && !force) {
    # keep existing for resume mode
  } else if (dir.exists(output_dir) && force) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  fs::dir_create(output_dir)
  subdirs <- c(
    "logs", "temp", "01_kleborate", "01_kleborate/kleborate_raw",
    "02_mobsuite", "03_blast_virulence", "04_ani", "05_annotation",
    "06_hgt", "summary"
  )
  purrr::walk(file.path(output_dir, subdirs), fs::dir_create)
  invisible(output_dir)
}

write_sample_manifest <- function(manifest_df, output_dir) {
  out <- file.path(output_dir, "sample_manifest.tsv")
  readr::write_tsv(manifest_df, out)
  out
}
