#' Run Kp-VirTracer Pipeline
#'
#' Top-level pipeline runner. Stage-1 skeleton only.
#'
#' @param input Input FASTA directory.
#' @param output Output directory.
#' @param threads Number of threads.
#' @param config Optional path to YAML config.
#' @param virulence_db BLAST nucleotide database prefix for virulence genes.
#' @param force Overwrite outputs if existing.
#' @param resume Resume from previous outputs.
#' @param keep_temp Keep temporary files.
#' @param verbose Verbose logging mode.
#' @return Invisibly returns a list describing the run context.
#' @export
run_kp_virtracer <- function(input,
                             output,
                             threads = 4L,
                             config = NULL,
                             virulence_db = NULL,
                             force = FALSE,
                             resume = FALSE,
                             keep_temp = FALSE,
                             verbose = FALSE) {
  cfg <- load_kp_virtracer_config(config = config, threads = threads)
  if (is.null(cfg$blast)) {
    cfg$blast <- list()
  }
  if (!is.null(virulence_db) && nzchar(virulence_db)) {
    cfg$blast$virulence_db <- virulence_db
  }
  cfg$runtime$force <- force
  cfg$runtime$resume <- resume
  cfg$runtime$keep_temp <- keep_temp
  cfg$runtime$verbose <- verbose

  assert_dir_exists(input, "input directory")
  manifest <- discover_fasta_files(input)
  if (nrow(manifest) < 2) {
    stop("At least 2 FASTA samples are required.", call. = FALSE)
  }
  init_output_structure(output, force = force)
  write_sample_manifest(manifest, output)

  tools <- check_required_tools(cfg)
  klebo <- run_kleborate(manifest, cfg, file.path(output, "01_kleborate"), tools)
  mob <- run_mobsuite(manifest, cfg, file.path(output, "02_mobsuite"), tools)
  vir <- run_virulence_blast(manifest, cfg, file.path(output, "03_blast_virulence"), tools)
  ani <- run_pairwise_ani(manifest, cfg, file.path(output, "04_ani"), tools)
  hgt <- detect_hgt_events(NULL, vir, ani, cfg, file.path(output, "06_hgt"))

  sample_summary <- classify_hvkp_status(klebo, vir, cfg) %>%
    dplyr::left_join(
      mob %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::summarise(
          Plasmid_Count = sum(!is.na(plasmid_id)),
          Replicon_Types = paste(unique(stats::na.omit(replicon_type)), collapse = ","),
          PTU_Types = paste(unique(stats::na.omit(PTU)), collapse = ","),
          .groups = "drop"
        ),
      by = "Sample_ID"
    ) %>%
    dplyr::mutate(
      Plasmid_Count = dplyr::coalesce(Plasmid_Count, 0L),
      Replicon_Types = dplyr::coalesce(Replicon_Types, ""),
      PTU_Types = dplyr::coalesce(PTU_Types, "")
    ) %>%
    dplyr::select(
      Sample_ID, hvKp_Status, Virulence_Genes, Chromosomal_Genes, Plasmid_Genes,
      Virulence_Location, Plasmid_Count, Replicon_Types, PTU_Types, Kleborate_Virulence_Score
    )

  readr::write_tsv(ani, file.path(output, "04_ani", "ani_results.tsv"))
  readr::write_tsv(hgt, file.path(output, "06_hgt", "hgt_events.tsv"))
  write_pipeline_summaries(sample_summary, vir, ani, hgt, output)

  invisible(list(
    manifest = manifest,
    sample_summary = sample_summary,
    virulence_hits = vir,
    ani_results = ani,
    hgt_events = hgt
  ))
}
