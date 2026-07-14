#' Run Kp-VirTracer Pipeline
#'
#' Top-level pipeline runner. Stage-1 skeleton only.
#'
#' @param input Input FASTA directory.
#' @param output Output directory.
#' @param threads Number of threads.
#' @param config Optional path to YAML config.
#' @param virulence_db BLAST nucleotide database prefix for virulence genes.
#' @param blast_task BLAST task override, e.g. `blastn` or `blastn-short`.
#' @param blast_evalue BLAST e-value cutoff override.
#' @param blast_min_identity BLAST minimum identity (%) override.
#' @param blast_min_coverage BLAST minimum coverage (%) override.
#' @param ani_relatedness ANI relatedness cutoff (%) override.
#' @param ani_min_fraction ANI minimum mapped fraction override.
#' @param ani_frag_len ANI fragment length override.
#' @param ani_kmer ANI k-mer size override.
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
                             blast_task = NULL,
                             blast_evalue = NULL,
                             blast_min_identity = NULL,
                             blast_min_coverage = NULL,
                             ani_relatedness = NULL,
                             ani_min_fraction = NULL,
                             ani_frag_len = NULL,
                             ani_kmer = NULL,
                             force = FALSE,
                             resume = FALSE,
                             keep_temp = FALSE,
                             verbose = FALSE) {
  log_info("Kp-VirTracer pipeline started.")
  cfg <- load_kp_virtracer_config(config = config, threads = threads)
  if (is.null(cfg$blast)) {
    cfg$blast <- list()
  }
  if (!is.null(virulence_db) && nzchar(virulence_db)) {
    cfg$blast$virulence_db <- virulence_db
  }
  if (!is.null(blast_task) && nzchar(blast_task)) {
    cfg$blast$task <- blast_task
  }
  if (!is.null(blast_evalue)) {
    cfg$blast$evalue <- as.numeric(blast_evalue)
  }
  if (!is.null(blast_min_identity)) {
    cfg$blast$min_identity <- as.numeric(blast_min_identity)
  }
  if (!is.null(blast_min_coverage)) {
    cfg$blast$min_coverage <- as.numeric(blast_min_coverage)
  }
  if (is.null(cfg$ani)) {
    cfg$ani <- list()
  }
  if (!is.null(ani_relatedness)) {
    cfg$ani$relatedness_threshold <- as.numeric(ani_relatedness)
    cfg$ani$related_cutoff <- as.numeric(ani_relatedness)
  }
  if (!is.null(ani_min_fraction)) {
    cfg$ani$min_fraction <- as.numeric(ani_min_fraction)
  }
  if (!is.null(ani_frag_len)) {
    cfg$ani$frag_len <- as.integer(ani_frag_len)
  }
  if (!is.null(ani_kmer)) {
    cfg$ani$kmer <- as.integer(ani_kmer)
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
  log_info(glue::glue("Discovered {nrow(manifest)} samples."))
  init_output_structure(output, force = force)
  write_sample_manifest(manifest, output)
  log_info(glue::glue("Output initialized at: {output}"))

  tools <- check_required_tools(cfg)
  log_info("Stage 1/6: running Kleborate.")
  klebo <- run_kleborate(manifest, cfg, file.path(output, "01_kleborate"), tools)
  log_info("Stage 2/6: running MOB-suite.")
  mob <- run_mobsuite(manifest, cfg, file.path(output, "02_mobsuite"), tools)
  log_info("Stage 3/6: running virulence BLAST.")
  vir <- run_virulence_blast(
    manifest, cfg, file.path(output, "03_blast_virulence"), tools,
    mobsuite_dir = file.path(output, "02_mobsuite")
  )

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
      Sample_ID, Virulence_Type, hvKp_Status, Virulence_Genes, Chromosomal_Genes, Plasmid_Genes,
      Virulence_Location, Plasmid_Count, Replicon_Types, PTU_Types, Kleborate_Virulence_Score
    )

  log_info("Stage 4/6: running ANI on virulent-strain chromosomes only.")
  virulent_ids <- sample_summary %>%
    dplyr::filter(Virulence_Type != "nKp") %>%
    dplyr::pull(Sample_ID)
  vir_manifest <- manifest %>%
    dplyr::filter(Sample_ID %in% virulent_ids) %>%
    dplyr::select(Sample_ID, Input_FASTA)
  chr_manifest <- build_chromosome_manifest(
    vir_manifest,
    mobsuite_dir = file.path(output, "02_mobsuite"),
    out_dir = file.path(output, "04_ani")
  )
  ani <- run_pairwise_ani(chr_manifest, cfg, file.path(output, "04_ani"), tools)
  log_info("Stage 5/6: running annotation.")
  ann <- run_prodigal_annotation(
    sample_manifest = manifest,
    virulence_hits = vir,
    cfg = cfg,
    out_dir = file.path(output, "05_annotation"),
    tools = tools,
    mobsuite_dir = file.path(output, "02_mobsuite")
  )
  log_info("Stage 6/6: inferring HGT events and writing summaries.")
  hgt <- detect_hgt_events(
    NULL, vir, ani, cfg, file.path(output, "06_hgt"),
    annotation_df = ann,
    annotation_dir = file.path(output, "05_annotation"),
    tools = tools
  )

  readr::write_tsv(ani, file.path(output, "04_ani", "ani_results.tsv"))
  readr::write_tsv(hgt, file.path(output, "06_hgt", "hgt_events.tsv"))
  write_pipeline_summaries(sample_summary, vir, ani, hgt, output)
  log_info("Kp-VirTracer pipeline completed.")

  invisible(list(
    manifest = manifest,
    sample_summary = sample_summary,
    virulence_hits = vir,
    ani_results = ani,
    annotation = ann,
    hgt_events = hgt
  ))
}
