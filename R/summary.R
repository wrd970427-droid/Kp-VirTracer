# Summary output module.

write_pipeline_summaries <- function(sample_summary, virulence_hits, ani_results, hgt_events, out_dir) {
  summary_dir <- file.path(out_dir, "summary")
  fs::dir_create(summary_dir)
  readr::write_tsv(sample_summary, file.path(summary_dir, "sample_summary.tsv"))
  virulence_type_summary <- sample_summary %>%
    dplyr::count(Virulence_Type, name = "n_samples") %>%
    dplyr::arrange(dplyr::desc(n_samples), Virulence_Type)
  readr::write_tsv(virulence_type_summary, file.path(summary_dir, "virulence_type_summary.tsv"))
  readr::write_tsv(virulence_hits, file.path(summary_dir, "virulence_hits.tsv"))
  final_summary <- tibble::tibble(
    n_samples = dplyr::n_distinct(sample_summary$Sample_ID),
    n_virulence = sum(sample_summary$Virulence_Type != "nKp", na.rm = TRUE),
    n_ani_pairs = nrow(ani_results),
    n_hgt_events = nrow(hgt_events)
  )
  readr::write_tsv(final_summary, file.path(summary_dir, "final_summary.tsv"))
  invisible(final_summary)
}
