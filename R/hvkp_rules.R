# hvKp/nKp rule engine.

classify_hvkp_status <- function(kleborate_summary, virulence_hits, cfg) {
  sample_ids <- unique(kleborate_summary$Sample_ID)
  hit_summary <- virulence_hits %>%
    dplyr::group_by(Sample_ID) %>%
    dplyr::summarise(
      Virulence_Genes = paste(sort(unique(Gene)), collapse = ","),
      Chromosomal_Genes = paste(sort(unique(Gene[Location == "chromosome"])), collapse = ","),
      Plasmid_Genes = paste(sort(unique(Gene[Location == "plasmid"])), collapse = ","),
      .groups = "drop"
    )
  out <- tibble::tibble(Sample_ID = sample_ids) %>%
    dplyr::left_join(kleborate_summary %>% dplyr::select(Sample_ID, Kleborate_Virulence_Score), by = "Sample_ID") %>%
    dplyr::left_join(hit_summary, by = "Sample_ID") %>%
    dplyr::mutate(
      Virulence_Genes = dplyr::coalesce(Virulence_Genes, ""),
      Chromosomal_Genes = dplyr::coalesce(Chromosomal_Genes, ""),
      Plasmid_Genes = dplyr::coalesce(Plasmid_Genes, ""),
      Virulence_Location = dplyr::case_when(
        Chromosomal_Genes != "" & Plasmid_Genes != "" ~ "both",
        Chromosomal_Genes != "" ~ "chromosome",
        Plasmid_Genes != "" ~ "plasmid",
        TRUE ~ "absent"
      ),
      Virulence_Type = dplyr::case_when(
        Virulence_Location == "both" ~ "pc-hvKp",
        Virulence_Location == "chromosome" ~ "c-hvKp",
        Virulence_Location == "plasmid" ~ "p-hvKp",
        TRUE ~ "nKp"
      ),
      hvKp_Status = dplyr::if_else(Virulence_Type == "nKp", "nKp", "hvKp")
    ) %>%
    dplyr::select(
      Sample_ID,
      Virulence_Type,
      hvKp_Status,
      Virulence_Location,
      Virulence_Genes,
      Chromosomal_Genes,
      Plasmid_Genes,
      Kleborate_Virulence_Score
    ) %>%
    dplyr::arrange(
      dplyr::desc(Virulence_Type != "nKp"),
      Sample_ID
    )
  out
}
