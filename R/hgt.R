# HGT module (minimal scaffold with stable schema).

empty_hgt_events <- function() {
  tibble::tibble(
    Sample_A = character(), Sample_B = character(), Gene = character(),
    Location_A = character(), Location_B = character(),
    Plasmid_A = character(), Plasmid_B = character(),
    HGT_Type = character(), Confirmed = logical(),
    pident = numeric(), qcov = numeric(), Synteny_Index = integer(),
    Evidence_Files = character()
  )
}

detect_hgt_events <- function(sample_pairs, virulence_hits, ani_results, cfg, out_dir) {
  events <- empty_hgt_events()
  readr::write_tsv(events, file.path(out_dir, "hgt_events.tsv"))
  readr::write_tsv(tibble::tibble(), file.path(out_dir, "highly_similar_plasmids.tsv"))
  readr::write_tsv(tibble::tibble(), file.path(out_dir, "synteny_results.tsv"))
  events
}
