# HGT module.

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

detect_hgt_events <- function(sample_pairs, virulence_hits, ani_results, cfg, out_dir, annotation_df = NULL) {
  fs::dir_create(out_dir)
  if (nrow(virulence_hits) == 0) {
    events <- empty_hgt_events()
    readr::write_tsv(events, file.path(out_dir, "hgt_events.tsv"))
    readr::write_tsv(tibble::tibble(), file.path(out_dir, "highly_similar_plasmids.tsv"))
    readr::write_tsv(tibble::tibble(), file.path(out_dir, "synteny_results.tsv"))
    return(events)
  }

  by_sample_gene <- virulence_hits %>%
    dplyr::group_by(Sample_ID, Gene) %>%
    dplyr::summarise(
      Location = dplyr::case_when(
        any(Location == "chromosome") & any(Location == "plasmid") ~ "both",
        any(Location == "plasmid") ~ "plasmid",
        TRUE ~ "chromosome"
      ),
      Mean_Identity = mean(Identity, na.rm = TRUE),
      Mean_Coverage = mean(Coverage, na.rm = TRUE),
      .groups = "drop"
    )

  related_pairs <- ani_results %>%
    dplyr::transmute(
      Sample_A, Sample_B,
      ANI = ANI,
      Relatedness = Relatedness,
      Pair_Related = Relatedness == "Related"
    )
  if (nrow(related_pairs) == 0) {
    sids <- unique(by_sample_gene$Sample_ID)
    if (length(sids) < 2) {
      events <- empty_hgt_events()
      readr::write_tsv(events, file.path(out_dir, "hgt_events.tsv"))
      readr::write_tsv(tibble::tibble(), file.path(out_dir, "highly_similar_plasmids.tsv"))
      readr::write_tsv(tibble::tibble(), file.path(out_dir, "synteny_results.tsv"))
      return(events)
    }
    pairs <- t(utils::combn(sids, 2))
    related_pairs <- tibble::tibble(
      Sample_A = pairs[, 1],
      Sample_B = pairs[, 2],
      ANI = NA_real_,
      Relatedness = "Unknown",
      Pair_Related = FALSE
    )
  }

  events <- purrr::map_dfr(seq_len(nrow(related_pairs)), function(i) {
    a <- related_pairs$Sample_A[i]
    b <- related_pairs$Sample_B[i]
    ag <- by_sample_gene %>% dplyr::filter(Sample_ID == a)
    bg <- by_sample_gene %>% dplyr::filter(Sample_ID == b)
    shared <- intersect(ag$Gene, bg$Gene)
    if (length(shared) == 0) {
      return(tibble::tibble())
    }
    purrr::map_dfr(shared, function(g) {
      ga <- ag %>% dplyr::filter(Gene == g) %>% dplyr::slice(1)
      gb <- bg %>% dplyr::filter(Gene == g) %>% dplyr::slice(1)
      loc_a <- ga$Location[1]
      loc_b <- gb$Location[1]
      pair_related <- isTRUE(related_pairs$Pair_Related[i])
      hgt_type <- dplyr::case_when(
        loc_a != loc_b ~ "location_switch",
        loc_a == "plasmid" & loc_b == "plasmid" ~ "shared_plasmid_gene",
        TRUE ~ "shared_chromosomal_gene"
      )
      confirmed <- pair_related || hgt_type == "location_switch"
      tibble::tibble(
        Sample_A = a,
        Sample_B = b,
        Gene = g,
        Location_A = loc_a,
        Location_B = loc_b,
        Plasmid_A = ifelse(loc_a == "plasmid", "plasmid", NA_character_),
        Plasmid_B = ifelse(loc_b == "plasmid", "plasmid", NA_character_),
        HGT_Type = hgt_type,
        Confirmed = confirmed,
        pident = mean(c(ga$Mean_Identity[1], gb$Mean_Identity[1]), na.rm = TRUE),
        qcov = mean(c(ga$Mean_Coverage[1], gb$Mean_Coverage[1]), na.rm = TRUE),
        Synteny_Index = ifelse(is.null(annotation_df), NA_integer_, 1L),
        Evidence_Files = "03_blast_virulence/summary/virulence_hits.tsv;04_ani/ani_results.tsv"
      )
    })
  })

  if (nrow(events) == 0) {
    events <- empty_hgt_events()
  }
  highly_similar_plasmids <- events %>%
    dplyr::filter(Location_A == "plasmid" | Location_B == "plasmid") %>%
    dplyr::transmute(Sample_A, Sample_B, Gene, pident, qcov, HGT_Type)
  synteny_results <- events %>%
    dplyr::transmute(Sample_A, Sample_B, Gene, Synteny_Index, Confirmed)

  readr::write_tsv(events, file.path(out_dir, "hgt_events.tsv"))
  readr::write_tsv(highly_similar_plasmids, file.path(out_dir, "highly_similar_plasmids.tsv"))
  readr::write_tsv(synteny_results, file.path(out_dir, "synteny_results.tsv"))
  events
}
