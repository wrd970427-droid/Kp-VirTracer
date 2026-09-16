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

detect_hgt_events <- function(sample_pairs,
                              virulence_hits,
                              ani_results,
                              cfg,
                              out_dir,
                              annotation_df = NULL,
                              annotation_dir = NULL,
                              tools = NULL) {
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

  # Resolve annotation_dir from attribute if not passed explicitly.
  if (is.null(annotation_dir) && !is.null(annotation_df)) {
    annotation_dir <- attr(annotation_df, "annotation_dir", exact = TRUE)
  }

  syn_cfg <- cfg$synteny
  synteny_cutoff <- as.integer(dplyr::coalesce(syn_cfg$vertical_cutoff, syn_cfg$cutoff, 7L))
  annotation_ctx <- NULL
  if (!is.null(annotation_dir) && nzchar(annotation_dir) && dir.exists(annotation_dir)) {
    blastp_bin <- get_tool(tools, "blastp")
    if (is.null(blastp_bin)) {
      blastn_bin <- get_tool(tools, "blastn")
      if (!is.null(blastn_bin)) {
        # blastp usually sits next to blastn
        cand <- file.path(dirname(blastn_bin), "blastp")
        if (file.exists(cand)) {
          blastp_bin <- cand
        }
      }
    }
    if (is.null(blastp_bin)) {
      wh <- Sys.which("blastp")
      if (nzchar(wh)) {
        blastp_bin <- wh
      }
    }
    annotation_ctx <- list(
      annotation_dir = annotation_dir,
      virulence_hits = virulence_hits,
      blastp = blastp_bin
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

      pair_record <- list(
        Sample_A = a,
        Sample_B = b,
        Gene = g,
        Location_A = loc_a,
        Location_B = loc_b,
        Pair_Related = pair_related
      )
      # Synteny only for related pairs + chromosomal virulence genes.
      si <- if (is.null(annotation_ctx) || !pair_related) {
        NA_integer_
      } else if (!(loc_a %in% c("chromosome", "both") && loc_b %in% c("chromosome", "both"))) {
        NA_integer_
      } else {
        tryCatch(
          compute_synteny_index(pair_record, annotation_ctx, cfg),
          error = function(e) {
            log_warn(glue::glue("[Synteny] {a} vs {b} / {g}: {conditionMessage(e)}"))
            NA_integer_
          }
        )
      }

      # Chromosomal synteny (when computed) gates confirmation; plasmid events
      # still confirm by relatedness / location_switch without synteny.
      base_ok <- pair_related || hgt_type == "location_switch"
      if (!is.na(si)) {
        confirmed <- base_ok && si >= synteny_cutoff
      } else {
        confirmed <- base_ok
      }

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
        Synteny_Index = si,
        Evidence_Files = "03_blast_virulence/summary/virulence_hits.tsv;04_ani/ani_results.tsv;05_annotation/*.prodigal.gff"
      )
    })
  })

  if (nrow(events) == 0) {
    events <- empty_hgt_events()
  }
  highly_similar_plasmids <- events %>%
    dplyr::filter(Location_A == "plasmid" | Location_B == "plasmid") %>%
    dplyr::transmute(Sample_A, Sample_B, Gene, pident, qcov, HGT_Type, Synteny_Index)
  synteny_results <- events %>%
    dplyr::transmute(Sample_A, Sample_B, Gene, Synteny_Index, Confirmed, HGT_Type)

  readr::write_tsv(events, file.path(out_dir, "hgt_events.tsv"))
  readr::write_tsv(highly_similar_plasmids, file.path(out_dir, "highly_similar_plasmids.tsv"))
  readr::write_tsv(synteny_results, file.path(out_dir, "synteny_results.tsv"))
  events
}
