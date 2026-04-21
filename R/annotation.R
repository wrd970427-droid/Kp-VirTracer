# Annotation module.

read_contig_report_with_location <- function(sample_id, mobsuite_dir) {
  if (is.null(mobsuite_dir) || !nzchar(mobsuite_dir)) {
    return(tibble::tibble())
  }
  f <- file.path(mobsuite_dir, sample_id, "contig_report.txt")
  if (!file.exists(f) || file.info(f)$size == 0) {
    return(tibble::tibble())
  }
  d <- tryCatch(readr::read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) {
    return(tibble::tibble())
  }
  nms <- names(d)
  id_col <- dplyr::coalesce(
    dplyr::if_else("contig_id" %in% nms, "contig_id", NA_character_),
    dplyr::if_else("sample_id" %in% nms, "sample_id", NA_character_),
    dplyr::if_else("id" %in% nms, "id", NA_character_)
  )
  if (is.na(id_col)) {
    return(tibble::tibble())
  }
  if (!("molecule_type" %in% nms)) {
    return(
      d %>%
        dplyr::transmute(
          Contig_ID = normalize_contig_name(.data[[id_col]]),
          molecule_type = "chromosome"
        )
    )
  }
  d %>%
    dplyr::transmute(
      Contig_ID = normalize_contig_name(.data[[id_col]]),
      molecule_type = tolower(as.character(molecule_type))
    ) %>%
    dplyr::mutate(
      molecule_type = dplyr::case_when(
        molecule_type == "plasmid" ~ "plasmid",
        molecule_type == "chromosome" ~ "chromosome",
        TRUE ~ "chromosome"
      )
    )
}

run_prodigal_annotation <- function(sample_manifest, virulence_hits, cfg, out_dir, tools, mobsuite_dir = NULL) {
  fs::dir_create(out_dir)
  ann_rows <- purrr::map_dfr(seq_len(nrow(sample_manifest)), function(i) {
    sid <- sample_manifest$Sample_ID[i]
    fasta <- sample_manifest$Input_FASTA[i]
    sample_dir <- file.path(out_dir, sid)
    fs::dir_create(sample_dir)

    contigs <- read_contig_report_with_location(sid, mobsuite_dir)
    sample_hits <- virulence_hits %>% dplyr::filter(Sample_ID == sid)
    vir_by_contig <- sample_hits %>%
      dplyr::group_by(Contig_ID) %>%
      dplyr::summarise(
        Virulence_Genes = paste(sort(unique(Gene)), collapse = ","),
        n_virulence_hits = dplyr::n(),
        .groups = "drop"
      )
    ann <- contigs %>%
      dplyr::left_join(vir_by_contig, by = "Contig_ID") %>%
      dplyr::mutate(
        Sample_ID = sid,
        n_virulence_hits = dplyr::coalesce(n_virulence_hits, 0L),
        Virulence_Genes = dplyr::coalesce(Virulence_Genes, "")
      ) %>%
      dplyr::select(Sample_ID, Contig_ID, molecule_type, n_virulence_hits, Virulence_Genes)

    gff_file <- file.path(sample_dir, paste0(sid, ".prodigal.gff"))
    faa_file <- file.path(sample_dir, paste0(sid, ".prodigal.faa"))
    fna_file <- file.path(sample_dir, paste0(sid, ".prodigal.fna"))
    prodigal_status <- "not_run"
    if (!is.null(tools[["prodigal"]]) && nzchar(tools[["prodigal"]])) {
      res <- tryCatch(
        run_external_tool(
          tools[["prodigal"]],
          args = c("-i", fasta, "-a", faa_file, "-d", fna_file, "-f", "gff", "-o", gff_file, "-p", "meta")
        ),
        error = function(e) NULL
      )
      if (!is.null(res) && isTRUE(res$status == 0)) {
        prodigal_status <- "ok"
        readr::write_lines(res$stdout, file.path(sample_dir, "prodigal.stdout.log"))
        readr::write_lines(res$stderr, file.path(sample_dir, "prodigal.stderr.log"))
      } else if (!is.null(res)) {
        prodigal_status <- "failed"
        readr::write_lines(res$stdout, file.path(sample_dir, "prodigal.stdout.log"))
        readr::write_lines(res$stderr, file.path(sample_dir, "prodigal.stderr.log"))
      }
    }

    if (nrow(ann) == 0) {
      ann <- tibble::tibble(
        Sample_ID = sid, Contig_ID = NA_character_, molecule_type = NA_character_,
        n_virulence_hits = 0L, Virulence_Genes = ""
      )
    }
    ann %>% dplyr::mutate(Prodigal_Status = prodigal_status)
  })

  readr::write_tsv(ann_rows, file.path(out_dir, "annotation_manifest.tsv"))
  summary_tbl <- ann_rows %>%
    dplyr::group_by(Sample_ID) %>%
    dplyr::summarise(
      n_contigs = dplyr::n_distinct(Contig_ID[!is.na(Contig_ID)]),
      n_virulence_contigs = sum(n_virulence_hits > 0, na.rm = TRUE),
      Prodigal_Status = dplyr::first(Prodigal_Status),
      .groups = "drop"
    )
  readr::write_tsv(summary_tbl, file.path(out_dir, "annotation_summary.tsv"))
  ann_rows
}
