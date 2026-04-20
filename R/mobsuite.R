# MOB-suite module.

get_first_present_col <- function(df, choices, fallback = NA_character_) {
  nm <- names(df)
  hit <- choices[choices %in% nm]
  if (length(hit) == 0) {
    return(rep(fallback, nrow(df)))
  }
  as.character(df[[hit[1]]])
}

read_table_if_exists <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(NULL)
  }
  tryCatch(readr::read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
}

read_mob_reports <- function(sample_dir) {
  mobtyper_candidates <- c(
    "mobtyper_results.txt",
    "mobtyper_results",
    "mobtyper_report.txt"
  )
  mobtyper <- NULL
  for (f in mobtyper_candidates) {
    mobtyper <- read_table_if_exists(file.path(sample_dir, f))
    if (!is.null(mobtyper)) {
      break
    }
  }
  contig <- read_table_if_exists(file.path(sample_dir, "contig_report.txt"))
  list(mobtyper = mobtyper, contig = contig)
}

build_mob_manifest_rows <- function(sid, reports, status) {
  mty <- reports$mobtyper
  ctr <- reports$contig

  # Preferred source: mobtyper_results(.txt), which is plasmid-level output.
  if (!is.null(mty) && nrow(mty) > 0) {
    return(
      tibble::tibble(
        Sample_ID = sid,
        MOB_Status = status,
        plasmid_id = get_first_present_col(mty, c("sample_id", "plasmid_id")),
        replicon_type = get_first_present_col(mty, c("rep_type(s)", "replicon_type", "rep_type")),
        relaxase = get_first_present_col(mty, c("relaxase_type(s)", "relaxase_type")),
        mobility = get_first_present_col(mty, c("predicted_mobility", "mobility")),
        predicted_host_range = get_first_present_col(
          mty,
          c("predicted_host_range_overall_rank", "predicted_host_range", "host_range")
        ),
        PTU = get_first_present_col(
          mty,
          c("primary_cluster_id", "secondary_cluster_id", "mash_nearest_neighbor")
        )
      )
    )
  }

  # Fallback source: contig_report.txt (filter to plasmid molecules when present).
  if (!is.null(ctr) && nrow(ctr) > 0) {
    plasmid_rows <- ctr
    if ("molecule_type" %in% names(plasmid_rows)) {
      plasmid_rows <- plasmid_rows %>%
        dplyr::mutate(molecule_type = tolower(as.character(molecule_type))) %>%
        dplyr::filter(molecule_type == "plasmid")
    }
    if (nrow(plasmid_rows) > 0) {
      return(
        tibble::tibble(
          Sample_ID = sid,
          MOB_Status = status,
          plasmid_id = get_first_present_col(plasmid_rows, c("sample_id", "plasmid_id", "contig_id")),
          replicon_type = get_first_present_col(plasmid_rows, c("rep_type(s)", "replicon_type", "rep_type")),
          relaxase = get_first_present_col(plasmid_rows, c("relaxase_type(s)", "relaxase_type")),
          mobility = get_first_present_col(plasmid_rows, c("predicted_mobility", "mobility")),
          predicted_host_range = get_first_present_col(
            plasmid_rows,
            c("predicted_host_range_overall_rank", "predicted_host_range", "mash_neighbor_identification")
          ),
          PTU = get_first_present_col(
            plasmid_rows,
            c("primary_cluster_id", "secondary_cluster_id", "mash_nearest_neighbor")
          )
        )
      )
    }
  }

  tibble::tibble(
    Sample_ID = sid, MOB_Status = status, plasmid_id = NA_character_,
    replicon_type = NA_character_, relaxase = NA_character_,
    mobility = NA_character_, predicted_host_range = NA_character_, PTU = NA_character_
  )
}

run_mobsuite <- function(manifest_df, cfg, out_dir, tools) {
  rows <- purrr::map_dfr(seq_len(nrow(manifest_df)), function(i) {
    sid <- manifest_df$Sample_ID[i]
    fasta <- manifest_df$Input_FASTA[i]
    sample_dir <- file.path(out_dir, sid)
    fs::dir_create(sample_dir)
    status <- "not_run"
    warn_reason <- NA_character_
    log_info(glue::glue("[MOB-suite] ({i}/{nrow(manifest_df)}) {sid}"))
    if (!is.null(tools[["mob_recon"]]) && nzchar(tools[["mob_recon"]])) {
      res <- tryCatch(
        run_external_tool(
          tools[["mob_recon"]],
          args = c("-i", fasta, "-o", sample_dir, "-n", as.character(cfg$threads))
        ),
        error = function(e) NULL
      )
      status <- if (!is.null(res) && isTRUE(res$status == 0)) "ok" else "failed"
      if (!is.null(res)) {
        readr::write_lines(res$stdout, file.path(sample_dir, "mob_recon.stdout.log"))
        readr::write_lines(res$stderr, file.path(sample_dir, "mob_recon.stderr.log"))
      }
      if (status == "failed") {
        warn_reason <- "mob_recon command failed; check mob_recon.stderr.log."
      }
    } else {
      warn_reason <- "mob_recon executable not available."
    }

    reports <- read_mob_reports(sample_dir)
    if (is.null(reports$mobtyper) && is.null(reports$contig) && is.na(warn_reason)) {
      warn_reason <- "No MOB-suite report found (expected mobtyper_results(.txt) or contig_report.txt)."
    }
    if (!is.na(warn_reason)) {
      log_warn(glue::glue("[MOB-suite] {sid}: {warn_reason}"))
    }
    build_mob_manifest_rows(sid, reports, status)
  })
  out_file <- file.path(dirname(out_dir), "plasmid_manifest.tsv")
  readr::write_tsv(rows, out_file)
  rows
}
