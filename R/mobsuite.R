# MOB-suite module.

read_mob_report <- function(sample_dir) {
  candidates <- c(
    "plasmid_report.txt",
    "mobtyper_results.txt",
    "mobtyper_report.txt",
    "contig_report.txt"
  )
  for (f in candidates) {
    p <- file.path(sample_dir, f)
    if (file.exists(p) && file.info(p)$size > 0) {
      df <- tryCatch(readr::read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(df)) {
        return(df)
      }
    }
  }
  tibble::tibble()
}

get_first_present_col <- function(df, choices, fallback = NA_character_) {
  nm <- names(df)
  hit <- choices[choices %in% nm]
  if (length(hit) == 0) {
    return(rep(fallback, nrow(df)))
  }
  as.character(df[[hit[1]]])
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
    rep <- read_mob_report(sample_dir)
    if (nrow(rep) == 0 && is.na(warn_reason)) {
      warn_reason <- "No MOB-suite report found (expected plasmid_report/mobtyper/contig report)."
    }
    if (!is.na(warn_reason)) {
      log_warn(glue::glue("[MOB-suite] {sid}: {warn_reason}"))
    }
    if (nrow(rep) == 0) {
      tibble::tibble(
        Sample_ID = sid, MOB_Status = status, plasmid_id = NA_character_,
        replicon_type = NA_character_, relaxase = NA_character_,
        mobility = NA_character_, predicted_host_range = NA_character_, PTU = NA_character_
      )
    } else {
      tibble::tibble(
        Sample_ID = sid,
        MOB_Status = status,
        plasmid_id = get_first_present_col(rep, c("plasmid_id", "sample_id", names(rep)[1])),
        replicon_type = get_first_present_col(rep, c("rep_type(s)", "replicon_type", "rep_type")),
        relaxase = get_first_present_col(rep, c("relaxase_type(s)", "relaxase_type")),
        mobility = get_first_present_col(rep, c("predicted_mobility", "mobility")),
        predicted_host_range = get_first_present_col(rep, c("predicted_host_range_overall_rank", "predicted_host_range")),
        PTU = get_first_present_col(rep, c("mash_nearest_neighbor", "ptu", "PTU"))
      )
    }
  })
  out_file <- file.path(dirname(out_dir), "plasmid_manifest.tsv")
  readr::write_tsv(rows, out_file)
  rows
}
