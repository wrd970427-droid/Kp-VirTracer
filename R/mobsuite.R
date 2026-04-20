# MOB-suite module.

run_mobsuite <- function(manifest_df, cfg, out_dir, tools) {
  rows <- purrr::map_dfr(seq_len(nrow(manifest_df)), function(i) {
    sid <- manifest_df$Sample_ID[i]
    fasta <- manifest_df$Input_FASTA[i]
    sample_dir <- file.path(out_dir, sid)
    fs::dir_create(sample_dir)
    status <- "not_run"
    if (!is.null(tools[["mob_recon"]]) && nzchar(tools[["mob_recon"]])) {
      res <- tryCatch(
        run_external_tool(
          tools[["mob_recon"]],
          args = c("-i", fasta, "-o", sample_dir, "-n", as.character(cfg$threads), "-s", sid)
        ),
        error = function(e) NULL
      )
      status <- if (!is.null(res) && isTRUE(res$status == 0)) "ok" else "failed"
    }
    rep_file <- file.path(sample_dir, "plasmid_report.txt")
    rep <- if (file.exists(rep_file)) {
      tryCatch(readr::read_tsv(rep_file, show_col_types = FALSE), error = function(e) tibble::tibble())
    } else tibble::tibble()
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
        plasmid_id = as.character(rep[[1]]),
        replicon_type = if ("rep_type(s)" %in% names(rep)) as.character(rep[["rep_type(s)"]]) else NA_character_,
        relaxase = if ("relaxase_type(s)" %in% names(rep)) as.character(rep[["relaxase_type(s)"]]) else NA_character_,
        mobility = if ("predicted_mobility" %in% names(rep)) as.character(rep[["predicted_mobility"]]) else NA_character_,
        predicted_host_range = if ("predicted_host_range_overall_rank" %in% names(rep)) as.character(rep[["predicted_host_range_overall_rank"]]) else NA_character_,
        PTU = if ("mash_nearest_neighbor" %in% names(rep)) as.character(rep[["mash_nearest_neighbor"]]) else NA_character_
      )
    }
  })
  out_file <- file.path(dirname(out_dir), "plasmid_manifest.tsv")
  readr::write_tsv(rows, out_file)
  rows
}
