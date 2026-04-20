# Kleborate module.

run_kleborate <- function(manifest_df, cfg, out_dir, tools) {
  raw_dir <- file.path(out_dir, "kleborate_raw")
  fs::dir_create(raw_dir)
  out_rows <- purrr::map_dfr(seq_len(nrow(manifest_df)), function(i) {
    sid <- manifest_df$Sample_ID[i]
    fasta <- manifest_df$Input_FASTA[i]
    sample_dir <- file.path(raw_dir, sid)
    fs::dir_create(sample_dir)
    status <- "not_run"
    species <- NA_character_
    lineage <- NA_character_
    virulence_score <- NA_real_
    aerobactin <- NA_character_
    salmochelin <- NA_character_
    rmpA <- NA_character_
    rmpA2 <- NA_character_

    if (!is.null(tools[["kleborate"]]) && nzchar(tools[["kleborate"]])) {
      res <- tryCatch(
        run_external_tool(
          tools[["kleborate"]],
          args = c("-a", fasta, "-o", sample_dir, "-p", "kpsc")
        ),
        error = function(e) NULL
      )
      if (!is.null(res) && isTRUE(res$status == 0)) {
        status <- "ok"
        txt <- list.files(sample_dir, pattern = "\\.txt$", full.names = TRUE)
        if (length(txt) > 0) {
          df <- tryCatch(readr::read_tsv(txt[1], show_col_types = FALSE), error = function(e) NULL)
          if (!is.null(df) && nrow(df) > 0) {
            nm <- names(df)
            get_col <- function(x) if (x %in% nm) as.character(df[[x]][1]) else NA_character_
            species <- dplyr::coalesce(get_col("species"), get_col("Species"))
            lineage <- dplyr::coalesce(get_col("st"), get_col("ST"), get_col("lineage"))
            virulence_score <- suppressWarnings(as.numeric(dplyr::coalesce(get_col("virulence_score"), get_col("virulence_score/kleb"), NA_character_)))
            aerobactin <- dplyr::coalesce(get_col("iuc"), get_col("aerobactin"))
            salmochelin <- dplyr::coalesce(get_col("iro"), get_col("salmochelin"))
            rmpA <- get_col("rmpA")
            rmpA2 <- get_col("rmpA2")
          }
        }
      } else {
        status <- "failed"
      }
    }
    tibble::tibble(
      Sample_ID = sid,
      Kleborate_Status = status,
      Species = species,
      ST_or_Lineage = lineage,
      Kleborate_Virulence_Score = virulence_score,
      Aerobactin = aerobactin,
      Salmochelin = salmochelin,
      rmpA = rmpA,
      rmpA2 = rmpA2
    )
  })
  readr::write_tsv(out_rows, file.path(out_dir, "kleborate_summary.tsv"))
  out_rows
}
