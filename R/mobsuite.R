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

normalize_record_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  stringr::str_split_fixed(x, "\\s+", 2)[, 1]
}

read_fasta_records_mob <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(list(headers = character(), seqs = character()))
  }
  ln <- readLines(path, warn = FALSE)
  hdr_idx <- grep("^>", ln)
  if (length(hdr_idx) == 0) {
    return(list(headers = character(), seqs = character()))
  }
  starts <- hdr_idx
  ends <- c(hdr_idx[-1] - 1L, length(ln))
  headers <- sub("^>", "", ln[hdr_idx])
  seqs <- vapply(
    seq_along(starts),
    function(i) paste0(ln[(starts[i] + 1L):ends[i]], collapse = ""),
    character(1)
  )
  list(headers = headers, seqs = seqs)
}

write_fasta_records_mob <- function(headers, seqs, path) {
  if (length(headers) == 0 || length(seqs) == 0) {
    return(FALSE)
  }
  lines <- unlist(purrr::map2(headers, seqs, ~ c(paste0(">", .x), .y)), use.names = FALSE)
  readr::write_lines(lines, path)
  TRUE
}

read_mob_reports <- function(sample_dir) {
  mobtyper_candidates <- c(
    "mobtyper_results.txt",
    "mobtyper_results",
    "mobtyper_results.tsv",
    "mobtyper_report.txt",
    "mobtyper_report.tsv"
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
        PTU = NA_character_,
        PTU_Score = NA_real_,
        PTU_Host_Range = NA_character_,
        PTU_Notes = NA_character_
      )
    )
  }

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
          PTU = NA_character_,
          PTU_Score = NA_real_,
          PTU_Host_Range = NA_character_,
          PTU_Notes = NA_character_
        )
      )
    }
  }

  tibble::tibble(
    Sample_ID = sid, MOB_Status = status, plasmid_id = NA_character_,
    replicon_type = NA_character_, relaxase = NA_character_,
    mobility = NA_character_, predicted_host_range = NA_character_,
    PTU = NA_character_, PTU_Score = NA_real_, PTU_Host_Range = NA_character_, PTU_Notes = NA_character_
  )
}

extract_plasmid_fasta <- function(sample_dir, sample_fasta) {
  contig_report <- file.path(sample_dir, "contig_report.txt")
  if (!file.exists(contig_report) || file.info(contig_report)$size == 0) {
    return(NULL)
  }
  ctr <- tryCatch(readr::read_tsv(contig_report, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(ctr) || nrow(ctr) == 0 || !("molecule_type" %in% names(ctr))) {
    return(NULL)
  }
  id_col <- dplyr::coalesce(
    dplyr::if_else("contig_id" %in% names(ctr), "contig_id", NA_character_),
    dplyr::if_else("sample_id" %in% names(ctr), "sample_id", NA_character_)
  )
  if (is.na(id_col)) {
    return(NULL)
  }
  plasmid_ids <- ctr %>%
    dplyr::mutate(molecule_type = tolower(as.character(molecule_type))) %>%
    dplyr::filter(molecule_type == "plasmid") %>%
    dplyr::pull(.data[[id_col]]) %>%
    normalize_record_id() %>%
    unique()
  if (length(plasmid_ids) == 0) {
    return(NULL)
  }
  rec <- read_fasta_records_mob(sample_fasta)
  if (length(rec$headers) == 0) {
    return(NULL)
  }
  keep <- normalize_record_id(rec$headers) %in% plasmid_ids
  if (!any(keep)) {
    return(NULL)
  }
  out <- file.path(sample_dir, "copla_input_plasmids.fasta")
  ok <- write_fasta_records_mob(rec$headers[keep], rec$seqs[keep], out)
  if (!ok) {
    return(NULL)
  }
  out
}

resolve_copla_runtime <- function(cfg) {
  c0 <- cfg$copla
  if (is.null(c0)) {
    return(list(ready = FALSE, reason = "copla section missing in config."))
  }
  enabled <- isTRUE(c0$enabled)
  if (!enabled) {
    return(list(ready = FALSE, reason = "copla.enabled is false."))
  }

  script_path <- c0$script_path
  pickle_path <- c0$pickle_path
  fofn_path <- c0$fofn_path
  conda_env <- dplyr::coalesce(c0$conda_env, "copla")
  conda_bin <- dplyr::coalesce(c0$conda_bin, "conda")
  python_bin <- dplyr::coalesce(c0$python_bin, "python3")
  topology <- dplyr::coalesce(c0$topology, "linear")
  project_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = FALSE)

  if (is.null(script_path) || !nzchar(script_path) || !file.exists(script_path)) {
    return(list(ready = FALSE, reason = "copla script_path is not configured or not found."))
  }
  if (!dir.exists(project_root)) {
    return(list(ready = FALSE, reason = "copla project root cannot be resolved from script_path."))
  }
  if (!file.exists(file.path(project_root, "bin", "get_ani_identity.pl"))) {
    return(list(ready = FALSE, reason = "COPLA helper script bin/get_ani_identity.pl not found under project root."))
  }
  if (is.null(pickle_path) || !nzchar(pickle_path) || !file.exists(pickle_path)) {
    return(list(ready = FALSE, reason = "copla pickle_path is not configured or not found."))
  }
  if (is.null(fofn_path) || !nzchar(fofn_path) || !file.exists(fofn_path)) {
    return(list(ready = FALSE, reason = "copla fofn_path is not configured or not found."))
  }

  conda_resolved <- Sys.which(conda_bin)
  if (is.null(conda_resolved) || !nzchar(conda_resolved)) {
    return(list(ready = FALSE, reason = paste0("conda binary not found: ", conda_bin)))
  }

  list(
    ready = TRUE,
    conda = conda_resolved,
    conda_env = conda_env,
    python_bin = python_bin,
    script_path = script_path,
    project_root = project_root,
    pickle_path = pickle_path,
    fofn_path = fofn_path,
    topology = topology
  )
}

find_copla_prediction_file <- function(out_dir) {
  if (!dir.exists(out_dir)) {
    return(NA_character_)
  }
  files <- list.files(out_dir, pattern = "ptu_prediction\\.tsv$", full.names = TRUE)
  if (length(files) == 0) {
    return(NA_character_)
  }
  files[1]
}

read_copla_prediction <- function(path, sid) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(
      tibble::tibble(
        Sample_ID = character(),
        plasmid_id = character(),
        PTU = character(),
        PTU_Score = numeric(),
        PTU_Host_Range = character(),
        PTU_Notes = character()
      )
    )
  }
  d <- tryCatch(readr::read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) {
    return(
      tibble::tibble(
        Sample_ID = character(),
        plasmid_id = character(),
        PTU = character(),
        PTU_Score = numeric(),
        PTU_Host_Range = character(),
        PTU_Notes = character()
      )
    )
  }
  nms <- names(d)
  ptu_col <- dplyr::coalesce(
    dplyr::if_else(any(grepl("^#?predicted$", nms, ignore.case = TRUE)), nms[grepl("^#?predicted$", nms, ignore.case = TRUE)][1], NA_character_),
    dplyr::if_else(any(grepl("^ptu$", nms, ignore.case = TRUE)), nms[grepl("^ptu$", nms, ignore.case = TRUE)][1], NA_character_),
    dplyr::if_else(any(grepl("ptu", nms, ignore.case = TRUE)), nms[grepl("ptu", nms, ignore.case = TRUE)][1], NA_character_)
  )
  score_col <- dplyr::if_else(any(grepl("score", nms, ignore.case = TRUE)), nms[grepl("score", nms, ignore.case = TRUE)][1], NA_character_)
  host_range_col <- dplyr::if_else(any(grepl("host[_ ]?range", nms, ignore.case = TRUE)), nms[grepl("host[_ ]?range", nms, ignore.case = TRUE)][1], NA_character_)
  notes_col <- dplyr::if_else(any(grepl("notes?", nms, ignore.case = TRUE)), nms[grepl("notes?", nms, ignore.case = TRUE)][1], NA_character_)
  id_col <- dplyr::if_else(any(grepl("query|plasmid|seq|accession", nms, ignore.case = TRUE)), nms[grepl("query|plasmid|seq|accession", nms, ignore.case = TRUE)][1], NA_character_)
  if (is.na(ptu_col)) {
    return(
      tibble::tibble(
        Sample_ID = character(),
        plasmid_id = character(),
        PTU = character(),
        PTU_Score = numeric(),
        PTU_Host_Range = character(),
        PTU_Notes = character()
      )
    )
  }
  if (is.na(id_col)) {
    ptu_val <- as.character(d[[ptu_col]][1])
    if (!is.na(ptu_val) && ptu_val %in% c("-", "NA", "N/A", "")) {
      ptu_val <- NA_character_
    }
    return(
      tibble::tibble(
        Sample_ID = sid,
        plasmid_id = NA_character_,
        PTU = ptu_val,
        PTU_Score = if (!is.na(score_col)) suppressWarnings(as.numeric(d[[score_col]][1])) else NA_real_,
        PTU_Host_Range = if (!is.na(host_range_col)) as.character(d[[host_range_col]][1]) else NA_character_,
        PTU_Notes = if (!is.na(notes_col)) as.character(d[[notes_col]][1]) else NA_character_
      )
    )
  }
  tibble::tibble(
    Sample_ID = sid,
    plasmid_id = normalize_record_id(d[[id_col]]),
    PTU = as.character(d[[ptu_col]]),
    PTU_Score = if (!is.na(score_col)) suppressWarnings(as.numeric(d[[score_col]])) else NA_real_,
    PTU_Host_Range = if (!is.na(host_range_col)) as.character(d[[host_range_col]]) else NA_character_,
    PTU_Notes = if (!is.na(notes_col)) as.character(d[[notes_col]]) else NA_character_
  ) %>%
    dplyr::mutate(
      PTU = dplyr::if_else(PTU %in% c("-", "NA", "N/A", ""), NA_character_, PTU)
  ) %>%
    dplyr::filter(!(is.na(PTU) & (is.na(PTU_Notes) | PTU_Notes == "")))
}

run_copla_for_sample <- function(sid, sample_fasta, sample_dir, copla_runtime, force = FALSE) {
  if (!isTRUE(copla_runtime$ready)) {
    return(tibble::tibble(Sample_ID = character(), plasmid_id = character(), PTU = character(), PTU_Score = numeric(), PTU_Host_Range = character(), PTU_Notes = character()))
  }
  plasmid_fa <- extract_plasmid_fasta(sample_dir, sample_fasta)
  if (is.null(plasmid_fa)) {
    return(tibble::tibble(Sample_ID = character(), plasmid_id = character(), PTU = character(), PTU_Score = numeric(), PTU_Host_Range = character(), PTU_Notes = character()))
  }

  copla_dir <- file.path(sample_dir, "copla")
  fs::dir_create(copla_dir)
  if (isTRUE(force) && dir.exists(file.path(copla_dir, "output"))) {
    unlink(file.path(copla_dir, "output"), recursive = TRUE, force = TRUE)
  }
  out_dir <- file.path(copla_dir, "output")
  args <- c(
    "run", "-n", copla_runtime$conda_env,
    copla_runtime$python_bin,
    copla_runtime$script_path,
    plasmid_fa,
    copla_runtime$pickle_path,
    copla_runtime$fofn_path,
    out_dir,
    "-t", copla_runtime$topology
  )
  res <- tryCatch(
    run_external_tool(copla_runtime$conda, args = args, wd = copla_runtime$project_root),
    error = function(e) list(status = 1L, stdout = "", stderr = conditionMessage(e))
  )
  readr::write_lines(as.character(res$stdout), file.path(copla_dir, "copla.stdout.log"))
  readr::write_lines(as.character(res$stderr), file.path(copla_dir, "copla.stderr.log"))
  if (!isTRUE(res$status == 0)) {
    return(tibble::tibble(Sample_ID = character(), plasmid_id = character(), PTU = character(), PTU_Score = numeric(), PTU_Host_Range = character(), PTU_Notes = character()))
  }

  pred <- find_copla_prediction_file(out_dir)
  if (is.na(pred)) {
    return(tibble::tibble(Sample_ID = character(), plasmid_id = character(), PTU = character(), PTU_Score = numeric(), PTU_Host_Range = character(), PTU_Notes = character()))
  }
  out <- read_copla_prediction(pred, sid)
  out
}

merge_ptu_from_copla <- function(rows, ptu_rows) {
  if (nrow(rows) == 0) {
    return(rows)
  }
  if (nrow(ptu_rows) == 0) {
    return(rows)
  }
  rows2 <- rows %>%
    dplyr::mutate(plasmid_id_norm = normalize_record_id(plasmid_id))
  ptu2 <- ptu_rows %>%
    dplyr::mutate(plasmid_id_norm = normalize_record_id(plasmid_id))

  out <- rows2 %>%
    dplyr::left_join(
      ptu2 %>% dplyr::select(Sample_ID, plasmid_id_norm, PTU_copla = PTU, PTU_Score_copla = PTU_Score),
      by = c("Sample_ID", "plasmid_id_norm")
    ) %>%
    dplyr::left_join(
      ptu2 %>% dplyr::select(
        Sample_ID,
        plasmid_id_norm,
        PTU_Host_Range_copla = PTU_Host_Range,
        PTU_Notes_copla = PTU_Notes
      ),
      by = c("Sample_ID", "plasmid_id_norm")
    ) %>%
    dplyr::mutate(
      PTU = dplyr::coalesce(PTU_copla, PTU),
      PTU_Score = dplyr::coalesce(PTU_Score_copla, PTU_Score),
      PTU_Host_Range = dplyr::coalesce(PTU_Host_Range_copla, PTU_Host_Range),
      PTU_Notes = dplyr::coalesce(PTU_Notes_copla, PTU_Notes)
    ) %>%
    dplyr::select(
      -plasmid_id_norm,
      -PTU_copla,
      -PTU_Score_copla,
      -PTU_Host_Range_copla,
      -PTU_Notes_copla
    )

  fallback <- ptu_rows %>% dplyr::filter(is.na(plasmid_id) | plasmid_id == "")
  if (nrow(fallback) > 0) {
    first_ptu <- fallback$PTU[1]
    first_score <- fallback$PTU_Score[1]
    first_host <- fallback$PTU_Host_Range[1]
    first_notes <- fallback$PTU_Notes[1]
    out <- out %>%
      dplyr::mutate(
        PTU = dplyr::if_else(is.na(PTU) | PTU == "", first_ptu, PTU),
        PTU_Score = dplyr::if_else(is.na(PTU_Score), first_score, PTU_Score),
        PTU_Host_Range = dplyr::if_else(is.na(PTU_Host_Range) | PTU_Host_Range == "", first_host, PTU_Host_Range),
        PTU_Notes = dplyr::if_else(is.na(PTU_Notes) | PTU_Notes == "", first_notes, PTU_Notes)
      )
  }
  out
}

run_mobsuite <- function(manifest_df, cfg, out_dir, tools) {
  fs::dir_create(out_dir)
  emit_warn <- isTRUE(cfg$runtime$verbose)
  copla_runtime <- resolve_copla_runtime(cfg)
  if (emit_warn && !isTRUE(copla_runtime$ready)) {
    log_warn(glue::glue("[MOB-suite] COPLA environment not configured correctly ({copla_runtime$reason}); PTU will be NA."))
  }

  rows <- purrr::map_dfr(seq_len(nrow(manifest_df)), function(i) {
    sid <- manifest_df$Sample_ID[i]
    fasta <- manifest_df$Input_FASTA[i]
    sample_dir <- file.path(out_dir, sid)
    status <- "not_run"
    warn_reason <- NA_character_
    launch_error <- NA_character_
    log_info(glue::glue("[MOB-suite] ({i}/{nrow(manifest_df)}) {sid}"))
    if (!is.null(tools[["mob_recon"]]) && nzchar(tools[["mob_recon"]])) {
      mob_args <- c("-i", fasta, "-o", sample_dir, "-n", as.character(cfg$threads), "-s", sid)
      if (isTRUE(cfg$runtime$force)) {
        mob_args <- c(mob_args, "--force")
      }
      res <- tryCatch(
        run_external_tool(tools[["mob_recon"]], args = mob_args),
        error = function(e) {
          launch_error <<- conditionMessage(e)
          NULL
        }
      )
      status <- if (!is.null(res) && isTRUE(res$status == 0)) "ok" else "failed"
      if (!is.null(res)) {
        log_dir <- if (dir.exists(sample_dir)) sample_dir else out_dir
        readr::write_lines(res$stdout, file.path(log_dir, paste0(sid, ".mob_recon.stdout.log")))
        readr::write_lines(res$stderr, file.path(log_dir, paste0(sid, ".mob_recon.stderr.log")))
      }
      if (status == "failed") {
        warn_reason <- "mob_recon command failed; check *.mob_recon.stderr.log in 02_mobsuite."
      }
      if (!is.na(launch_error)) {
        log_dir <- if (dir.exists(sample_dir)) sample_dir else out_dir
        readr::write_lines(launch_error, file.path(log_dir, paste0(sid, ".mob_recon.launch_error.log")))
        warn_reason <- "mob_recon failed to launch; check *.mob_recon.launch_error.log in 02_mobsuite."
      }
    } else {
      warn_reason <- "mob_recon executable not available."
    }

    reports <- read_mob_reports(sample_dir)
    if (is.null(reports$mobtyper) && is.null(reports$contig) && is.na(warn_reason)) {
      warn_reason <- "No MOB-suite report found (expected mobtyper_results(.txt) or contig_report.txt)."
    }
    if (emit_warn && !is.na(warn_reason)) {
      log_warn(glue::glue("[MOB-suite] {sid}: {warn_reason}"))
    }
    per_sample <- build_mob_manifest_rows(sid, reports, status)
    ptu_rows <- run_copla_for_sample(
      sid = sid,
      sample_fasta = fasta,
      sample_dir = sample_dir,
      copla_runtime = copla_runtime,
      force = isTRUE(cfg$runtime$force)
    )
    merge_ptu_from_copla(per_sample, ptu_rows)
  })

  out_file_stage <- file.path(out_dir, "plasmid_manifest.tsv")
  out_file_root <- file.path(dirname(out_dir), "plasmid_manifest.tsv")
  readr::write_tsv(rows, out_file_stage)
  readr::write_tsv(rows, out_file_root)
  rows
}
