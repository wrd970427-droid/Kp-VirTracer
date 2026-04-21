# Virulence BLAST module.

empty_virulence_hits <- function() {
  tibble::tibble(
    Sample_ID = character(), Gene = character(), Location = character(),
    Contig_ID = character(), Start = numeric(), End = numeric(), Strand = character(),
    Identity = numeric(), Coverage = numeric(), Evalue = numeric(), Bitscore = numeric(),
    Plasmid_ID = character(), Source_FASTA = character()
  )
}

empty_raw_blast_hits <- function() {
  tibble::tibble(
    qseqid = character(), sseqid = character(), pident = numeric(), length = numeric(),
    mismatch = numeric(), gapopen = numeric(), qstart = numeric(), qend = numeric(),
    sstart = numeric(), send = numeric(), evalue = numeric(), bitscore = numeric()
  )
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

has_blast_db_component <- function(db_prefix, ext) {
  exact <- paste0(db_prefix, ext)
  if (file.exists(exact)) {
    return(TRUE)
  }
  db_dir <- dirname(db_prefix)
  db_base <- basename(db_prefix)
  if (!dir.exists(db_dir)) {
    return(FALSE)
  }
  ext_no_dot <- substring(ext, 2)
  pattern <- paste0("^", escape_regex(db_base), "\\.[0-9]+\\.", escape_regex(ext_no_dot), "$")
  any(grepl(pattern, list.files(db_dir)))
}

normalize_contig_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  stringr::str_split_fixed(x, "\\s+", 2)[, 1]
}

read_mob_contig_locations <- function(sample_id, mobsuite_dir) {
  if (is.null(mobsuite_dir) || !nzchar(mobsuite_dir)) {
    return(NULL)
  }
  contig_report <- file.path(mobsuite_dir, sample_id, "contig_report.txt")
  if (!file.exists(contig_report) || file.info(contig_report)$size == 0) {
    return(NULL)
  }
  df <- tryCatch(readr::read_tsv(contig_report, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  nm <- names(df)
  if (!("molecule_type" %in% nm)) {
    return(NULL)
  }
  id_col <- dplyr::coalesce(
    dplyr::if_else("contig_id" %in% nm, "contig_id", NA_character_),
    dplyr::if_else("id" %in% nm, "id", NA_character_),
    dplyr::if_else("sample_id" %in% nm, "sample_id", NA_character_)
  )
  if (is.na(id_col)) {
    return(NULL)
  }
  df %>%
    dplyr::transmute(
      Contig_ID = normalize_contig_name(.data[[id_col]]),
      molecule_type = tolower(as.character(.data[["molecule_type"]]))
    ) %>%
    dplyr::mutate(
      Location = dplyr::case_when(
        molecule_type == "plasmid" ~ "plasmid",
        molecule_type == "chromosome" ~ "chromosome",
        TRUE ~ "chromosome"
      )
    ) %>%
    dplyr::distinct(Contig_ID, .keep_all = TRUE) %>%
    dplyr::select(Contig_ID, Location)
}

resolve_virulence_db_prefix <- function(cfg) {
  db_prefix <- cfg$blast$virulence_db
  if (is.null(db_prefix) || !nzchar(db_prefix)) {
    stop(
      paste0(
        "Missing virulence BLAST database prefix. ",
        "Provide --virulence-db with the prefix from `makeblastdb -out`.\n",
        "Example:\n",
        "  makeblastdb -in /path/to/virulence_genes.fasta -dbtype nucl -out /path/to/db/virulence_db\n",
        "  Rscript run_kpvirtracer.R --input <dir> --output <dir> --virulence-db /path/to/db/virulence_db"
      ),
      call. = FALSE
    )
  }

  db_prefix <- normalizePath(db_prefix, winslash = "/", mustWork = FALSE)
  if (file.exists(db_prefix) && grepl("\\.(fa|fasta|fna)$", db_prefix, ignore.case = TRUE)) {
    stop(
      paste0(
        "--virulence-db must be a BLAST database prefix, not a FASTA file: ", db_prefix, "\n",
        "Build DB first, for example:\n",
        "  makeblastdb -in ", db_prefix, " -dbtype nucl -out /path/to/db/virulence_db"
      ),
      call. = FALSE
    )
  }

  required_ext <- c(".nin", ".nhr", ".nsq")
  missing_ext <- required_ext[!vapply(required_ext, function(ext) has_blast_db_component(db_prefix, ext), logical(1))]
  if (length(missing_ext) > 0) {
    stop(
      paste0(
        "Invalid BLAST DB prefix: ", db_prefix, ". Missing index component(s): ",
        paste(missing_ext, collapse = ", "), ".\n",
        "Expected files like `", db_prefix, ".nin/.nhr/.nsq` ",
        "or split indexes like `", basename(db_prefix), ".00.nin`."
      ),
      call. = FALSE
    )
  }
  db_prefix
}

read_blast_table <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(empty_raw_blast_hits())
  }
  d <- tryCatch(readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) {
    return(empty_raw_blast_hits())
  }
  names(d) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  d
}

run_virulence_blast <- function(sample_manifest, cfg, out_dir, tools, mobsuite_dir = NULL) {
  min_identity <- dplyr::coalesce(cfg$blast$min_identity, cfg$blast$identity, 95)
  min_coverage <- dplyr::coalesce(cfg$blast$min_coverage, cfg$blast$coverage, 95)
  evalue_cutoff <- dplyr::coalesce(cfg$blast$evalue, 1e-5)
  blast_task <- dplyr::coalesce(cfg$blast$task, "blastn")
  db_prefix <- resolve_virulence_db_prefix(cfg)
  if (is.null(tools[["blastn"]]) || !nzchar(tools[["blastn"]])) {
    stop("Missing required tool: blastn", call. = FALSE)
  }

  logs_dir <- file.path(out_dir, "logs")
  raw_dir <- file.path(out_dir, "raw")
  filtered_dir <- file.path(out_dir, "filtered")
  by_loc_dir <- file.path(out_dir, "by_location")
  summary_dir <- file.path(out_dir, "summary")
  purrr::walk(c(logs_dir, raw_dir, filtered_dir, by_loc_dir, summary_dir), fs::dir_create)

  filtered_hits <- purrr::map_dfr(seq_len(nrow(sample_manifest)), function(i) {
    sid <- sample_manifest$Sample_ID[i]
    fasta <- sample_manifest$Input_FASTA[i]
    log_info(glue::glue("[BLAST] ({i}/{nrow(sample_manifest)}) {sid}"))

    raw_out <- file.path(raw_dir, paste0(sid, ".blast.raw.tsv"))
    filtered_out <- file.path(filtered_dir, paste0(sid, ".blast.filtered.tsv"))
    chrom_out <- file.path(by_loc_dir, paste0(sid, ".chromosome_hits.tsv"))
    plasmid_out <- file.path(by_loc_dir, paste0(sid, ".plasmid_hits.tsv"))

    run_blast_task <- function(task_name) {
      run_external_tool(
        tools[["blastn"]],
        c(
          "-task", task_name,
          "-query", fasta,
          "-db", db_prefix,
          "-evalue", as.character(evalue_cutoff),
          "-outfmt", "6",
          "-out", raw_out
        )
      )
    }

    res <- run_blast_task(blast_task)
    readr::write_lines(res$stdout, file.path(logs_dir, paste0(sid, ".blast.stdout.log")))
    readr::write_lines(res$stderr, file.path(logs_dir, paste0(sid, ".blast.stderr.log")))
    if (!isTRUE(res$status == 0)) {
      stop(
        paste0("blastn failed for sample ", sid, " with status ", res$status, ": ", res$stderr),
        call. = FALSE
      )
    }

    if ((!file.exists(raw_out) || file.info(raw_out)$size == 0) && identical(blast_task, "blastn-short")) {
      log_warn(glue::glue("[BLAST] {sid}: no hits with blastn-short, retrying with blastn."))
      res_retry <- run_blast_task("blastn")
      readr::write_lines(res_retry$stdout, file.path(logs_dir, paste0(sid, ".blast.retry.stdout.log")))
      readr::write_lines(res_retry$stderr, file.path(logs_dir, paste0(sid, ".blast.retry.stderr.log")))
      if (!isTRUE(res_retry$status == 0)) {
        stop(
          paste0("blastn retry failed for sample ", sid, " with status ", res_retry$status, ": ", res_retry$stderr),
          call. = FALSE
        )
      }
    }

    raw <- read_blast_table(raw_out)
    if (nrow(raw) == 0) {
      readr::write_tsv(empty_virulence_hits(), filtered_out)
      readr::write_tsv(empty_virulence_hits(), chrom_out)
      readr::write_tsv(empty_virulence_hits(), plasmid_out)
      return(empty_virulence_hits())
    }

    loc_map <- read_mob_contig_locations(sid, mobsuite_dir)
    raw2 <- raw %>%
      dplyr::mutate(qseqid_norm = normalize_contig_name(qseqid))
    if (!is.null(loc_map)) {
      raw2 <- raw2 %>%
        dplyr::left_join(loc_map, by = c("qseqid_norm" = "Contig_ID")) %>%
        dplyr::mutate(Location = dplyr::coalesce(Location, "chromosome"))
    } else {
      raw2 <- raw2 %>%
        dplyr::mutate(Location = "chromosome")
    }

    filtered <- raw2 %>%
      dplyr::mutate(
        Coverage = (length / (abs(send - sstart) + 1)) * 100,
        Sample_ID = sid,
        Gene = sseqid,
        Contig_ID = qseqid_norm,
        Start = qstart,
        End = qend,
        Strand = dplyr::if_else(qstart <= qend, "+", "-"),
        Identity = pident,
        Evalue = evalue,
        Bitscore = bitscore,
        Plasmid_ID = NA_character_,
        Source_FASTA = fasta
      ) %>%
      dplyr::filter(
        Identity >= min_identity,
        Coverage >= min_coverage,
        Evalue <= evalue_cutoff
      ) %>%
      dplyr::select(Sample_ID, Gene, Location, Contig_ID, Start, End, Strand, Identity, Coverage, Evalue, Bitscore, Plasmid_ID, Source_FASTA)

    readr::write_tsv(filtered, filtered_out)
    readr::write_tsv(filtered %>% dplyr::filter(Location == "chromosome"), chrom_out)
    readr::write_tsv(filtered %>% dplyr::filter(Location == "plasmid"), plasmid_out)
    filtered
  })

  sample_profile <- sample_manifest %>%
    dplyr::select(Sample_ID) %>%
    dplyr::left_join(
      filtered_hits %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::summarise(
          Virulence_Genes_All = paste(sort(unique(Gene)), collapse = ","),
          Chromosomal_Genes = paste(sort(unique(Gene[Location == "chromosome"])), collapse = ","),
          Plasmid_Genes = paste(sort(unique(Gene[Location == "plasmid"])), collapse = ","),
          n_chrom_hits = sum(Location == "chromosome", na.rm = TRUE),
          n_plasmid_hits = sum(Location == "plasmid", na.rm = TRUE),
          .groups = "drop"
        ),
      by = "Sample_ID"
    ) %>%
    dplyr::mutate(
      Virulence_Genes_All = dplyr::coalesce(Virulence_Genes_All, ""),
      Chromosomal_Genes = dplyr::coalesce(Chromosomal_Genes, ""),
      Plasmid_Genes = dplyr::coalesce(Plasmid_Genes, ""),
      n_chrom_hits = dplyr::coalesce(n_chrom_hits, 0L),
      n_plasmid_hits = dplyr::coalesce(n_plasmid_hits, 0L),
      Virulence_Location = dplyr::case_when(
        n_chrom_hits > 0 & n_plasmid_hits > 0 ~ "both",
        n_chrom_hits > 0 ~ "chromosome",
        n_plasmid_hits > 0 ~ "plasmid",
        TRUE ~ "absent"
      ),
      Virulence_Type = dplyr::case_when(
        Virulence_Location == "both" ~ "pc-hvKp",
        Virulence_Location == "chromosome" ~ "c-hvKp",
        Virulence_Location == "plasmid" ~ "p-hvKp",
        TRUE ~ "nKp"
      )
    )

  readr::write_tsv(filtered_hits, file.path(summary_dir, "virulence_hits.tsv"))
  readr::write_tsv(sample_profile, file.path(summary_dir, "sample_virulence_profile.tsv"))
  readr::write_tsv(
    sample_profile %>% dplyr::select(Sample_ID, Virulence_Location, Virulence_Type, n_chrom_hits, n_plasmid_hits),
    file.path(summary_dir, "virulence_type_calls.tsv")
  )

  filtered_hits
}
