# Virulence BLAST module.

empty_virulence_hits <- function() {
  tibble::tibble(
    Sample_ID = character(), Gene = character(), Location = character(),
    Contig_ID = character(), Start = numeric(), End = numeric(), Strand = character(),
    Identity = numeric(), Coverage = numeric(), Evalue = numeric(), Bitscore = numeric(),
    Plasmid_ID = character(), Source_FASTA = character()
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

run_virulence_blast <- function(sample_manifest, cfg, out_dir, tools) {
  min_identity <- dplyr::coalesce(cfg$blast$min_identity, cfg$blast$identity, 95)
  min_coverage <- dplyr::coalesce(cfg$blast$min_coverage, cfg$blast$coverage, 95)
  evalue_cutoff <- dplyr::coalesce(cfg$blast$evalue, 1e-5)
  db_prefix <- resolve_virulence_db_prefix(cfg)
  if (is.null(tools[["blastn"]]) || !nzchar(tools[["blastn"]])) {
    stop("Missing required tool: blastn", call. = FALSE)
  }

  purrr::map_dfr(seq_len(nrow(sample_manifest)), function(i) {
    sid <- sample_manifest$Sample_ID[i]
    fasta <- sample_manifest$Input_FASTA[i]
    out_file <- file.path(out_dir, paste0(sid, "_blast.tsv"))
    res <- run_external_tool(
      tools[["blastn"]],
      c(
        "-task", cfg$blast$task,
        "-query", fasta,
        "-db", db_prefix,
        "-evalue", as.character(evalue_cutoff),
        "-outfmt", "6",
        "-out", out_file
      )
    )
    if (!isTRUE(res$status == 0)) {
      stop(
        paste0("blastn failed for sample ", sid, " with status ", res$status, ": ", res$stderr),
        call. = FALSE
      )
    }
    if (!file.exists(out_file) || file.info(out_file)$size == 0) return(empty_virulence_hits())
    raw <- tryCatch(readr::read_tsv(out_file, col_names = FALSE, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(raw) || nrow(raw) == 0) return(empty_virulence_hits())
    names(raw) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
    raw %>%
      dplyr::mutate(
        Coverage = (length / (abs(send - sstart) + 1)) * 100,
        Sample_ID = sid, Gene = sseqid, Location = "chromosome", Contig_ID = qseqid,
        Start = qstart, End = qend, Strand = dplyr::if_else(qstart <= qend, "+", "-"),
        Identity = pident, Evalue = evalue, Bitscore = bitscore, Plasmid_ID = NA_character_,
        Source_FASTA = fasta
      ) %>%
      dplyr::filter(
        Identity >= min_identity,
        Coverage >= min_coverage,
        Evalue <= evalue_cutoff
      ) %>%
      dplyr::select(Sample_ID, Gene, Location, Contig_ID, Start, End, Strand, Identity, Coverage, Evalue, Bitscore, Plasmid_ID, Source_FASTA)
  })
}
