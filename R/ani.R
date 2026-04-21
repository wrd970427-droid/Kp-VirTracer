# ANI module.

empty_ani_results <- function() {
  tibble::tibble(
    Sample_A = character(),
    Sample_B = character(),
    ANI = numeric(),
    FragmentsMapped = integer(),
    FragmentsTotal = integer(),
    Relatedness = character()
  )
}

normalize_fasta_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  stringr::str_split_fixed(x, "\\s+", 2)[, 1]
}

read_fasta_records <- function(path) {
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

write_fasta_records <- function(headers, seqs, path) {
  if (length(headers) == 0 || length(seqs) == 0) {
    return(FALSE)
  }
  lines <- unlist(purrr::map2(headers, seqs, ~ c(paste0(">", .x), .y)), use.names = FALSE)
  readr::write_lines(lines, path)
  TRUE
}

read_sample_chromosome_ids <- function(sample_id, mobsuite_dir) {
  if (is.null(mobsuite_dir) || !nzchar(mobsuite_dir)) {
    return(character())
  }
  f <- file.path(mobsuite_dir, sample_id, "contig_report.txt")
  if (!file.exists(f) || file.info(f)$size == 0) {
    return(character())
  }
  d <- tryCatch(readr::read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) {
    return(character())
  }
  nm <- names(d)
  if (!("molecule_type" %in% nm) || !("contig_id" %in% nm)) {
    return(character())
  }
  d %>%
    dplyr::mutate(molecule_type = tolower(as.character(molecule_type))) %>%
    dplyr::filter(molecule_type == "chromosome") %>%
    dplyr::pull(contig_id) %>%
    normalize_fasta_id() %>%
    unique()
}

build_chromosome_manifest <- function(sample_manifest, mobsuite_dir, out_dir) {
  chr_dir <- file.path(out_dir, "chromosome_fasta")
  fs::dir_create(chr_dir)
  out <- purrr::map_dfr(seq_len(nrow(sample_manifest)), function(i) {
    sid <- sample_manifest$Sample_ID[i]
    in_fa <- sample_manifest$Input_FASTA[i]
    chr_ids <- read_sample_chromosome_ids(sid, mobsuite_dir)
    if (length(chr_ids) == 0) {
      log_warn(glue::glue("[ANI] {sid}: no chromosome contigs detected from MOB-suite; skipped."))
      return(tibble::tibble())
    }
    rec <- read_fasta_records(in_fa)
    if (length(rec$headers) == 0) {
      log_warn(glue::glue("[ANI] {sid}: empty FASTA; skipped."))
      return(tibble::tibble())
    }
    hdr_norm <- normalize_fasta_id(rec$headers)
    keep <- hdr_norm %in% chr_ids
    if (!any(keep)) {
      log_warn(glue::glue("[ANI] {sid}: chromosome IDs not found in FASTA headers; skipped."))
      return(tibble::tibble())
    }
    out_fa <- file.path(chr_dir, paste0(sid, ".chromosome.fna"))
    ok <- write_fasta_records(rec$headers[keep], rec$seqs[keep], out_fa)
    if (!ok || !file.exists(out_fa) || file.info(out_fa)$size == 0) {
      log_warn(glue::glue("[ANI] {sid}: failed to write chromosome FASTA; skipped."))
      return(tibble::tibble())
    }
    tibble::tibble(Sample_ID = sid, Input_FASTA = out_fa)
  })
  out
}

run_pairwise_ani <- function(sample_manifest, cfg, out_dir, tools) {
  related_cutoff <- dplyr::coalesce(cfg$ani$relatedness_threshold, cfg$ani$related_cutoff, 99)
  min_fraction <- dplyr::coalesce(cfg$ani$min_fraction, 0.2)
  frag_len <- dplyr::coalesce(cfg$ani$frag_len, 3000)
  kmer <- dplyr::coalesce(cfg$ani$kmer, 16)
  if (nrow(sample_manifest) < 2) return(empty_ani_results())
  pairs <- t(utils::combn(sample_manifest$Sample_ID, 2))
  res <- purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
    a <- pairs[i, 1]
    b <- pairs[i, 2]
    log_info(glue::glue("[ANI] ({i}/{nrow(pairs)}) {a} vs {b}"))
    fa <- sample_manifest$Input_FASTA[match(a, sample_manifest$Sample_ID)]
    fb <- sample_manifest$Input_FASTA[match(b, sample_manifest$Sample_ID)]
    out_file <- file.path(out_dir, paste0(a, "__", b, ".tsv"))
    if (!is.null(tools[["fastani"]]) && nzchar(tools[["fastani"]])) {
      run_external_tool(
        tools[["fastani"]],
        c(
          "--query", fa,
          "--ref", fb,
          "--minFraction", as.character(min_fraction),
          "--fragLen", as.character(frag_len),
          "--kmer", as.character(kmer),
          "--threads", "1",
          "--output", out_file
        )
      )
      if (file.exists(out_file) && file.info(out_file)$size > 0) {
        d <- tryCatch(readr::read_tsv(out_file, col_names = FALSE, show_col_types = FALSE), error = function(e) NULL)
        if (!is.null(d) && nrow(d) > 0) {
          return(tibble::tibble(
            Sample_A = a, Sample_B = b,
            ANI = as.numeric(d[[3]][1]),
            FragmentsMapped = as.integer(d[[4]][1]),
            FragmentsTotal = as.integer(d[[5]][1]),
            Relatedness = ifelse(as.numeric(d[[3]][1]) >= related_cutoff, "Related", "Unrelated")
          ))
        }
      }
    }
    tibble::tibble(
      Sample_A = a, Sample_B = b, ANI = NA_real_, FragmentsMapped = NA_integer_, FragmentsTotal = NA_integer_, Relatedness = "Unknown"
    )
  })
  readr::write_tsv(res, file.path(out_dir, "ani_results.tsv"))
  res
}
