# Synteny module: chromosomal flanking-CDS order around virulence genes
# (related strain pairs only).

read_prodigal_gff_cds <- function(gff_path) {
  empty <- tibble::tibble(
    Contig_ID = character(), Start = integer(), End = integer(), Strand = character(),
    CDS_ID = character(), FAA_ID = character()
  )
  if (is.null(gff_path) || !nzchar(gff_path) || !file.exists(gff_path) || file.info(gff_path)$size == 0) {
    return(empty)
  }
  lines <- readr::read_lines(gff_path)
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  if (length(lines) == 0) {
    return(empty)
  }
  parts <- stringr::str_split_fixed(lines, "\t", 9)
  keep <- parts[, 3] == "CDS"
  parts <- parts[keep, , drop = FALSE]
  if (nrow(parts) == 0) {
    return(empty)
  }
  attrs <- parts[, 9]
  cds_id <- stringr::str_match(attrs, "(?:^|;)ID=([^;]+)")[, 2]
  contig <- normalize_contig_name(parts[, 1])
  faa_id <- paste0(contig, "_", stringr::str_match(cds_id, "_([0-9]+)$")[, 2])
  tibble::tibble(
    Contig_ID = contig,
    Start = as.integer(parts[, 4]),
    End = as.integer(parts[, 5]),
    Strand = parts[, 7],
    CDS_ID = cds_id,
    FAA_ID = faa_id
  ) %>%
    dplyr::arrange(Contig_ID, Start, End)
}

read_prodigal_faa <- function(faa_path) {
  if (is.null(faa_path) || !nzchar(faa_path) || !file.exists(faa_path) || file.info(faa_path)$size == 0) {
    return(list())
  }
  lines <- readr::read_lines(faa_path)
  seqs <- list()
  cur_id <- NULL
  cur_seq <- character()
  flush <- function() {
    if (!is.null(cur_id)) {
      seqs[[cur_id]] <<- paste(cur_seq, collapse = "")
    }
  }
  for (ln in lines) {
    if (startsWith(ln, ">")) {
      flush()
      hdr <- sub("^>", "", ln)
      cur_id <- strsplit(hdr, "\\s+", perl = TRUE)[[1]][1]
      cur_seq <- character()
    } else if (!is.null(cur_id)) {
      cur_seq <- c(cur_seq, ln)
    }
  }
  flush()
  seqs
}

resolve_annotation_paths <- function(sample_id, annotation_dir) {
  sample_dir <- file.path(annotation_dir, sample_id)
  list(
    gff = file.path(sample_dir, paste0(sample_id, ".prodigal.gff")),
    faa = file.path(sample_dir, paste0(sample_id, ".prodigal.faa"))
  )
}

pick_chromosomal_virulence_hit <- function(virulence_hits, sample_id, gene) {
  hits <- virulence_hits %>%
    dplyr::filter(Sample_ID == sample_id, Gene == gene, Location == "chromosome")
  if (nrow(hits) == 0) {
    return(NULL)
  }
  hits %>%
    dplyr::arrange(dplyr::desc(Bitscore), dplyr::desc(Identity), dplyr::desc(Coverage)) %>%
    dplyr::slice(1)
}

interval_overlap <- function(a_start, a_end, b_start, b_end) {
  left <- pmax(pmin(a_start, a_end), pmin(b_start, b_end))
  right <- pmin(pmax(a_start, a_end), pmax(b_start, b_end))
  pmax(0L, right - left + 1L)
}

find_overlapping_cds <- function(cds_tbl, contig_id, hit_start, hit_end) {
  on_ctg <- cds_tbl %>% dplyr::filter(Contig_ID == contig_id)
  if (nrow(on_ctg) == 0) {
    return(NULL)
  }
  ov <- interval_overlap(on_ctg$Start, on_ctg$End, hit_start, hit_end)
  if (all(ov <= 0)) {
    mid_hit <- (hit_start + hit_end) / 2
    mid_cds <- (on_ctg$Start + on_ctg$End) / 2
    idx <- which.min(abs(mid_cds - mid_hit))
    return(list(row = on_ctg[idx, , drop = FALSE], index = idx))
  }
  idx <- which.max(ov)
  list(row = on_ctg[idx, , drop = FALSE], index = idx)
}

#' Ordered flanking CDS around a chromosomal virulence CDS.
#' Order is genomic (by Start). Upstream = lower coordinates; Downstream = higher.
#' Returned table keeps Position: -n_up .. -1, +1 .. +n_down.
extract_ordered_flanking_cds <- function(cds_tbl, contig_id, center_index, n_up, n_down) {
  on_ctg <- cds_tbl %>% dplyr::filter(Contig_ID == contig_id)
  if (nrow(on_ctg) == 0 || is.null(center_index) || center_index < 1 || center_index > nrow(on_ctg)) {
    return(tibble::tibble())
  }
  up_idx <- rev(seq(max(1L, center_index - n_up), center_index - 1L))
  # reverse so Position goes -1 (nearest), -2, ... then we re-order by Position ascending
  down_idx <- seq(center_index + 1L, min(nrow(on_ctg), center_index + n_down))

  up_rows <- if (length(up_idx) == 0) {
    tibble::tibble()
  } else {
    on_ctg[up_idx, , drop = FALSE] %>%
      dplyr::mutate(
        Flank_Role = "upstream",
        Position = -seq_len(dplyr::n())
      )
  }
  down_rows <- if (length(down_idx) == 0) {
    tibble::tibble()
  } else {
    on_ctg[down_idx, , drop = FALSE] %>%
      dplyr::mutate(
        Flank_Role = "downstream",
        Position = seq_len(dplyr::n())
      )
  }
  dplyr::bind_rows(up_rows, down_rows) %>%
    dplyr::arrange(Position)
}

write_flank_faa <- function(flank_cds, faa_seqs, out_path) {
  if (nrow(flank_cds) == 0) {
    return(FALSE)
  }
  chunks <- character()
  for (i in seq_len(nrow(flank_cds))) {
    fid <- flank_cds$FAA_ID[i]
    seq <- faa_seqs[[fid]]
    if (is.null(seq) || !nzchar(seq)) {
      next
    }
    # encode genomic order position into the FASTA id for later mapping
    pos <- flank_cds$Position[i]
    chunks <- c(chunks, paste0(">", fid, "|pos=", pos), gsub("\\*$", "", seq))
  }
  if (length(chunks) == 0) {
    return(FALSE)
  }
  readr::write_lines(chunks, out_path)
  TRUE
}

blastp_best_hits <- function(query_faa, subject_faa, identity_cutoff, coverage_cutoff, blastp_bin) {
  empty <- tibble::tibble(qseqid = character(), sseqid = character(), pident = numeric(), bitscore = numeric())
  if (!file.exists(query_faa) || !file.exists(subject_faa) || file.info(query_faa)$size == 0 || file.info(subject_faa)$size == 0) {
    return(empty)
  }
  if (is.null(blastp_bin) || !nzchar(blastp_bin)) {
    blastp_bin <- Sys.which("blastp")
  }
  if (!nzchar(blastp_bin)) {
    return(NULL)
  }
  out_tsv <- tempfile(fileext = ".blastp.tsv")
  on.exit(unlink(out_tsv), add = TRUE)
  res <- tryCatch(
    run_external_tool(
      blastp_bin,
      args = c(
        "-query", query_faa,
        "-subject", subject_faa,
        "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
        "-evalue", "1e-5",
        "-max_target_seqs", "5",
        "-out", out_tsv
      )
    ),
    error = function(e) NULL
  )
  if (is.null(res) || !isTRUE(res$status == 0) || !file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
    return(empty)
  }
  raw <- tryCatch(
    readr::read_tsv(
      out_tsv,
      col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(raw) || nrow(raw) == 0) {
    return(empty)
  }
  raw %>%
    dplyr::mutate(
      qcov = (length / qlen) * 100,
      scov = (length / slen) * 100
    ) %>%
    dplyr::filter(
      pident >= identity_cutoff,
      qcov >= coverage_cutoff,
      scov >= coverage_cutoff
    ) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(qseqid, sseqid, pident, bitscore)
}

parse_pos_from_faa_id <- function(x) {
  as.integer(stringr::str_match(x, "\\|pos=(-?[0-9]+)$")[, 2])
}

strip_pos_from_faa_id <- function(x) {
  sub("\\|pos=-?[0-9]+$", "", x)
}

#' Count order-preserving homologous flanking CDS between A and B.
#'
#' Uses blastp to map A flank CDS -> B flank CDS, then scores how many
#' homologs preserve genomic order (via LCS of mapped B positions along A order).
#' Also accepts the reverse-orientation LCS (B positions reversed).
score_flank_order_synteny <- function(flank_a, flank_b, faa_a, faa_b, identity_cutoff, coverage_cutoff, blastp_bin) {
  tmp_dir <- tempfile("kpvir_synteny_ord_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  path_a <- file.path(tmp_dir, "A_flank.faa")
  path_b <- file.path(tmp_dir, "B_flank.faa")
  if (!write_flank_faa(flank_a, faa_a, path_a) || !write_flank_faa(flank_b, faa_b, path_b)) {
    return(NA_integer_)
  }

  hits <- blastp_best_hits(path_a, path_b, identity_cutoff, coverage_cutoff, blastp_bin)
  if (is.null(hits)) {
    return(NA_integer_)
  }
  if (nrow(hits) == 0) {
    return(0L)
  }

  mapped <- hits %>%
    dplyr::mutate(
      q_pos = parse_pos_from_faa_id(qseqid),
      s_pos = parse_pos_from_faa_id(sseqid),
      q_id = strip_pos_from_faa_id(qseqid),
      s_id = strip_pos_from_faa_id(sseqid)
    ) %>%
    dplyr::filter(!is.na(q_pos), !is.na(s_pos)) %>%
    dplyr::arrange(q_pos)

  if (nrow(mapped) == 0) {
    return(0L)
  }

  # Order along A genomic flank; B positions of matched homologs.
  # Longest increasing/decreasing subsequence = order-preserving count
  # (decreasing allows local neighborhood inversion).
  b_along_a <- mapped$s_pos
  lis <- longest_monotonic_subseq_length(b_along_a, increasing = TRUE)
  lds <- longest_monotonic_subseq_length(b_along_a, increasing = FALSE)
  as.integer(max(lis, lds))
}

longest_monotonic_subseq_length <- function(x, increasing = TRUE) {
  if (length(x) == 0) {
    return(0L)
  }
  n <- length(x)
  dp <- rep(1L, n)
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {
      ok <- if (increasing) x[j] < x[i] else x[j] > x[i]
      # allow equal? no — positions should be unique
      if (ok) {
        dp[i] <- max(dp[i], dp[j] + 1L)
      }
    }
  }
  as.integer(max(dp))
}

build_chromosomal_flank_context <- function(sample_id, gene, virulence_hits, annotation_dir, n_up, n_down) {
  hit <- pick_chromosomal_virulence_hit(virulence_hits, sample_id, gene)
  if (is.null(hit) || nrow(hit) == 0) {
    return(NULL)
  }
  paths <- resolve_annotation_paths(sample_id, annotation_dir)
  cds <- read_prodigal_gff_cds(paths$gff)
  faa <- read_prodigal_faa(paths$faa)
  if (nrow(cds) == 0 || length(faa) == 0) {
    return(NULL)
  }
  ov <- find_overlapping_cds(cds, hit$Contig_ID[1], as.integer(hit$Start[1]), as.integer(hit$End[1]))
  if (is.null(ov)) {
    return(NULL)
  }
  flank <- extract_ordered_flanking_cds(cds, hit$Contig_ID[1], ov$index, n_up, n_down)
  list(
    hit = hit,
    center = ov$row,
    flank = flank,
    faa = faa,
    paths = paths
  )
}

#' Compute chromosomal flanking-CDS order synteny for related strains.
#'
#' Only runs when:
#' - pair_record$Pair_Related is TRUE
#' - both samples have a **chromosome** hit for the gene
#'
#' Synteny_Index = number of order-preserving homologous flanking CDS
#' (upstream + downstream Prodigal neighborhoods), allowing local inversion.
#'
#' @param pair_record Named list with Sample_A, Sample_B, Gene, Location_A,
#'   Location_B, and Pair_Related.
#' @param annotation_ctx List with annotation_dir, virulence_hits, blastp.
#' @param cfg Pipeline config.
#' @return Integer score, or NA_integer_ if not applicable / not computable.
compute_synteny_index <- function(pair_record, annotation_ctx, cfg) {
  if (is.null(annotation_ctx) || is.null(annotation_ctx$annotation_dir) || is.null(annotation_ctx$virulence_hits)) {
    return(NA_integer_)
  }
  # Related strains only.
  if (!isTRUE(pair_record[["Pair_Related"]])) {
    return(NA_integer_)
  }

  loc_a <- as.character(pair_record[["Location_A"]])
  loc_b <- as.character(pair_record[["Location_B"]])
  # Chromosomal virulence gene on both sides (chromosome or both).
  chrom_a <- loc_a %in% c("chromosome", "both")
  chrom_b <- loc_b %in% c("chromosome", "both")
  if (!chrom_a || !chrom_b) {
    return(NA_integer_)
  }

  annotation_dir <- annotation_ctx$annotation_dir
  virulence_hits <- annotation_ctx$virulence_hits
  blastp_bin <- annotation_ctx$blastp
  if (is.null(blastp_bin) || !nzchar(blastp_bin)) {
    blastp_bin <- Sys.which("blastp")
  }

  syn <- cfg$synteny
  n_up <- as.integer(dplyr::coalesce(syn$upstream_cds, syn$upstream_window, 5L))
  n_down <- as.integer(dplyr::coalesce(syn$downstream_cds, syn$downstream_window, 5L))
  identity_cutoff <- as.numeric(dplyr::coalesce(syn$protein_identity, 70))
  coverage_cutoff <- as.numeric(dplyr::coalesce(syn$protein_coverage, 70))

  sample_a <- as.character(pair_record[["Sample_A"]])
  sample_b <- as.character(pair_record[["Sample_B"]])
  gene <- as.character(pair_record[["Gene"]])

  ctx_a <- build_chromosomal_flank_context(sample_a, gene, virulence_hits, annotation_dir, n_up, n_down)
  ctx_b <- build_chromosomal_flank_context(sample_b, gene, virulence_hits, annotation_dir, n_up, n_down)
  if (is.null(ctx_a) || is.null(ctx_b) || nrow(ctx_a$flank) == 0 || nrow(ctx_b$flank) == 0) {
    return(NA_integer_)
  }

  score_flank_order_synteny(
    ctx_a$flank, ctx_b$flank, ctx_a$faa, ctx_b$faa,
    identity_cutoff, coverage_cutoff, blastp_bin
  )
}
