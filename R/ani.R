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

run_pairwise_ani <- function(sample_manifest, cfg, out_dir, tools) {
  related_cutoff <- dplyr::coalesce(cfg$ani$relatedness_threshold, cfg$ani$related_cutoff, 99)
  if (nrow(sample_manifest) < 2) return(empty_ani_results())
  pairs <- t(utils::combn(sample_manifest$Sample_ID, 2))
  res <- purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
    a <- pairs[i, 1]
    b <- pairs[i, 2]
    fa <- sample_manifest$Input_FASTA[match(a, sample_manifest$Sample_ID)]
    fb <- sample_manifest$Input_FASTA[match(b, sample_manifest$Sample_ID)]
    out_file <- file.path(out_dir, paste0(a, "__", b, ".tsv"))
    if (!is.null(tools[["fastani"]]) && nzchar(tools[["fastani"]])) {
      run_external_tool(
        tools[["fastani"]],
        c("--query", fa, "--ref", fb, "--minFraction", as.character(cfg$ani$min_fraction), "--fragLen", as.character(cfg$ani$frag_len), "--kmer", as.character(cfg$ani$kmer), "--threads", "1", "--output", out_file)
      )
      if (file.exists(out_file) && file.info(out_file)$size > 0) {
        d <- tryCatch(readr::read_tsv(out_file, col_names = FALSE, show_col_types = FALSE), error = function(e) NULL)
        if (!is.null(d) && nrow(d) > 0) {
          return(tibble::tibble(
            Sample_A = a, Sample_B = b,
            ANI = as.numeric(d[[3]][1]),
            FragmentsMapped = as.integer(d[[4]][1]),
            FragmentsTotal = as.integer(d[[5]][1]),
            Relatedness = ifelse(as.numeric(d[[3]][1]) > related_cutoff, "Related", "Unrelated")
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
