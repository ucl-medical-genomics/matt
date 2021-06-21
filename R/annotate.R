generate_annotations <- function(ref_cpgs, cache_key = "hg19_annot",
                                from_cache = TRUE, update_cache = TRUE) {
  fst_file <- get_fname(glue("ref_cpgs_annotations-{tolower(cache_key)}"),
                        "cache", ".fst")
  if (file.exists(fst_file) && identical(from_cache, TRUE)) {
    log_info("=> Reading from Cache: {fst_file}")
    dt <- read_fst(fst_file, as.data.table = TRUE)
    return(dt)
  }
  levels <- c("genes_intergenic", "genes_1to5kb", "genes_promoters",
    "genes_exons", "genes_firstexons", "genes_introns",
    "genes_exonintronboundaries", "genes_intronexonboundaries", "genes_cds",
    "genes_5UTRs", "genes_3UTRs", "lncrna_gencode", "enhancers_fantom",
    "cpg_inter", "cpg_shelves", "cpg_shores", "cpg_islands")
  levels <- paste0("hg19_", levels)

  annotations <- build_annotations(genome = "hg19", annotations = levels)
  ref_cpgs[, `:=`(end = start + 1, chr = paste0("chr", chr))]
  ref_cpgs_gr <- makeGRangesFromDataFrame(ref_cpgs)
  ref_cpgs[, `:=`(end = NULL, chr = str_sub(chr, 4L))]
  dt <- as.data.table(annotate_regions(ref_cpgs_gr, annotations,
          ignore.strand = TRUE, quiet = FALSE))

  rm_cols <- setdiff(names(dt), c("seqnames", "start", "annot.type"))
  set(dt, j = rm_cols, value = NULL)
  setnames(dt, c("seqnames", "annot.type"), c("chr", "type"))
  dt[, chr := str_sub(chr, 4L)]
  dt <- unique(dt)

  cpg_labels <- c("Inter CGI", "CpG Shelves", "CpG Shores", "CpG Islands")
  labels <- c("Intergenic", "1-5kb Upstream", "Promoter", "Exon", "First Exon",
    "Introns", "Exon-Intron\nBoundary", "Intron-Exon\nBoundary", "CDS", "5 UTR",
    "3 UTR", "LNC RNA", "Enhancer", cpg_labels)
  set(dt, j = "type", value = factor(dt[["type"]], levels = levels,
      labels = labels))
  dt[, annotation := fifelse(type %in% cpg_labels, "CpG", "Genomic")]
  setkey(dt, chr, start)
  if (identical(update_cache, TRUE)) {
    log_info("=> Writing to Cache: {fst_file}")
    write_fst(dt, path = fst_file)
  }
  return(dt)
}

annotate_dt <- function(dt, ref_cpgs_annot) {
  dt <- dt[!is.na(abs_diff)]
  merged <- merge(dt, ref_cpgs_annot, all.x = TRUE, by = c("chr", "start"))
  merged <- merged[!is.na(type)]
  merged_summary <- merged[, .(
    N = .N,
    y0 = min(abs_diff, na.rm = TRUE),
    y25 = quantile(abs_diff, 0.25, na.rm = TRUE),
    y50 = median(abs_diff, na.rm = TRUE),
    y75 = quantile(abs_diff, 0.75, na.rm = TRUE),
    y100 = max(abs_diff, na.rm = TRUE)
  ), by = .(type, annotation)]

  merged_summary[, `:=`(
    ymin = pmax(y0, (y25 - (y75 - y25) * 1.5)),
    ymax = pmin(y100, (y75 + (y75 - y25) * 1.5))
  ), by = .(type, annotation)]


  return(merged_summary)
}
