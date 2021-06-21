
# pipeline == begraph or gemBS
# for customs pipeline, pipeline needs to be a list
read_cpg_data <- function(fname, dataset_id = NULL, pipeline = "bedgraph",
  zero_based = TRUE, collapse_strands = TRUE,
  upper_cov_cutoff = list(method = "percentile", value = 0.999),
  reference = "hg19", align_to_reference = TRUE,
  use_cache = TRUE, update_cache = TRUE) {

  dt <- read_cache_file(cache_id = dataset_id, use_cache = use_cache)
  if (is.null(dt)) {
    pipeline <- process_pipeline(pipeline)
    dt <- read_data(fname, pipeline)
    dt <- process_data(dt, pipeline, collapse_strands, upper_cov_cutoff,
      reference, align_to_reference)
    write_cache_file(dt, dataset_id)
  }

  return(dt)
}


#### Internal Functions

# e.g.
# read_data(fname, "bedgraph")
# read_data(fname, "bedgraph", cmd = glue::glue("zcat {fname} | grep -v '^track'"))
read_data <- function(fname, pipeline, cmd = NULL, ...) {
  headers <- get_header_names(pipeline$cols)
  select_cols <- unlist(pipeline$cols, use.names = FALSE)

  logger::log_info("=> Reading from Source: {fname}")
  if (identical(cmd, NULL)) {
    dt <- data.table::fread(fname, col.names = headers, select = select_cols,
      ...)
  } else {
    dt <- data.table::fread(cmd = cmd, col.names = headers, 
      select = select_cols, ...)
  }
  return(dt)
}

process_data <- function(dt, pipeline, collapse_strands, upper_cov_cutoff,
  reference, align_to_reference) {
  remove_leading_chr(dt)
  convert_meth_to_frac(dt, pipeline)
  m_dt <- merge_strand_information(dt, pipeline, collapse_strands)

  if (! "cov" %in% names(m_dt)) {
    m_dt[, cov := meth_cov + unmeth_cov]
    m_dt[, .(meth_cov, unmeth_cov) := NULL]
  }

  mask_technical_noise(m_dt, upper_cov_cutoff)
  m_dt_all <- align_to_ref(m_dt, reference, align_to_reference)

  setcolorder(m_dt_all, c("chr", "start", "cov", "meth"))
  setkey(m_dt_all, chr, start)
  return(m_dt_all)
}







# custom format recognises the following keys:
# chr_col
# start_col
# cov_col
# beta_col
# strand_col
# n_meth
# n_unmeth
# fraction
process_pipeline <- function(format) {
  if (format == "bedgraph") {
    pipeline <- list(
      cols = list(
        chr = 1,
        start = 2,
        cov = 10,
        beta = 11,
        strand = 6
      ),
      fraction = TRUE
    )
  } else if (format == "gembs") {
    pipeline <- list(
      cols = list(
        chr = 1,
        start = 2,
        cov = 5,
        beta = 7
      ),
      fraction = TRUE
    )
  } else if (is.list(format)) {
    pipeline <- list(
      cols = list(
        chr = format$chr_col,
        start = format$start_col,
        cov = format$cov_col,
        beta = format$beta_col,
        strand = format$strand_col,
        meth_cov = format$meth_cov,
        unmeth_cov = format$unmeth_cov
      ),
      fraction = format$fraction
    )
  }
  # remove keys with a NULL value
  compact_pipeline <- pipeline[lengths(pipeline) != 0]
  return(compact_pipeline)
}

get_header_names <- function(pipeline) {
  standard_header_names <- list(
    chr = "chr",
    start = "start",
    cov = "cov",
    beta = "meth",
    strand = "strand",
    meth_cov = "meth_cov",
    unmeth_cov = "unmeth_cov"
  )
  headers <- sapply(names(pipeline), function (x) {
    standard_header_names[[x]]
  }, USE.NAMES = FALSE)
  # remove keys with a NULL value
  compact_headers <- headers[lengths(headers) != 0]
  return(compact_headers)
}

remove_leading_chr <- function(dt) {
  tmp_dt <- c(utils::head(dt)[["chr"]], utils::tail(dt)[["chr"]])
  presence_of_chr <- stringr::str_starts(tmp_dt, stringr::coll("chr"))
  if (any(presence_of_chr)) {
    logger::log_info("==> Removing `chr` from chromosome column...")
    dt[, chr := stringr::str_replace(chr, stringr::coll("chr"), "")]
  }
}

convert_meth_to_frac <- function(dt, pipeline) {
  if ((dt[, max(meth)] > 1) | identical(pipeline$fraction, FALSE)) {
    logger::log_info("==> Converting methylation to fraction...")
    dt[, meth := meth / 100]
  }
}

merge_strand_information <- function(dt, pipeline, collapse_strands) {
  if (!"strand" %in% names(dt) & !identical(collapse_strands, TRUE)) {
    return(dt)
  }
  logger::log_info("==> Merging Strands...")

  # move from single cov col to meth_cov + unmeth_cov
  if (! "meth_cov" %in% names(dt)) {
    dt[, meth_cov := round(cov * meth)]
    dt[, unmeth_cov := cov - meth_cov]
    set(dt, j = c("meth", "cov"), value = NULL)
  }

  dt[, start := fifelse(strand == "-", start, start - 1)]

  plus_dt <- dt[strand == "+"]
  neg_dt <- dt[strand == "-"]
  setkey(plus_dt, chr, start)
  setkey(neg_dt, chr, start)

  m_dt <- merge(plus_dt, neg_dt, all = TRUE)
  m_dt[is.na(meth_cov.x), meth_cov.x := 0]
  m_dt[is.na(meth_cov.y), meth_cov.y := 0]
  m_dt[is.na(unmeth_cov.x), unmeth_cov.x := 0]
  m_dt[is.na(unmeth_cov.y), unmeth_cov.y := 0]

  m_dt[, meth_cov := meth_cov.x + meth_cov.y]
  m_dt[, cov := meth_cov + unmeth_cov.x + unmeth_cov.y]
  m_dt[, meth := meth_cov / cov]

  set(m_dt, NULL, setdiff(names(m_dt), c("chr", "start", "cov", "meth")), NULL)

  return(m_dt[!is.na(cov) & !is.na(meth) & cov > 0, ])
}

# upper_cov_cutoff list(method = "percentile", value = 0.999)
#' @importFrom stats mad
mask_technical_noise <- function(dt, upper_cov_cutoff) {
  if (upper_cov_cutoff$method == "percentile") {
    max_cov <- dt[!is.na(cov), stats::quantile(cov, upper_cov_cutoff$value)]
  } else if (upper_cov_cutoff$method == "mad") {
    max_cov <- dt[!is.na(cov), mean(cov) + (3 * mad(cov))]
  }
  logger::log_info("==> Masking technical noise - Coverage >= {max_cov}")
  dt[cov >= max_cov, `:=`(cov = NA, meth = NA)]
}

align_to_ref <- function(dt, ref_name, align_to_reference) {
  if (identical(align_to_reference, TRUE)) {
    logger::log_info("==> Aligning data to reference dataset: {ref_name}")
    # check if reference is included in dataset 
    ref_cpg_data <- load_ref_data(ref_name)
    setkey(dt, chr, start)
    dt_all <- dt[ref_cpg_data]
    return(dt_all)
  } else {
    return(dt)
  }
}
