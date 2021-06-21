
load_ref_data <- function(key, ...) {
  ref_file <- system.file("extdata", glue::glue("reference-{key}.fst"),
    package = "matt")

  # check local cache
  if (!fs::file_exists(ref_file)) {
    ref_file <- fs::path(.matt_env[["cache_dir"]], "reference-{key}.fst")
  }

  # check for file again
  if (fs::file_exists(ref_file)) {
    dt <- fst::read_fst(ref_file, as.data.table = TRUE)
  } else {
    dt <- get_ref_cpgs(...)
  }
  return(dt)
}

get_ref_cpgs <- function(key, ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                         all_chromosomes = c(seq(1, 22), "X", "Y"),
                         dir = .matt_env[["cache_dir"]],
                         update_cache = TRUE) {
  all_cpgs <- extract_cpg_sites(ref_genome = ref_genome)
  ref_cpgs <- all_cpgs$cpgs
  ref_cpgs[, chr := stringr::str_sub(chr, start = 4L)]
  ref_cpgs <- ref_cpgs[chr %in% all_chromosomes, .(chr, start)]
  setkey(ref_cpgs, chr, start)
  if (identical(update_cache, TRUE)) {
    write_cache_file(ref_cpgs, glue::glue("reference-{key}.fst"), dir)
  }
  return(ref_cpgs)
}

# From Methrix::extract_CPGs
#' Extracts all CpGs from a genome
#' @param ref_genome BSgenome object or name of the installed BSgenome package.
#'    Example: BSgenome.Hsapiens.UCSC.hg19
#' @param contigs A vector showing the contigs in the BSgenome object
#' @importFrom BSgenome installed.genomes getBSgenome
#' @return a list of data.table containing number of CpG's and contig lengths
extract_cpg_sites <- function(ref_genome = NULL,
  contigs = c(glue::glue("chr{1:22}"), "chrX", "chrY")) {
  pkgname <- chr <- NULL
  gnoms_installed <- BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = gnoms_installed)

  if (nrow(gnoms_installed) == 0) {
    logger::log_info("Could not find any installed BSgenomes.\n",
                  "Use BSgenome::available.genomes() for options.")
    stop()
  } else if (nrow(gnoms_installed[pkgname %in% ref_genome]) == 0) {
    logger::log_info("Could not find BSgenome ", ref_genome)
    logger::log_info("Found following BSgenome installations.",
                    "Provide the correct 'pkgname'.")
    print(gnoms_installed)
    stop()
  }

  requireNamespace(ref_genome, quietly = TRUE)
  ref_genome <- BSgenome::getBSgenome(genome = ref_genome)

  logger::log_info("=> Extracting CpGs")

  # Code borrwed from from: https://support.bioconductor.org/p/95239/
  cgs <- lapply(contigs, function(x) {
    start(Biostrings::matchPattern("CG", ref_genome[[x]]))
  })
  cpgs <- do.call(c, lapply(seq_along(contigs), function(x) {
    GenomicRanges::GRanges(names(ref_genome)[x], IRanges::IRanges(cgs[[x]],
      width = 2))
  }))

  cpgs <- data.table::as.data.table(as.data.frame(cpgs,
    stringsAsFactors = FALSE))

  colnames(cpgs) <- c("chr", "start", "end", "width", "strand")
  cpgs[, `:=`(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end),
    width = as.numeric(width)
  )]
  data.table::setkey(cpgs, "chr", "start")

  num_of_cpgs <- format(nrow(cpgs), big.mark = ",")
  logger::log_info("=> Extracted {num_of_cpgs} CpG Sites")
  return(list(cpgs = cpgs))
}
