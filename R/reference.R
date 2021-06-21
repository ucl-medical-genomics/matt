
load_ref_data <- function(key, ...) {
  ref_file <- system.file("extdata", glue::glue("reference-{key}.fst"),
    package = "matt")

  # check local cache
  if (!fs::file_exists(ref_file)) {
    ref_file <- fs::path(.matt_env[["cache_dir"]], "reference-{key}.fst")
  }

  # check for file again
  if (!fs::file_exists(ref_file)) {
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
  all_cpgs <- methrix::extract_CPGs(ref_genome = ref_genome)
  ref_cpgs <- all_cpgs$cpgs
  ref_cpgs[, chr := stringr::str_sub(chr, start = 4L)]
  ref_cpgs <- ref_cpgs[chr %in% all_chromosomes, .(chr, start)]
  setkey(ref_cpgs, chr, start)
  if (identical(update_cache, TRUE)) {
    write_cache_file(ref_cpgs, glue::glue("reference-{key}.fst"), dir)
  }
  return(ref_cpgs)
}
