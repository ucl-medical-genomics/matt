
read_cache_file <- function (cache_id, use_cache, dir = .matt_env[["cache_dir"]]) {
  fst_path <- fs::path(dir, cache_id)
  if (file.exists(fst_path) && identical(use_cache, TRUE)) {
    logger::log_info("=> Reading from Cache: {fst_path}")
    dt <- fst::read_fst(fst_path, as.data.table = TRUE)
    return(dt)
  }
  return(NULL)
}

write_cache_file <- function (dt, cache_id, dir = .matt_env[["cache_dir"]],
  compress = 50) {
  if (! fs::dir_exists(dir)) {
    logger::log_info("=> Cache Directory not found. Creating directory: {dir}")
    fs::dir_create(dir)
  }
  fst_path <- fs::path(dir, cache_id)
  logger::log_info("=> Writing to Cache: {fst_path}")
  fst::write_fst(dt, fst_path, compress)
}
