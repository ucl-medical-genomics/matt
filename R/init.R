#' @import data.table

.matt_env <- new.env(parent = emptyenv())

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching MATT version ",
                        utils::packageDescription("matt")$Version, ".")
  .matt_env[["cache_dir"]] <- getOption("matt.cache_directory",
    default = "./matt_cache_store")
}
