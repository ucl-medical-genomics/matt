
add_neighbouring_data_fast <- function(dt, grouping = c("dt_key", "chr")) {
  setkey(dt, chr, start)
  dt[!is.na(meth.y), start_idx := start]
  dt[, `:=`(
    b_meth = nafill(meth.y, type = "locf"),
    b_start = nafill(start_idx, type = "locf"),

    f_meth = nafill(meth.y, type = "nocb"),
    f_start = nafill(start_idx, type = "nocb")
  ), by = grouping]
  dt[, `:=`(
    b_dist = start - b_start,
    f_dist = f_start - start
  )]

  set(dt, NULL, "start_idx", NULL)
  return(dt)
}

generate_corr_matrix <- function(dt, grouping = c("dt_key", "chr")) {
  train_dt <- dt[to_impute == FALSE]
  train_dt[, `:=`(
    f_dist = data.table::shift(start, 1L, type = "lead") - start,
    f_meth = data.table::shift(meth.y, 1L, type = "lead")
  ), by = grouping]
  meth_corr <- train_dt[, .(dist = f_dist, meth_0 = meth.y, meth_1 = f_meth)]
  meth_corr <- meth_corr[!is.na(meth_0) & !is.na(meth_1)]
  meth_corr[, round_dist := fcase(dist <= 200, as.numeric(dist),
                                 dist <= 500, round_any(dist, 10),
                                 dist <= 2000, round_any(dist, 100))]

  corr_matrix <- meth_corr[!is.na(round_dist),
    .(corr = stats::cor(meth_0, meth_1)), by = round_dist]
  corr_matrix <- corr_matrix[!is.na(corr), .(dist = round_dist, corr)]
  return(corr_matrix)
}

annotate_training_data <- function(dt, grouping = c("dt_key", "chr")) {
    dt[to_impute == FALSE, `:=`(
        f_start = data.table::shift(start, 1L, type = "lead"),
        f_meth = data.table::shift(meth.y, 1L, type = "lead"),
        b_start = data.table::shift(start, 1L, type = "lag"),
        b_meth = data.table::shift(meth.y, 1L, type = "lag")
    ), by = grouping]
    dt[to_impute == FALSE, `:=`(
        b_dist = start - b_start,
        f_dist = f_start - start
    )]
}

annotate_corr_data <- function(dt){
  corr_matrix <- generate_corr_matrix(dt, grouping = c("chr"))
  dt[, `:=`(
    round_f_dist = fcase(f_dist <= 200, as.double(f_dist),
                         f_dist <= 500, round_any(f_dist, 10),
                         f_dist <= 2000, round_any(f_dist, 100)),
    round_b_dist = fcase(b_dist <= 200, as.double(b_dist),
                         b_dist <= 500, round_any(b_dist, 10),
                         b_dist <= 2000, round_any(b_dist, 100))
  )]
  dt <- merge(dt, corr_matrix, by.x = "round_f_dist", by.y = "dist",
              all.x = TRUE)
  dt[is.na(corr), corr := 0]
  setnames(dt, "corr", "up_corr")

  dt <- merge(dt, corr_matrix, by.x = "round_b_dist", by.y = "dist",
              all.x = TRUE)
  dt[is.na(corr), corr := 0]
  setnames(dt, "corr", "down_corr")
  set(dt, NULL, c("round_b_dist", "round_f_dist"), NULL)
  return(dt)
}


# add_neighbouring_info <- function(dt, num_neighbours = 2, meth_field = "meth.y",
#                                   update_original = TRUE,
#                                   additional_grouping = c()) {
#   j <- seq_len(num_neighbours)
#   if (identical(update_original, TRUE)) tmp <- dt
#   if (!identical(update_original, TRUE)) tmp <- copy(dt)
#   setorder(tmp, "chr", "start")

#   tmp[to_impute == FALSE, c(
#     glue::glue("down_{meth_field}_{j}"),
#     glue::glue("up_{meth_field}_{j}"),
#     glue::glue("down_start_{j}"),
#     glue::glue("up_start_{j}")
#   ) := c(
#     lapply(j, function(x) data.table::shift(get(meth_field), n = x)),
#     lapply(j, function(x) {
#       data.table::shift(get(meth_field), n = x, type = "lead")
#     }),
#     lapply(j, function(x) data.table::shift(start, n = x)),
#     lapply(j, function(x) data.table::shift(start, n = x, type = "lead"))
#   ), by = c("chr", additional_grouping)]

#   tmp[, c(
#     glue::glue("down_dist_{j}"),
#     glue::glue("up_dist_{j}")
#   ) := c(
#     lapply(j, function(x) get(glue::glue("down_start_{x}")) - start),
#     lapply(j, function(x) get(glue::glue("up_start_{x}")) - start)
#   )]

#   if (!identical(update_original, TRUE)) return(tmp)
# }

# # dt
# # num_neighbours = 1
# # meth_field = "meth.y"
# # additional_grouping = c("dt_key")
# annotate_with_closest_neighbours <- function(dt, num_neighbours = 1,
#                                              meth_field = "meth.y",
#                                              additional_grouping = c()) {
#   j <- seq_len(num_neighbours)
#   meth_fields <- c(glue::glue("down_{meth_field}_{j}"), glue::glue("up_{meth_field}_{j}"))
#   all_fields <- c(meth_fields, glue::glue("down_start_{j}"), glue::glue("up_start_{j}"),
#                   glue::glue("down_dist_{j}"), glue::glue("up_dist_{j}"))

#   lapply(meth_fields, function(field) {
#     dt[to_impute == FALSE & is.na(get(field)), to_impute := TRUE]
#     ""
#   })
#   dt[to_impute == TRUE, c(all_fields) := NA]

#   dt_ref <- dt[to_impute == FALSE,
#               c("chr", "start", ..additional_grouping, ..meth_field, ..all_fields)]
#   # add initial start to work out where the joining row is coming from
#   dt_ref[, nearest_start := start]
#   dt_ref[, ref_dt_key := dt_key]

#   setkey(dt_ref, dt_key, chr, start)
#   setkey(dt, dt_key, chr, start)

#   mdt <- dt_ref[dt, roll = "nearest"]
#   mdt_orig <- copy(mdt)
#   mdt <- copy(mdt_orig)

#   # mdt[, all_fields] are all from dt_ref (with repeated rows nearest)
#   setnames(mdt, all_fields, glue::glue("n.{all_fields}"))
#   # mdt[, i.{all_fields}] are all from dt and contain orig values for to_impute
#   setnames(mdt, glue::glue("i.{all_fields}"), all_fields)

#   # switch meth.y and i.meth.y
#   mdt[, c("meth.y", "n.meth.y") := .SD, .SDcols = c("i.meth.y", "meth.y")]
#   setcolorder(mdt, names(dt))

#   mdt[start < n.down_start_1,
#     c(glue::glue("n.{all_fields}"), "nearest_start") := NA]
#   mdt[start > n.up_start_1, c(glue::glue("n.{all_fields}"), "nearest_start") := NA]

#   mdt[start < nearest_start, c(
#     "up_start_1",
#     glue::glue("up_{meth_field}_1")
#   ) := c(
#     .(as.integer(nearest_start)),
#     .(get(glue::glue("n.{meth_field}")))
#   )]
#   mdt[start < nearest_start, c(
#     glue::glue("up_start_{tail(j, -1)}"),
#     glue::glue("up_{meth_field}_{tail(j, -1)}"),
#     glue::glue("down_start_{j}"),
#     glue::glue("down_{meth_field}_{j}")
#   ) := c(
#     lapply(head(j, -1), function(x) get(glue::glue("n.up_start_{x}"))),
#     lapply(head(j, -1), function(x) get(glue::glue("n.up_{meth_field}_{x}"))),
#     lapply(j, function(x) get(glue::glue("n.down_start_{x}"))),
#     lapply(j, function(x) get(glue::glue("n.down_{meth_field}_{x}")))
#   )]

#   # start > start.x -> When is up_stream of the closest neighbour
#   mdt[start > nearest_start, c(
#     "down_start_1",
#     glue::glue("down_{meth_field}_1")
#   ) := c(
#     .(as.integer(nearest_start)),
#     .(get(glue::glue("n.{meth_field}")))
#   )]
#   mdt[start > nearest_start, c(
#     glue::glue("up_start_{j}"),
#     glue::glue("up_{meth_field}_{j}"),
#     glue::glue("down_start_{tail(j, -1)}"),
#     glue::glue("down_{meth_field}_{tail(j, -1)}")
#   ) := c(
#     lapply(j, function(x) get(glue::glue("n.up_start_{x}"))),
#     lapply(j, function(x) get(glue::glue("n.up_{meth_field}_{x}"))),
#     lapply(head(j, -1), function(x) get(glue::glue("n.down_start_{x}"))),
#     lapply(head(j, -1), function(x) get(glue::glue("n.down_{meth_field}_{x}")))
#   )]

#   set(mdt, NULL, subset(names(mdt), startsWith(names(mdt), "n.")), NULL)

#   # recalculate all distances
#   mdt[, c(
#     glue::glue("down_dist_{j}"),
#     glue::glue("up_dist_{j}")
#   ) := c(
#     lapply(j, function(x) ifelse(get(glue::glue("down_start_{x}")) >= start, NA, get(glue::glue("down_start_{x}")) - start)),
#     lapply(j, function(x) ifelse(get(glue::glue("up_start_{x}")) <= start, NA, get(glue::glue("up_start_{x}")) - start))
#   )]

#   # Remove all cases where downstream is actually upstream or vice versa
#   lapply(j, function(x) {
#     mdt[is.na(get(glue::glue("up_dist_{x}"))) | is.na(get(glue::glue("down_dist_{x}"))), 
#       c(glue::glue("up_start_{x}"), glue::glue("up_{meth_field}_{x}"),
#         glue::glue("down_start_{x}"), glue::glue("down_{meth_field}_{x}"),
#         glue::glue("down_dist_{x}"), glue::glue("up_dist_{x}")
#         ) := NA]
#   })
#   return(mdt)
# }

# get_neighbors <- function(dt) {
#   .Call("_boostme_neighbors", PACKAGE = "boostme", dt)
# }

# getNeighbouringCpGsites <- function(orig_dt, update_original = TRUE) {
#   if (identical(update_original, TRUE)) dt <- orig_dt
#   if (!identical(update_original, TRUE)) dt <- copy(orig_dt)

#   neighb_in <- dt[, .(chr = chr, pos = start, meth = meth.y)]
#   neighb_in[chr == "X", chr := 23]
#   neighb_in[chr == "Y", chr := 24]
#   neighb_in[, chr := as.integer(chr)]
#   neighb <- get_neighbors(as.matrix(neighb_in))
#   dt[, upBeta := neighb$upstreamBetas]
#   dt[, downBeta := rev(neighb$downstreamBetas)]
#   dt[, upDist := neighb$upstreamDistances]
#   dt[, downDist := rev(neighb$downstreamDistances)]
#   return(dt)
# }


# getMultiNeighbouringCpGsites <- function(orig_dt, update_original = TRUE) {
#   if (identical(update_original, TRUE)) dt <- orig_dt
#   if (!identical(update_original, TRUE)) dt <- copy(orig_dt)

#   neighb_in <- dt[, .(chr = chrNum, pos = start, meth = meth.y)]
#   neighb <- get_neighbors(as.matrix(neighb_in))
#   dt[, upBeta := neighb$upstreamBetas]
#   dt[, downBeta := rev(neighb$downstreamBetas)]
#   dt[, upDist := neighb$upstreamDistances]
#   dt[, downDist := rev(neighb$downstreamDistances)]
#   return(dt)
# }
