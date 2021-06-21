summarise_data <- function(dt, field = "value", round_data = TRUE,
                            round_by = 1, keep_field_name = TRUE,
                            by_c = c("sample", "name")) {
  if (identical(round_data, TRUE)) {
    dt[, glue("rounded.{field}") := round_any(get(field), round_by)]
    orig_field <- field
    field <- glue("rounded.{field}")
  }

  summarised_dt <- dt[, .N, by = c(field, by_c)]
  if (identical(round_data, TRUE) && identical(keep_field_name, TRUE)) {
    setnames(summarised_dt, field, orig_field)
  }
  return(summarised_dt)
}

compare_datasets_from_strings <- function(d1_string, d2_string,
                                          subset_sites = FALSE,
                                          from_cache = TRUE,
                                          cache_columns = NULL) {
  dt1_cache_key <- tolower(d1_string)
  dt2_cache_key <- tolower(d2_string)
  cached_fname <- get_compare_cache_fname(dt1_cache_key, dt2_cache_key)
  if (from_cache & file.exists(cached_fname)) {
    log_info("=> Reading from Cache: {cached_fname}")
    dt <- read_fst(cached_fname, as.data.table = TRUE, columns = cache_columns)
    return(dt)
  }
  dt1 <- get(glue("d_{dt1_cache_key}"))
  dt2 <- get(glue("d_{dt2_cache_key}"))
  if (!identical(subset_sites, FALSE)) {
    dt1 <- merge(subset_sites, dt1, all.x = TRUE, by = c("chr", "start"))
    dt2 <- merge(subset_sites, dt2, all.x = TRUE, by = c("chr", "start"))

    dt1_cache_key <- glue("sub_{dt1_cache_key}")
    dt2_cache_key <- glue("sub_{dt2_cache_key}")
  }
  dt <- compare_data(dt1, dt2, dt1_cache_key, dt2_cache_key,
                     from_cache = from_cache)
  return(dt)
}

compare_data <- function(x, y, x_cache_key, y_cache_key,
                         by_index = c("chr", "start"), from_cache = TRUE,
                         update_cache = TRUE, cache_columns = NULL) {
  cached_fname <- get_compare_cache_fname(x_cache_key, y_cache_key)
  if (from_cache & file.exists(cached_fname)) {
    log_info("=> Reading from Cache: {cached_fname}")
    dt <- read_fst(cached_fname, as.data.table = TRUE, columns = cache_columns)
    return(dt)
  }
  dt <- merge(x, y, all = TRUE, by = by_index)
  dt[, `:=`(
    size_x = x[, .N],
    size_y = y[, .N],
    size_xy = dt[, .N],
    size_fxy = dt[!(is.na(meth.x) | is.na(meth.y)), .N],
    diff = meth.x - meth.y
  )]
  dt[, `:=`(abs_diff = abs(diff), s_diff = diff ^ 2)]
  if (update_cache) {
    write_fst(dt, path = cached_fname)
  }
  return(dt)
}

generate_dt_summary <- function(dt_orig, x_field = "meth.x", y_field = "meth.y",
                                cov_filter = 10, extra_fields = c(),
                                missing_penalty = 0.5) {
  if (dt_orig[, max(get(x_field), na.rm = TRUE)] > 1 & missing_penalty < 1) {
    missing_penalty <- missing_penalty * 100
  }

  dt <- copy(dt_orig)
  dt_alt <- dt[!is.na(get(x_field)) & !is.na(get(y_field)) &
               cov.x >= cov_filter & cov.y >= cov_filter]
  corr <- cor(dt_alt[, get(x_field)], dt_alt[, get(y_field)])
  data <- dt[cov.x >= cov_filter & cov.y >= cov_filter, .(
    overlap = .N,
    mae = mean(abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(s_diff, na.rm = TRUE)),
    r = corr,
    r2 = corr ^ 2,
    size_fxy = dt[1, size_fxy],
    size_xy = dt[1, size_xy],
    size_x = dt[1, size_x],
    size_y = dt[1, size_y]
  )]
  data[, corrected_mae := {
    total_size = if (dt[is.na(diff), .N] == 0) size_fxy else size_xy
    uncommon_sites_x = fifelse(size_x > total_size, size_x - total_size, 0)
    uncommon_sites_y = fifelse(size_y > total_size, size_y - total_size, 0)
    total_uncommon_num = (uncommon_sites_x + uncommon_sites_y)
    corrected_mae = (mae * (total_uncommon_num / total_size))
    .(fifelse(corrected_mae == 0, mae, corrected_mae))
  }]

  dt_alt <- dt[is.na(diff), `:=`(
    diff = missing_penalty,
    abs_diff = missing_penalty,
    s_diff = missing_penalty * missing_penalty
  )]
  data_with_penalty <- dt_alt[, .(
    mae_with_penalty = mean(abs_diff, na.rm = TRUE),
    rmse_with_penalty = sqrt(mean(s_diff, na.rm = TRUE))
  )]

  data <- cbind(data, data_with_penalty, dt[1][, ..extra_fields])
  return(data)
}

generate_dt_summary_string <- function(dt, x_field = "meth.x",
                                       y_field = "meth.y", extra_fields = c(),
                                       missing_penalty = 0.5) {
  if ("dt_key" %in% names(dt)) {
    sum_dt <-  dt[,
      generate_dt_summary(.SD, x_field, y_field, extra_fields = extra_fields,
                          missing_penalty = missing_penalty),
      by = dt_key]
    data <- sum_dt[, lapply(.SD, mean), .SDcols = is.numeric]
  } else {
    data <- generate_dt_summary(dt, x_field, y_field,
                                extra_fields = extra_fields,
                                missing_penalty = missing_penalty)

  }
  return(
    glue("{overlap} CpG sites; R: {r}; R^2: {r2}; MAE: {mae} ; RMSE: {rmse}\n",
      "(With missing penalty: MAE: {p_mae}; RMSE: {p_rmse})",
      overlap = format(data[, overlap], big.mark = ","),
      r = round(data[, r], 3), r2 = round(data[, r2], 3),
      mae = round(data[, mae], 3), rmse = round(data[, rmse], 3),
      p_mae = round(data[, mae_with_penalty], 3),
      p_rmse = round(data[, rmse_with_penalty], 3)
    )
  )
}

generate_single_dt_summary_string <- function(dt, field = "meth",
                                              cov_field = "cov", filter = 0) {
  num <- format(dt[, .N], big.mark = ",")
  mean <- round(dt[, mean(get(field), na.rm = TRUE)], 3)
  sd <- round(dt[, sd(get(field), na.rm = TRUE)], 3)
  cov <- round(dt[, mean(get(cov_field), na.rm = TRUE)], 3)
  q <- round(dt[, quantile(get(field), na.rm = TRUE)], 3)
  str <- glue("Number of CpG sites: {num}; Mean: {mean}; SD: {sd}; ",
              "Quantiles: [{q[1]}, {q[2]}, {q[3]}, {q[4]}, {q[5]}]; ",
              "Mean CpG Depth: {cov}.")
  return(str)
}

# produce all combinations from a vector
generate_dt_permutations <- function(dt_keys) {
  perm <- as.data.table(t(combn(dt_keys, 2)))
  setnames(perm, c("V1", "V2"), c("d_name_1", "d_name_2"))
  return(perm)
}

get_mean_cov <- function(dataset_key, cov_filter = 10) {
  if (length(dataset_key) == 1) {
    dt <- get(glue("d_{tolower(dataset_key)}"))
    cov <- dt[cov >= cov_filter, round(mean(cov))]
  } else {
    cov <- round(mean(sapply(dataset_key, function (key) {
      dt <- get(glue("d_{tolower(key)}"))
      return(dt[cov >= cov_filter, mean(cov)])
    })))
  }

  return(glue("{cov}x"))
}
