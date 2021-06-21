
get_summary <- function(d1_string, d2_string, subset_sites = FALSE,
                        meth_display = "%", from_cache = TRUE) {
  dt <- compare_datasets_from_strings(d1_string, d2_string, subset_sites,
                                      from_cache)
  if (meth_display == "%") {
    # RMSE + R are scale invariant
    dt[, `:=`(
      meth.x = meth.x * 100,
      meth.y = meth.y * 100
    )]
    dt[, diff := meth.x - meth.y]
    dt[, `:=`(abs_diff = abs(diff), s_diff = diff ^ 2)]
  }
  summary <- generate_dt_summary(dt)
  print(summary)
  return(summary)
}

all_vs_all_comparisons <- function(samples_summary, cache_key, x_dt_keys,
                                   y_dt_keys, subset_sites = FALSE,
                                   meth_display = "%", from_cache = TRUE,
                                   x_sub_start = 0, y_sub_start = 0,
                                   x_title_cov = TRUE, y_title_cov = TRUE) {
  fst_file <- get_fname(cache_key, "cache", ".fst")

  if (file.exists(fst_file) && identical(from_cache, TRUE)) {
    log_info("=> Reading from Cache: {fst_file}")
    dt <- read_fst(fst_file, as.data.table = TRUE)
  } else {
    all_samples <- generate_dt_permutations(samples_summary$dt_key)
    log_info("=> Found Combinations : {all_samples[, .N]}")

    f_all_samples <- all_samples[d_name_1 %in% x_dt_keys &
      d_name_2 %in% y_dt_keys]
    log_info("=> Filtered Combinations : {f_all_samples[, .N]}")

    dt <- f_all_samples[,
      get_summary(d_name_1, d_name_2, subset_sites, meth_display = meth_display),
      by = .(d_name_1, d_name_2)
    ]

    fwrite(dt, file = glue("{CSV_DATA_DIR}/{cache_key}.csv"))
    write_fst(dt, fst_file)
  }

  # Add the reverse -Add info on samples B-A using data from samples A-B
  dt_rev <- copy(dt)
  setnames(dt_rev, c("d_name_1", "d_name_2"), c("d_name_2", "d_name_1"))
  dt <- rbindlist(list(original = dt, reverse = dt_rev), idcol = "origin",
                  use.names = TRUE )

  dt[, `:=`(
    display_mae = round(mae, 1),
    display_rmse = round(rmse, 1),
    display_r = signif(r, 2),
    display_corrected_mae = round(corrected_mae, 1),
    display_mae_with_penalty = round(mae_with_penalty, 1),
    display_rmse_with_penalty = round(rmse_with_penalty, 1),
    display_overlap = f2si(signif(overlap, 2))
  )]

  f_dt <- dt
  f_dt <- reformat_dt_axis_titles(dt, x_dt_keys, y_dt_keys, x_sub_start,
                                  y_sub_start, x_title_cov, y_title_cov)

  return(f_dt)
}

#subset to just the required fields
reformat_dt_axis_titles <- function(dt, x_dt_keys = FALSE, y_dt_keys = FALSE,
                                    x_sub_start = 0, y_sub_start = 0,
                                    x_title_cov = TRUE, y_title_cov = TRUE) {
  dt[, c("d_title_1", "d_title_2") := .(d_name_1, d_name_2)]

  if (!identical(x_dt_keys, FALSE)) {
    dt <- dt[d_title_1 %in% x_dt_keys]
    dt$d_title_1 <- factor(dt$d_title_1, levels = x_dt_keys)
    # Rename the levels to add coverage to the end
    levels(dt$d_title_1) <- sapply(x_dt_keys, function(x) {
      if (identical(x_title_cov, TRUE)) {
        glue("{str_sub(x, x_sub_start)}\n({get_mean_cov(x)})")
      } else {
        str_sub(x, x_sub_start)
      }
    })
  }

  if (!identical(y_dt_keys, FALSE)) {
    dt <- dt[d_title_2 %in% y_dt_keys]
    dt$d_title_2 <- factor(dt$d_title_2, levels = y_dt_keys)
    # Rename the levels to add coverage to the end
    levels(dt$d_title_2) <- sapply(y_dt_keys, function(x) {
      if (identical(y_title_cov, TRUE)) {
        glue("{str_sub(x, y_sub_start)}\n({get_mean_cov(x)})")
      } else {
        str_sub(x, y_sub_start)
      }
    })
  }

  return(dt)
}

# averag by coverage
get_average_comp <- function(dt_orig) {
  dt <- copy(dt_orig)
  dt[, `:=`(
    cov_1 = as.integer(str_match(d_title_1, "(\\d+)x")[, 2]),
    cov_2 = as.integer(str_match(d_title_2, "(\\d+)x")[, 2]),
    set_cov_1 = as.integer(str_match(d_name_1, "D(\\d+)_")[, 2]),
    set_cov_2 = as.integer(str_match(d_name_2, "D(\\d+)_")[, 2])
  )][is.na(set_cov_1), set_cov_1 := 100][is.na(set_cov_2), set_cov_2 := 100]

  dt_mean <- dt[,
    lapply(.SD, mean),
    by = .(set_cov_1, set_cov_2),
    .SDcols = is.numeric
  ]
  dt_mean <- dt_mean[, .SD, .SDcols = unique(names(dt_mean))]

  dt_mean[, `:=`(
    d_title_1 = glue("D{set_cov_1}\n({cov_1}x)", .envir = .SD),
    d_title_2 = glue("D{set_cov_2}\n({cov_2}x)", .envir = .SD)
  )]

  levels_1 <- dt_mean[order(set_cov_1), unique(d_title_1)]
  levels_2 <- dt_mean[order(set_cov_2), unique(d_title_2)]
  dt_mean$d_title_1 <- factor(dt_mean$d_title_1, levels = levels_1)
  dt_mean$d_title_2 <- factor(dt_mean$d_title_2, levels = levels_2)

  # reformat numbers
  dt_mean[, `:=`(
    display_mae = signif(mae, 3),
    display_rmse = round(rmse, 3),
    display_r = signif(r, 2),
    display_corrected_mae = signif(corrected_mae, 3),
    display_mae_with_penalty = signif(mae_with_penalty, 3),
    display_rmse_with_penalty = round(rmse_with_penalty, 3),
    display_overlap = f2si(signif(overlap, 2))
  )]

  return(dt_mean)
}
