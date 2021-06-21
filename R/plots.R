draw_all_heatmaps <- function(dt, fname_prefix, title = FALSE, subtitle = NULL,
                              xlab = "Dataset 1", ylab = "Dataset 2",
                              text_size = 4.5, width = 8, height = 8,
                              show_labels = TRUE,
                              legend_location = "top",
                              legend_justification = "right",
                              plots_dir = COMPARISONS_PLOTS_DIR) {
  title <- NULL
  fname <- glue("{fname_prefix}_mae")
  if (identical(title, FALSE)) title <- "Pairwise Comparisons: MAE"
  draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title, "MAE",
    text_field = "display_mae", text_size = text_size, subtitle = subtitle,
    width = width, height = height, show_labels = show_labels,
    legend_location = legend_location,
    legend_justification = legend_justification)

  if (!identical(dt$display_mae, dt$display_corrected_mae)) {
    fname <- glue("{fname_prefix}_corrected_mae")
    if (identical(title, FALSE)) title <- "Pairwise Comparisons: Corrected MAE"
    draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title,
      "Corrected MAE", fill_field = "corrected_mae",
      text_field = "display_corrected_mae", text_size = text_size,
      subtitle = subtitle, width = width, height = height,
      show_labels = show_labels, legend_location = legend_location,
      legend_justification = legend_justification)
  }
  if (!identical(dt$display_mae, dt$display_mae_with_penalty)) {
    fname <- glue("{fname_prefix}_mae_with_penalty")
    if (identical(title, FALSE)) {
      title <- "Pairwise Comparisons: MAE With Missing Penalty"
    }
    draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title,
      "MAE with\nMissing Penalty", fill_field = "mae_with_penalty",
      text_field = "display_mae_with_penalty", text_size = text_size,
      subtitle = subtitle, width = width, height = height,
      show_labels = show_labels, legend_location = legend_location,
      legend_justification = legend_justification)
  }
  if (!identical(dt$display_rmse, dt$display_rmse_with_penalty)) {
    fname <- glue("{fname_prefix}_rmse_with_penalty")
    if (identical(title, FALSE)) {
      title <- "Pairwise Comparisons: RMSE With Missing Penalty"
    }
    draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title,
      "RMSE with\nMissing Penalty", fill_field = "rmse_with_penalty",
      text_field = "display_rmse_with_penalty", text_size = text_size,
      subtitle = subtitle, width = width, height = height,
      show_labels = show_labels, legend_location = legend_location,
      legend_justification = legend_justification)
  }
  fname <- glue("{fname_prefix}_RMSE")
  if (identical(title, FALSE)) title <- "Pairwise Comparisons: RMSE"
  draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title, "RMSE",
    fill_field = "rmse", text_field = "display_rmse", text_size = text_size,
    subtitle = subtitle, width = width, height = height,
    show_labels = show_labels, legend_location = legend_location,
    legend_justification = legend_justification)
  fname <- glue("{fname_prefix}_R")
  if (identical(title, FALSE)) {
    title <- "Pairwise Comparisons: Pearson Correlation (R)"
  }
  draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title, "R",
    fill_field = "r", text_field = "display_r", rev_colour = TRUE,
    text_size = text_size, subtitle = subtitle, width = width, height = height,
    show_labels = show_labels, legend_location = legend_location,
    legend_justification = legend_justification
  )
  fname <- glue("{fname_prefix}_num")
  if (identical(title, FALSE)) {
    title <- "Pairwise Comparisons: Num of Overlapping CpG sites"
  }
  draw_comparison_heatmap(dt, plots_dir, fname, xlab, ylab, title,
    "Num of\nOverlapping\nCpG sites", fill_field = "overlap",
    text_field = "display_overlap", rev_colour = TRUE, si_labels = TRUE,
    text_size = text_size,
    subtitle = subtitle, width = width, height = height,
    show_labels = show_labels, legend_location = legend_location,
    legend_justification = legend_justification
  )
}

draw_comparison_heatmap <- function(dt, dir, fname, xlab, ylab, title,
                                    fill_label, x_field = "d_title_1",
                                    y_field = "d_title_2", subtitle = NULL,
                                    rev_colour = FALSE, fill_field = "mae",
                                    text_field = NULL, text_size = 2.5,
                                    reverse_y = FALSE, reverse_x = FALSE,
                                    show_labels = TRUE, si_labels = FALSE,
                                    legend_location = "right",
                                    legend_justification = "center",
                                    limits = NULL, width = 8, height = 8) {
  if (rev_colour) {
    colors <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(32)
  } else {
    colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(32)
  }
  if (is.null(text_field)) {
    text_field <- fill_field
  }
  g <- ggplot(dt, aes_string(x_field, y_field))
  g <- g + geom_raster(aes_string(fill = fill_field))
  if (identical(show_labels, TRUE)) {
    g <- g + geom_label(aes_string(label = text_field), size = text_size)
  }
  if (identical(si_labels, TRUE)) {
    g <- g + scale_fill_gradientn(colors = colors, limits = limits,
                                  label = label_number_si())
  } else {
    g <- g + scale_fill_gradientn(colors = colors, limits = limits)
  }
  g <- g + labs(x = xlab, y = ylab, title = title, fill = fill_label,
                subtitle = subtitle)
  g <- g + theme(legend.position = legend_location,
                 legend.justification = legend_justification,
                 axis.text = element_text(size = 12))

  g <- g + ggpubr::rotate_x_text()
  if (identical(reverse_y, TRUE)) g <- g + scale_y_discrete(limits = rev)
  if (identical(reverse_x, TRUE)) g <- g + scale_x_discrete(limits = rev)

  ggsave(glue("{dir}/hq/{fname}.pdf"), plot = g, device = "pdf", width = width,
         height = height)
  ggsave(glue("{dir}/{fname}.jpg"), plot = g, device = "jpg", dpi = 600,
         width = width, height = height)

  return(g)
}

draw_smooth_scatter <- function(dt_summary, dt_x, dt_y, fname, xlab, ylab,
                                x_field = "meth.x", y_field = "meth.y",
                                title = NULL, subtitle = NULL,
                                fill_label = "Frequency", trans = "log10",
                                interpolate = TRUE, width = 8, height = 8,
                                legend_max = NULL, separate_legend = FALSE,
                                base_font_size = 13, return_plot_object = FALSE,
                                tile_height = NULL,
                                dir = SMOOTHSCATTER_PLOTS_DIR) {
  if (is.null(legend_max)) {
    legend_max <- 10 ^ dt_summary[, ceiling(max(log10(N)))]
  }
  g <- ggplot(dt_summary)
  if (identical(tile_height, NULL)) {
    g <- g + geom_raster(aes(x = get(x_field), y = get(y_field), fill = N),
                          interpolate = interpolate)
  } else {
    g <- g + geom_tile(aes(x = get(x_field), y = get(y_field), fill = N),
                       height = tile_height)
  }
  g <- g + scale_fill_gradientn(colors = R_COLOURS, trans = trans,
                                label = label_number_si(),
                                limits = c(NA, legend_max))
  g <- g + scale_x_continuous(limits = c(-5, 105), expand = c(0.01, -2))
  g <- g + scale_y_continuous(limits = c(-5, 105), expand = c(0.01, -2))
  g <- g + labs(x = xlab, y = ylab, title = title, fill = fill_label,
                subtitle = subtitle)
  g <- g + theme(plot.subtitle = element_text(size = 14, colour = "#132B43"),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14))
  g <- g + geom_abline(slope = 1, intercept = 0, colour = "#132B43",
                       linetype = "dashed")

  if (identical(separate_legend, TRUE)) {
    g_legend <- cowplot::get_legend(g)
    ggsave(plot = g_legend, glue("{dir}/hq/{fname}-legend.pdf"), device = "pdf",
           width = 1.1, height = 2)
    ggsave(plot = g_legend, glue("{dir}/{fname}-legend.jpg"), device = "jpg",
           dpi = 600, width = 1.1, height = 2)

    g <- g + theme(legend.position = "none")
  }

  xg <- add_axis_bar_plot(g, dt_x, x_field, "N", "x")
  yg <- add_axis_bar_plot(g, dt_y, y_field, "N", "y", coord_flip = TRUE)
  g <- insert_xaxis_grob(g, xg, grid::unit(.1, "null"), position = "top")
  g <- insert_yaxis_grob(g, yg, grid::unit(.1, "null"), position = "right")

  if (identical(return_plot_object, TRUE)) return(g)

  ggsave(plot = g, glue("{dir}/hq/{fname}.pdf"), device = "pdf", width = width,
         height = height)
  ggsave(plot = g, glue("{dir}/{fname}.jpg"), device = "jpg", dpi = 600,
         width = width, height = height)
  return(glue("{dir}/{fname}.jpg"))
}

add_axis_bar_plot <- function(g, dt, x_field, y_field, axis, coord_flip = FALSE,
                            fill = "#132B43", alpha = 0.7) {
  ag <- axis_canvas(g, data = dt, axis = axis, coord_flip = coord_flip,
                     aes_string(x = x_field, y = y_field))
  ag <- ag + geom_bar(fill = fill, alpha = alpha, size = 0.1, stat = "identity")
  if (identical(coord_flip, TRUE)) {
    ag <- ag + coord_flip()
  }
  return(ag)
}

preprocess_smooth_scatter_data <- function(dt, x_field = "meth.x",
                                           y_field = "meth.y",
                                           round_data = FALSE, x_round_by = 0.1,
                                           y_round_by = 0.1) {
  if (identical(round_data, TRUE)) {
    dt[, glue("rounded.{x_field}") := round_any(get(x_field), x_round_by)]
    dt[, glue("rounded.{y_field}") := round_any(get(y_field), y_round_by)]
    dt_summary <- dt[, .N, by = c(glue("rounded.{c(x_field, y_field)}"))]
    setnames(dt_summary, c(glue("rounded.{c(x_field, y_field)}")),
             c(x_field, y_field))
  } else {
    dt_summary <- dt[, .N, by = c(x_field, y_field)]
  }
  return(dt_summary)
}

summarise_smoothscatter_marginal_data <- function(dt, field, round_data) {
  if (identical(round_data, TRUE)) {
    s_dt <- dt[, .N, by = c(glue("rounded.{field}"))]
    setnames(s_dt, glue("rounded.{field}"), field)
  } else {
    s_dt <- dt[, .N, by = field]
  }
  return(s_dt)
}

run_draw_smooth_scatter <- function(d1_string, d2_string, fname_prefix = "ss",
                                    x_field = "meth.x", y_field = "meth.y",
                                    cov_filter = 10, meth_display = "%",
                                    subset_sites = FALSE, round_data = TRUE,
                                    x_round_by = 5, y_round_by = 5,
                                    xlab_prefix = NULL, ylab_prefix = NULL,
                                    subtitle = NULL,
                                    show_coverage_in_label = TRUE,
                                    add_subtitle_summary = FALSE,
                                    legend_max = NULL, separate_legend = FALSE,
                                    dir = SMOOTHSCATTER_PLOTS_DIR,
                                    return_plot_object = FALSE,
                                    from_cache = TRUE) {
  if (length(d1_string) == 1) {
    dt <- compare_datasets_from_strings(d1_string, d2_string, subset_sites,
                                        from_cache)
  } else {
    dts <- lapply(seq(1, length(d1_string)), function(x) {
      compare_datasets_from_strings(d1_string[x], d2_string[x], subset_sites,
                                    from_cache)
    })
    names(dts) <- d1_string
    dt <- rbindlist(dts, idcol = "dt_key")
  }

  if (identical(meth_display, "%")) {
    dt[, `:=`(meth.x = meth.x * 100, meth.y = meth.y * 100)]
    dt[, `:=`(diff = meth.x - meth.y)]
    dt[, `:=`(abs_diff = abs(diff), s_diff = diff ^ 2)]
  }

  fdt <- dt[cov.x >= cov_filter & cov.y >= cov_filter]
  dt_summary <- preprocess_smooth_scatter_data(fdt, x_field, y_field,
                                               round_data, x_round_by,
                                               y_round_by)

  dt_x <- summarise_smoothscatter_marginal_data(fdt, x_field, round_data)
  dt_y <- summarise_smoothscatter_marginal_data(fdt, y_field, round_data)

  if (is.null(xlab_prefix)) xlab_prefix <- d1_string
  if (is.null(ylab_prefix)) ylab_prefix <- d2_string
  if (identical(show_coverage_in_label, TRUE)) {
    xlab_prefix <- glue("{xlab_prefix} ({get_mean_cov(d1_string)})")
    ylab_prefix <- glue("{ylab_prefix} ({get_mean_cov(d2_string)})")
  }

  xlab <- glue("{xlab_prefix} - Percent methylation")
  ylab <- glue("{ylab_prefix} - Percent methylation")
  fname <- glue("{fname_prefix}-{d1_string[1]}-vs-{d2_string[1]}")

  if (identical(add_subtitle_summary, TRUE)) {
    subtitle <- generate_dt_summary_string(dt, x_field, y_field)
  }

  out <- draw_smooth_scatter(dt_summary, dt_x, dt_y, fname, xlab, ylab, x_field,
                             y_field, subtitle = subtitle,
                             legend_max = legend_max,
                             separate_legend = separate_legend, dir = dir,
                             return_plot_object = return_plot_object)
  return(out)
}

draw_venn_diagram_plot <- function(venn_stats, fname, sample_titles) {
  v <- as.list(venn_stats$V1)
  names(v) <- venn_stats$variable
  opts <- list(category = sample_titles, lwd = 2, lty = "blank", sigdigs = 2,
    fill = brewer.pal(5, "Set1"), print.mode = c("raw", "percent"),
    ind = FALSE, margin = 0.06)
  g <- do.call(draw.quintuple.venn, c(v, opts))

  idx <- sapply(g, function(i) grepl("text", i$name))
  for (i in seq_along(g[idx])) {
    label_first_line <- strsplit(g[idx][[i]]$label, "\n")[[1]][1]
    label_second_line <- strsplit(g[idx][[i]]$label, "\n")[[1]][2]

    num <- as.integer(label_first_line)
    percent <- as.numeric(gsub("[\\(\\)\\%]", "", label_second_line))
    if (!is.na(num)) {
      if (i == 31) {
        # i.e. the center most field
        g[idx][[i]]$gp$cex <- 1.5
      } else if (percent <= 1) {
        g[idx][[i]]$gp$cex <- 0.75
      } else if (percent <= 5) {
        g[idx][[i]]$gp$cex <- 0.9
      }

      formated_num <- paste(format(num, big.mark = ",", scientific = FALSE),
        label_second_line, sep = "\n")
      g[idx][[i]]$label <- formated_num
    }
  }

  jpeg(glue("{fname}.jpg"), width = 1200, height = 1200, res = 600,
    pointsize = 2.7)
  grid.draw(g)
  dev.off()

  # pdf(file = glue("{fname}.pdf"), pointsize = 8.5)
  # grid.draw(g)
  # dev.off()
}

draw_venn_diagram <- function(dt, fname, cache_key = "venn_data",
                              sample_list = c("S3_R1", "S3_R2", "S3_R3",
                                "S3_R4", "S3_R5"),
                              sample_titles = c("Agilent (S3_R1)", "Roche (S3_R2)", "Illumina (S3_R3)",
                                "Diagenode (S3_R4)", "NuGen (S3_R5)"),
                              from_cache = TRUE, update_cache = TRUE) {
  cached_fname <- get_fname(cache_key, "cache", ".fst")
  if (file.exists(cached_fname) && from_cache) {
    log_info("=> Reading from Cache: {cached_fname}")
    venn_stats <- read_fst(cached_fname, as.data.table = TRUE)
  } else {
    venn_stats <- calculate_venn_params(dt, sample_list)
    if (update_cache) {
      write_fst(venn_stats, path = cached_fname)
    }
  }
  draw_venn_diagram_plot(venn_stats, fname, sample_titles)
}

calculate_venn_params <- function(dt, sample_list = c("S3_R1", "S3_R2", "S3_R3",
                                    "S3_R4", "S3_R5"),
                                  remove.char.from.sample = -4L) {
  dt <- dt[, .(chr, start, dt_key)]
  dt[, dt_key := str_sub(dt_key, end = remove.char.from.sample)]
  dt_cast <- dcast(dt, chr + start ~ dt_key)
  for (j in sample_list) {
    set(dt_cast, j = j, value = fifelse(dt_cast[[j]] == 0, "0", "1"))
  }
  dt_cast[, sum := do.call(paste0, .SD), .SDcols = sample_list]

  summary <- dt_cast[, .N, by = sum]

  summary[, `:=`(
    area1 = str_count(sum, pattern = "1\\d\\d\\d\\d"),
    area2 = str_count(sum, pattern = "\\d1\\d\\d\\d"),
    area3 = str_count(sum, pattern = "\\d\\d1\\d\\d"),
    area4 = str_count(sum, pattern = "\\d\\d\\d1\\d"),
    area5 = str_count(sum, pattern = "\\d\\d\\d\\d1"),
    n12 = str_count(sum, pattern = "11\\d\\d\\d"),
    n13 = str_count(sum, pattern = "1\\d1\\d\\d"),
    n14 = str_count(sum, pattern = "1\\d\\d1\\d"),
    n15 = str_count(sum, pattern = "1\\d\\d\\d1"),
    n23 = str_count(sum, pattern = "\\d11\\d\\d"),
    n24 = str_count(sum, pattern = "\\d1\\d1\\d"),
    n25 = str_count(sum, pattern = "\\d1\\d\\d1"),
    n34 = str_count(sum, pattern = "\\d\\d11\\d"),
    n35 = str_count(sum, pattern = "\\d\\d1\\d1"),
    n45 = str_count(sum, pattern = "\\d\\d\\d11"),
    n123 = str_count(sum, pattern = "111\\d\\d"),
    n124 = str_count(sum, pattern = "11\\d1\\d"),
    n125 = str_count(sum, pattern = "11\\d\\d1"),
    n134 = str_count(sum, pattern = "1\\d11\\d"),
    n135 = str_count(sum, pattern = "1\\d1\\d1"),
    n145 = str_count(sum, pattern = "1\\d\\d11"),
    n234 = str_count(sum, pattern = "\\d111\\d"),
    n235 = str_count(sum, pattern = "\\d11\\d1"),
    n245 = str_count(sum, pattern = "\\d1\\d11"),
    n345 = str_count(sum, pattern = "\\d\\d111"),
    n1234 = str_count(sum, pattern = "1111\\d"),
    n1235 = str_count(sum, pattern = "111\\d1"),
    n1245 = str_count(sum, pattern = "11\\d11"),
    n1345 = str_count(sum, pattern = "1\\d111"),
    n2345 = str_count(sum, pattern = "\\d1111"),
    n12345 = str_count(sum, pattern = "11111")
  )]

  ssummary <- melt(summary, id.vars = c("sum", "N"))
  ssummary <- ssummary[value == 1, sum(N), by = variable]
  return(ssummary)
}
