#
# vocabulary_summary_plot.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# number of n+1 vocabulary modes as a function of vocabulary size
plot_nmodes_summary <- function(modes_summary) {
  modes_summary |>
    ggplot(aes(x = nwords, y = nmodes, color = type, fill = type)) +
      labs(title = "Number of (n + 1)-word vocabularies",
           x     = "number of words (n + 1)",
           y     = "number of vocabularies",
           col   = "Precursor",
           fill  = "Precursor") +
      scale_colour_discrete(labels = c("de novo", "historical")) +
      scale_fill_discrete(labels = c("de novo", "historical")) +
      scale_x_continuous(breaks = 3:max(modes_summary$nwords)) +
      scale_y_continuous(breaks = c(1, seq(10, 40, 10))) +
      geom_smooth(method = "loess") +
      geom_jitter(width = 0.1, height = 0, size = 0.5) +
      stat_summary(fun.data = "mean_cl_boot") +
      theme(panel.grid.minor.x = element_blank())
}

# effective number of n+1 vocabulary modes as a function of vocabulary size
plot_eff_nmodes_summary <- function(modes_summary) {
  modes_summary |>
    ggplot(aes(x = nwords, y = eff_nmodes, color = type, fill = type)) +
      labs(title = "Effective number of (n + 1)-word vocabularies",
           x     = "number of words (n + 1)",
           y     = expression("effective number of vocabularies, " * 2^{H(V)}),
           col   = "Precursor",
           fill  = "Precursor") +
      scale_colour_discrete(labels = c("de novo", "historical")) +
      scale_fill_discrete(labels = c("de novo", "historical")) +
      scale_x_continuous(breaks = 3:max(modes_summary$nwords)) +
      scale_y_continuous(breaks = seq(1, 10, 1)) +
      geom_smooth(method = "loess") +
      geom_jitter(width = 0.1, height = 0, size = 0.5) +
      stat_summary(fun.data = "mean_cl_boot") +
      theme(panel.grid.minor.x = element_blank())
}

compute_postprocessed_nplus1_pairwise_distances <- function(
  wcs_words,
  nplus1_min_distances,
  nplus1_pairwise_distances
) {
  min_distances <- purrr::map_dfr(wcs_languages, function(lnum) {
    modes <- nplus1_pairwise_distances[[lnum]]$terms |>
      dplyr::select(type, mode, size) |>
      dplyr::distinct()
    nplus1_min_distances[[lnum]] |>
      dplyr::left_join(
        modes |> dplyr::rename(size_i = size),
        by = c("ti" = "type", "i" = "mode")
      ) |>
      dplyr::left_join(
        modes |> dplyr::rename(size_j = size),
        by = c("tj" = "type", "j" = "mode")
      ) |>
      tibble::add_column(language = lnum, .before = 1)
  }) |>
    dplyr::mutate(size = (size_i / 330) * (size_j / 330))

  # n + 1 word languages
  wcs_words |>
    dplyr::group_by(language) |>
    dplyr::summarize(nwords = as.integer(head(nwords, 1) + 1)) |>
    dplyr::left_join(min_distances, by = "language") |>
    dplyr::mutate(
      type = factor(paste(ti, tj, sep = "_"), levels = c(
        "denovo_denovo", "historical_historical", "historical_denovo"
      )),
      nwords = forcats::as_factor(nwords)
    ) |>
    dplyr::filter(!(ti == tj & i == j)) |>
    dplyr::filter(type != "historical_denovo") |>
    dplyr::select(-ti, -tj)
}

compute_pairwise_distance_percent_difference <- function(
  postprocessed_nplus1_pairwise_distances
) {
  postprocessed_nplus1_pairwise_distances |>
    dplyr::group_by(nwords, type) |>
    dplyr::summarize(avg = mean(min_distance)) |>
    tidyr::pivot_wider(names_from = "type", values_from = "avg") |>
    dplyr::mutate(
      nwords             = strtoi(nwords),
      percent_difference = (denovo_denovo / historical_historical - 1) * 100
    )
}

# pairwise minimum distance between n+1 word modes
plot_nplus1_min_distances <- function(
  postprocessed_nplus1_pairwise_distances
) {
  postprocessed_nplus1_pairwise_distances |>
    ggplot2::ggplot(aes(x = nwords, y = min_distance, fill = type)) +
      ggplot2::geom_boxplot() +
      labs(title = "Distance between (n + 1)-word vocabularies",
        x        = "number of words (n + 1)",
        y        = "pairwise distance",
        fill     = "Pairwise between"
      ) +
      scale_colour_discrete(labels = c("de novo", "historical")) +
      scale_fill_discrete(labels = c("de novo", "historical"))
}

# rate-distortion tradeoff across all n+1 word vocabulary modes
plot_modes_efficiency <- function(
  modes_efficiency
) {
  modes_efficiency |>
    ggplot(aes(
      y      = mutual_information,
      x      = average_distortion,
      shape  = type,
      fill   = type,
      color  = type,
      weight = weight,
      alpha  = weight
    )) +
      scale_colour_discrete(labels = c("de novo", "historical")) +
      scale_fill_discrete(labels = c("de novo", "historical")) +
      scale_shape_discrete(labels = c("de novo", "historical")) +
      guides(
        alpha  = FALSE,
        weight = FALSE,
        size   = FALSE
      ) +
      labs(
        title  = "Rate-distortion tradeoff",
        x      = "average distortion (MSE)",
        y      = "mutual information (bits)",
        shape  = "Precursor",
        fill   = "Precursor",
        col    = "Precursor",
        alpha  = "Frequency",
        weight = "Frequency"
      ) +
      geom_point(size = 0.75) +
      geom_smooth(method = "loess") +
      geom_curve(
        aes(x = 2250, y = 2.75, xend = 1100, yend = 3.75),
        arrow = arrow(
          length = unit(0.03, "npc"),
          type   = "closed" # Describes arrow head (open or closed)
        ),
        colour    = "gray70",
        size      = 0.5,
        angle     = 90,
        curvature = -0.25
      ) +
      annotate(
        geom   = "text",
        label  = "number of words",
        x      = (2250 + 1100) / 2 - 100,
        y      = (2.75 + 3.75) / 2,
        colour = "gray50",
        size   = 3.5,
        hjust  = 0
      )
}

# combined plot of vocabulary mode summaries
plot_modes_summary <- function(
  nmodes_summary_figure,
  eff_nmodes_summary_figure,
  modes_efficiency_figure,
  nplus1_min_distances_figure
) {
  ((nmodes_summary_figure + eff_nmodes_summary_figure) +
    plot_layout(guides = "collect")
  ) /
  (modes_efficiency_figure + nplus1_min_distances_figure) +
    plot_annotation(tag_levels = "a") &
    theme(
      legend.position = "bottom",
      plot.tag        = element_text(face = "bold")
    )
}
