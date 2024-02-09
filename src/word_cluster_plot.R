#
# word_cluster_plot.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# show the mappings for each word (cluster of terms)
plot_words <- function(wcs, word_clusters) {
  descriptions <- c(
    "red", "black", "yellow",
    "green-blue", "off-white", "brown",
    "blue", "white", "pink",
    "purple", "gray", "green",
    "orange", "red-orange", "light-brown"
  )

  map_dfr(seq_along(word_clusters$grp_num_langs), function(i) {
    wcs_color_chips(wcs,
      chip_colors_from_mapping(wcs,
        word_clusters$grp_avg_px_xhat[[i]],
        word_clusters$grp_avg_pxhat_x[[i]]
    )) |> mutate(
      description = descriptions[i],
      nlang       = paste0("N = ", word_clusters$grp_num_langs[i]),
      word        = i
    )
  }) |>
    plot_wcs_color_chips(point_size = 1.5, axis_text_size = 3) +
      theme(
        strip.background = element_blank(),
        strip.text.x     = element_blank()
      ) +
      facet_wrap(~ word, ncol = 3) +
      geom_text(aes(x = 41.5, y = 0.5, label = nlang),
        hjust = 1, vjust = 0.25, size = 3,
        fontface = "plain", family = "sans"
      ) +
      geom_text(aes(x = 40 / 2, y = -Inf, label = description),
        vjust = 0.25, size = 3,
        fontface = "plain", family = "sans"
      ) +
      coord_cartesian(clip = "off")
}

# histogram of vocabulary size for WCS languages
plot_distribution_of_nwords <- function(wcs_nword_frequencies) {
  wcs_nword_frequencies |>
    ggplot(aes(x = lang_nwords, y = nlangs)) +
      labs(title = "WCS languages",
        x        = "number of words",
        y        = "frequency"
      ) +
      geom_bar(stat = "identity", fill = "black") +
      scale_x_continuous(breaks = seq(
        min(wcs_nword_frequencies$lang_nwords),
        max(wcs_nword_frequencies$lang_nwords),
        1
      )) +
      theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
      )
}

# combine panels for full figure
plot_word_cluster_figure <- function(
  wcs_color_stimuli_plot,
  distribution_of_nwords_figure,
  word_plots
) {
  # panels
  p1 <- wcs_color_stimuli_plot + theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
  p2 <- word_plots + theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
  (
    (p1 | distribution_of_nwords_figure) +
    plot_layout(widths = c(2.5, 1))
  ) / p2 +
    plot_layout(heights = c(1, 3)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
}
