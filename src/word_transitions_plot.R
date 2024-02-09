#
# word_transitions_plot.R
#
# Plots for showing word transitions between stages.
# All code is open source under the GPL v3 (see LICENSE)
#

# show conditional transition probabilities as a heatmap over a 2d grid
plot_conditional_transition_panel <- function(
  word_clusters,
  probabilities,
  nwords,
  sfc,
  zlims
) {
  scale_fill_continuous <- sfc
  color_order <- c(0, plotting_color_order)

  probabilities <- probabilities[, c(color_order[-1], 16)]
  probabilities <- probabilities[color_order[-1], ]

  keep_n0       <- apply(probabilities, 2, function(x) all(is.finite(x)))
  probabilities <- probabilities[, keep_n0]

  keep_n1       <- rowSums(probabilities) > 0
  keep_n1[]     <- TRUE
  probabilities <- probabilities[keep_n1, ]

  b1 <- tibble(
    Var1 = names(keep_n0)[keep_n0],
    Var2 = names(keep_n1)[keep_n1][1],
    color = c(word_clusters$grp_avg_colors[color_order[-1]], NA)[keep_n0],
    probability = 1
  )

  b2 <- tibble(
    Var1 = names(keep_n0)[keep_n0][1],
    Var2 = names(keep_n1)[keep_n1],
    color = word_clusters$grp_avg_colors[color_order[-1]][keep_n1],
    probability = 1
  )

  m <- as.data.frame.table(t(probabilities),
    responseName = "probability"
  )
  m |>
    ggplot(aes(x = Var1, y = Var2, fill = probability)) +
      geom_tile() +
      coord_flip(clip = "off") +
      scale_fill_continuous(na.value = "black", limits = zlims) +
      labs(
        y = paste0(nwords + 1, " word colors"),
        x = paste0(nwords, " word colors")
      ) +
      theme(
        axis.ticks.x       = element_blank(),
        axis.text.x        = element_blank(),
        axis.line.x        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.text.y        = element_blank(),
        axis.line.y        = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.margin = margin(
          t = 5, r = 10, b = 10, l = 5
        )
      ) +
      guides(fill = guide_colorbar(
        label.position = "bottom",
        title.position = "left",
        title.vjust    = 0.85
      )) +
      geom_point(data = b1, color = b1$color, size = 3, shape = 15,
        position = position_nudge(y = -1)
      ) +
      geom_point(data = b2, color = b2$color, size = 3, shape = 15,
        position = position_nudge(x = -1)
      ) +
      annotate("text", y = 0, x = max(as.numeric(m$Var1)) + 0.2,
        label = "+", size = 4
      )
}

# helper function for plotting matrix panels
plot_conditional_transition_matricies_panel <- function(
  word_clusters,
  conditional_transition_matricies,
  type_num,
  nwords_num,
  nwords
) {
  probabilities <- conditional_transition_matricies[[type_num]][[nwords_num]]
  rownames(probabilities) <- seq_len(nrow(probabilities))
  colnames(probabilities) <- seq_len(ncol(probabilities))

  plot_conditional_transition_panel(word_clusters, probabilities, nwords,
    sfc = function(...) {
      scale_fill_viridis_c(..., option = "viridis")
    },
    zlims = c(0, 1)
  )
}

# show conditional transition panels for transitions between the first four
# stages (3->4, 4->5, 5->6, and 6->7)
plot_conditional_transitions <- function(
  word_clusters,
  conditional_transition_matricies
) {
  plot_conditional_transition_matricies_panel(
    word_clusters, conditional_transition_matricies, 2, 1, 3
  ) /
  plot_conditional_transition_matricies_panel(
    word_clusters, conditional_transition_matricies, 2, 2, 4
  ) /
  plot_conditional_transition_matricies_panel(
    word_clusters, conditional_transition_matricies, 2, 3, 5
  ) /
  plot_conditional_transition_matricies_panel(
    word_clusters, conditional_transition_matricies, 2, 4, 6
  ) +
  plot_layout(guides = "collect", ncol = 1, nrow = 4,
    heights  = c(5, 10, 12, 16)
  ) & theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
}

# sankey diagram of transitions between vocabulary size stages
plot_sankey_transitions_panel <- function(
  word_clusters,
  n_to_nplus1_transitions,
  transitions_type,
  panel_title
) {
  # order for displaying color words
  color_order <- plotting_color_order

  n_start <- min(n_to_nplus1_transitions$nplus1words) - 1
  n_end   <- max(n_to_nplus1_transitions$nplus1words) - 1

  n_to_nplus1_transitions |>
    tidyr::pivot_longer(
      cols      = !c(type, nplus1words, word_n1),
      names_to  = "word_n0",
      values_to = "count"
    ) |>
    dplyr::filter(type == transitions_type) |>
    dplyr::filter(!is.na(count)) |>
    dplyr::group_by(type, nplus1words) |>
    dplyr::mutate(
      freq = count / sum(count)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(word_n0 != "added") |>
    dplyr::filter(nplus1words <= 8) |>
    dplyr::mutate(nplus1words = max(nplus1words) - nplus1words) |>
    dplyr::mutate(
      word_n1 = factor(word_n1, levels = color_order, ordered = TRUE),
      word_n0 = factor(word_n0, levels = color_order, ordered = TRUE)
    ) |>
    ggplot(aes(
      x         = nplus1words + 1,
      next_x    = nplus1words,
      node      = word_n0,
      next_node = word_n1,
      fill      = word_n0,
      value     = freq
    )) +
      scale_x_continuous(labels = 7:3, expand = c(0, 0)) +
      scale_fill_manual(
        values = word_clusters$grp_avg_colors[color_order]
      ) +
      geom_sankey(type = "sankey",
        flow.alpha = 0.5, node.color = 1, space = 0.1
      ) +
      theme(
        axis.ticks.x       = element_blank(),
        axis.text.x        = element_blank(),
        axis.line.x        = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(
          t = 0, r = 0, b = 0, l = 0
        )
      ) +
      labs(
        title = panel_title,
        x     = "number of words (n)",
        fill  = "word"
      ) +
      theme(
        plot.title      = element_text(hjust = 0.5),
        legend.position = "none"
      ) +
      coord_flip()
}

# sankey diagram of transitions between vocabulary size stages
plot_sankey_transitions <- function(
  word_clusters,
  n_to_nplus1_transitions
) {
  plot_sankey_transitions_panel(
    word_clusters, n_to_nplus1_transitions,
    transitions_type = "historical",
    panel_title      = NULL
  ) +
  plot_layout(guides = "collect")
}

# show transition matricies and sankey diagram side-by-side
plot_transitions_figure <- function(
  conditional_transitions_figure,
  sankey_transitions_figure
) {
  (
    wrap_elements(conditional_transitions_figure) +
    sankey_transitions_figure
  ) +
    plot_layout(widths = c(3, 6.5)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
}

# show comparison of distribution over n+1 historical and de novo modes
plot_modes_figures <- function(
  wcs,
  nplus1_vocabularies,
  denovo_nplus1_vocabularies,
  lnum = 1,
  dnum = 1
) {
  h_modes  <- nplus1_vocabularies[[lnum]]
  d_modes  <- denovo_nplus1_vocabularies[[lnum]][[dnum]]
  nh_modes <- length(h_modes)
  nd_modes <- length(d_modes)

  delta <- lapply(seq_len(nh_modes), function(i) {
    sapply(seq_len(nd_modes), function(j) {
      best_mapped_delta_rd(
        h_modes[[i]]$mode,
        d_modes[[j]]$mode
      )
    })
  })
  delta <- do.call("rbind", delta)

  # group together identical modes
  matched <- (delta < same_mode_threshold) + 0

  # check for double matches
  matched <- sapply(seq_len(ncol(matched)), function(i) {
    matched_i <- matched[, i]
    delta_i   <- delta[, i]
    if (sum(matched_i) > 1) {
      best       <- head(which(delta_i == min(delta_i[matched_i == 1])), 1)
      best_match <- rep(0, length(matched_i))
      best_match[best] <- 1
      best_match
    } else {
      matched_i
    }
  })

  h_sizes <- purrr::map_int(h_modes, function(m) m$size)
  d_sizes <- purrr::map_int(d_modes, function(m) m$size)

  matched_modes <- cbind(
    seq_len(nh_modes),
    matched %*% seq_len(nd_modes)
  )
  unmatched <- setdiff(seq_len(nd_modes), matched_modes[, 2])
  if (length(unmatched) > 0) {
    matched_modes <- rbind(matched_modes, cbind(
      rep(0, length(unmatched)),
      rev(unmatched)
    ))
  }
  matched_modes <- cbind(seq_len(nrow(matched_modes)), matched_modes)
  matched_modes[matched_modes == 0] <- NA
  colnames(matched_modes) <- c("row", "historical", "denovo")
  matched_modes <- tidyr::as_tibble(matched_modes)

  # mode ids
  matched_modes_ids <- as_tibble(matched_modes)
  matched_hids      <- matched_modes_ids$historical[
    !is.na(matched_modes_ids$historical)
  ]
  matched_dids <- matched_modes_ids$denovo[
    is.na(matched_modes_ids$historical)
  ]

  # create thumbnail directory if needed
  thumbnail_dir <- paste0(fig_dir, "/mode_thumbnails")
  if (!file.exists(thumbnail_dir)) {
    dir.create(thumbnail_dir)
  }

  mode_thumbnail_file <- function(row_id) {
    file.path(thumbnail_dir, paste0("mode_", row_id, ".png"))
  }

  plot_mode <- function(modes, mnum, row_id) {
    mode_file <- mode_thumbnail_file(row_id)
    pl <- plot_rd(wcs,
      modes[[mnum]]$mode,
      point_size     = 0.55,
      axis_text_size = 0,
      boundry_lwd    = 0.1
    ) + theme(plot.margin = unit(c(0, 0, 0, 0), "null"))
    ggsave(pl, filename = mode_file,
      width = 5, height = 1, units = "cm", dpi = 300
    )
  }

  #-------------------------------------------------------
  # NOTE: very slow!
  for (row_id in seq_along(matched_hids)) {
    plot_mode(h_modes, matched_hids[row_id], row_id)
  }
  for (row_id in seq_along(matched_dids)) {
    plot_mode(d_modes, matched_dids[row_id], row_id + length(matched_hids))
  }
  #-------------------------------------------------------

  row_labels <- sapply(matched_modes_ids$row, function(row_id) {
    paste0("<img src='", mode_thumbnail_file(row_id), "' width='175' />")
  })

  # frequencies
  matched_modes$historical <- h_sizes[matched_modes$historical]
  matched_modes$denovo     <- d_sizes[matched_modes$denovo]
  matched_modes[is.na(matched_modes)] <- 0

  # reorder according to difference in frequency
  row_labels <- row_labels[
    order(matched_modes$historical - matched_modes$denovo)
  ]

  matched_modes |>
    dplyr::mutate(row = factor(row)) |>
    dplyr::mutate(row = fct_reorder(row, historical - denovo)) |>
    tidyr::pivot_longer(
      c("historical", "denovo"), names_to = "type", values_to = "frequency"
    ) |>
    ggplot(aes(x = row, y = frequency, fill = type)) +
      geom_col(position = position_dodge(), width = 0.75) +
      scale_x_discrete(name = NULL, labels = row_labels) + #rev(row_labels)) +
      scale_fill_discrete(labels = c("de novo", "historical")) +
      coord_flip() +
      theme(
        axis.ticks.length.y = unit(0, "cm"),
        axis.text.y = element_markdown(vjust = 0.5, margin = unit(-0.25, "cm")),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.justification = c(1, 0),
        legend.position      = c(1, 0.25),
        plot.margin          = margin(t = 0, l = 0, b = 0, r = 20)
      ) +
      geom_vline(
        xintercept = seq_len(nrow(matched_modes)) + 0.5,
        color     = "gray90"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(
        y    = "frequency",
        x    = "n + 1 vocabulary",
        fill = "Precursor"
      )
}

# show an example comparison between a WCS language's distribution over likely
# n+1 vocabularies and a related de novo n+1 distribution
plot_historical_denovo_comparison <- function(
  wcs,
  rd_optimal_vocabularies,
  denovo_vocabularies,
  nplus1_vocabularies,
  denovo_nplus1_vocabularies,
  word_clusters,
  n_to_nplus1_transitions
) {
  # example language and de novo comparison
  lnum   <- 11
  dnum   <- 1
  h_mode <- rd_optimal_vocabularies[[lnum]]
  d_mode <- denovo_vocabularies[[lnum]][[dnum]]$mode
  wrap_elements(plot_rd(wcs, h_mode,
    point_size     = 2.5,
    axis_text_size = 5,
    label = bquote(
      "historical " * .(h_mode$nxhat) * " word vocabulary"
    )
  ) + plot_rd(wcs, d_mode,
    point_size     = 2.5,
    axis_text_size = 5,
    label = bquote(
      "de novo " * .(d_mode$nxhat) *
      " word vocabulary with the same " * p(x) * " and " * beta
    )
  )) / wrap_elements(plot_modes_figures(
    wcs,
    nplus1_vocabularies,
    denovo_nplus1_vocabularies,
    lnum, dnum
  )) + plot_layout(heights = c(2.25, 4)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"),
      plot.margin  = margin(t = 0, r = 0, b = 0, l = 0)
    )
}
