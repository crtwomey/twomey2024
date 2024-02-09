#
# vocabulary_examples_plot.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# plots the extant n word language alongside the most likely and least likely
# ancestral languages
plot_nminus1_vocabulary <- function(
  lnum, wcs, wcs_words,
  rd_optimal_vocabularies,
  denovo_nminus1_vocabularies,
  denovo_nminus1plus1_likelihoods
) {
  point_size        <- 2.25
  axis_text_size    <- 5
  likelihood_format <- function(ll) {
    format(ll, digits = 2, nsmall = 2)
  }

  # language's nminus1 vocabularies ranked by aggregate likelihood
  lang_nminus1 <- get_language_likelihoods(
    lnum, denovo_nminus1plus1_likelihoods
  )

  # get the top two most likely precursor vocabularies
  best      <- lang_nminus1$nminus1[1]
  next_best <- lang_nminus1$nminus1[2]

  rd_empirical <- rd_optimal_vocabularies[[lnum]]
  p_empirical  <- plot_rd(
    wcs, rd_empirical,
    point_size     = point_size,
    axis_text_size = axis_text_size,
    label = paste0(
      "n = ", rd_empirical$nxhat, " word RD vocabulary"
    )
  )
  rd_ancestral <- denovo_nminus1_vocabularies[[lnum]][[best]]$mode
  p_nminus1_best <- plot_rd(
    wcs, rd_ancestral,
    point_size     = point_size,
    axis_text_size = axis_text_size,
    label = paste0(
      "n = ", rd_ancestral$nxhat, " word ancestral vocabulary ",
      "(p = ",
        likelihood_format(lang_nminus1$likelihood[1]),
      ")"
    )
  ) + theme(plot.margin = margin(t = 10, b = 10, r = 0, l = 0))
  rd_ancestral <- denovo_nminus1_vocabularies[[lnum]][[next_best]]$mode
  p_nminus1_next_best <- plot_rd(
    wcs, rd_ancestral,
    point_size     = point_size,
    axis_text_size = axis_text_size,
    label = paste0(
      "n = ", rd_ancestral$nxhat, " word ancestral vocabulary ",
      "(p = ",
        likelihood_format(lang_nminus1$likelihood[2]),
      ")"
    )
  )
  p_empirical / p_nminus1_best / p_nminus1_next_best
}

# show examples of inferred ancestral vocabularies
plot_nminus1_vocabularies <- function(
  wcs, wcs_words,
  rd_optimal_vocabularies,
  denovo_nminus1_vocabularies,
  denovo_nminus1plus1_likelihoods
) {
  kbmm_stages       <- file.path("ext", "KBMM_stages.png")
  kbmm_trajectories <- file.path("ext", "KBMM_trajectories.png")
  (
    (
      figpatch::fig(kbmm_stages) /
      figpatch::fig(kbmm_trajectories) +
      plot_layout(heights = c(2, 1))
    ) |
    wrap_elements(plot_nminus1_vocabulary(
        fig5_language_examples[1], wcs, wcs_words,
        rd_optimal_vocabularies,
        denovo_nminus1_vocabularies,
        denovo_nminus1plus1_likelihoods
    )) |
    wrap_elements(plot_nminus1_vocabulary(
        fig5_language_examples[2], wcs, wcs_words,
        rd_optimal_vocabularies,
        denovo_nminus1_vocabularies,
        denovo_nminus1plus1_likelihoods
    ))
  ) +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag    = element_text(face = "bold"),
      plot.margin = margin(t = 0, b = 0, l = 0, r = 0)
    )
}
