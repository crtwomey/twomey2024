#
# vocabulary_likelihoods_check.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# computed by summing over each n-1 vocabulary's +1 modes
precursor_likelihoods <- function(lls) {
  lls |>
    dplyr::group_by(language, nminus1) |>
    dplyr::summarize(likelihood = sum(likelihood)) |>
    dplyr::mutate(l_rank = rank(-likelihood)) |>
    dplyr::ungroup()
}

# check the sensitivity of the likelihood calculations to particular choices of
# scale parameter
compute_likelihood_sensitivity <- function(
  fd_table,
  used,
  d_scales = likelihood_distance_scales
) {
  # previously computed precursor likelihoods
  used <- used |> precursor_likelihoods()

  # correlations (pearson and spearman for likelihoods and ranked likelihoods,
  # respectively) between used likelihoods and test likelihoods computed with a
  # particular choice of scale parameter
  correlation <- lapply(d_scales, function(s) {
    test <- fd_table |>
      likelihoods_from_soft_threshold(s) |>
      precursor_likelihoods()
    tidyr::tibble(
      l = cor(used$likelihood, test$likelihood, method = "pearson"),
      r = cor(used$l_rank, test$l_rank, method = "spearman")
    )
  })
  do.call("rbind", correlation)
}

# example of varying scale parameter for one language
compute_sensitivity_example <- function(
  fd_table,
  d_scales = likelihood_distance_scales,
  lnum     = 1
) {
  options(dplyr.summarise.inform = FALSE)

  lexample <- fd_table |> dplyr::filter(language == lnum)
  t(sapply(d_scales, function(s) {
    (lexample |>
      likelihoods_from_soft_threshold(s) |>
      precursor_likelihoods()
    )$likelihood
  }))
}

# plots for assessing how sensitive our likelihoods estimates are in comparison
# to particular choices of the scale parameter
plot_likelihood_sensitivities <- function(
  likelihood_sensitivity,
  sensitivity_example,
  lb_distance_scales,
  d_scales = likelihood_distance_scales
) {
  dbs <- range(d_scales)
  lbs <- range(lb_distance_scales)

  structure(function() {
    op  <- par(mar = c(5, 5, 3, 2), mfrow = c(1, 3))

    plot(NA, xlim = dbs,
      log  = "x", bty = "l", las = 1, xaxs = "i", yaxs = "i",
      ylim = c(0, 1),
      xlab = expression(sigma),
      ylab = "likelihood of precursor vocabulary",
      main = "example language"
    )
    abline(v = lbs, col = "gray60", lty = 2)
    abline(h = seq(0, 1, 0.2), col = "gray80", lty = 3, lwd = 1.0)
    abline(v = 10^seq(-5, 5, 2), col = "gray80", lty = 3, lwd = 1.0)
    for (i in seq_len(ncol(sensitivity_example))) {
      lines(d_scales, sensitivity_example[, i], col = i)
    }

    plot(d_scales, likelihood_sensitivity$l,
      type = "l", log = "x", las = 1, bty = "l", xaxs = "i", yaxs = "i",
      ylim = c(0, 1),
      xlab = expression(sigma),
      ylab = "Pearson correlation",
      main = "likelihood (all languages)"
    )
    abline(v = lbs, col = "gray60", lty = 2)
    abline(h = seq(0, 1, 0.2), col = "gray80", lty = 3, lwd = 1.0)
    abline(v = 10^seq(-5, 5, 2), col = "gray80", lty = 3, lwd = 1.0)

    plot(d_scales, likelihood_sensitivity$r,
      type = "l", log = "x", las = 1, bty = "l", xaxs = "i", yaxs = "i",
      ylim = c(0, 1),
      xlab = expression(sigma),
      ylab = "Spearman correlation",
      main = "rank likelihood (all languages)"
    )
    abline(v = lbs, col = "gray60", lty = 2)
    abline(h = seq(0, 1, 0.2), col = "gray80", lty = 3, lwd = 1.0)
    abline(v = 10^seq(-5, 5, 2), col = "gray80", lty = 3, lwd = 1.0)

    par(op)
  }, class = "base_plot")
}
