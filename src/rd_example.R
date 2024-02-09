#
# rd_example.R
#
# Rate-distortion with finite cardinality (RDFC) based on Banerjee et al. 2005
# All code is open source under the GPL v3 (see LICENSE)
#

# generate RDFC results for a simple example with 9 points in a grid
gen_example_rd_tradeoff <- function() {
  xs <- tidyr::expand_grid(
    x = seq(0, 1, 0.5),
    y = seq(0, 1, 0.5)
  )
  xs <- as.matrix(xs)
  nx <- nrow(xs)
  px <- rep(1 / nx, nx)

  # run BA updates
  run_rd <- function(nxhat, beta, tmax = 500) {
    rd <- init_dirichlet_rd(xs,
      nxhat = nxhat,
      beta  = beta,
      px    = px,
      xhat  = matrix(runif(2 * nxhat), nxhat, 2)
    )
    for (i in 1:tmax) {
      rd <- update_rd(rd)
    }
    return(rd)
  }

  # run BA from many initializations and take the best
  best_run_rd <- function(nreps, ...) {
    best <- Inf
    opt  <- NULL
    for (i in seq_len(nreps)) {
      rd <- run_rd(...)
      d  <- average_distortion(rd)
      if (d < best) {
        best <- d
        opt  <- rd
      }
    }
    return(opt)
  }

  # get the best RD results across a range of betas
  scan_best_rds <- function(nxhat, betas, nreps = 10) {
    if (nxhat == 1) {
      xhat <- colSums(px * xs)
      d    <- squared_euclidean_distance(xs, matrix(xhat, ncol = 2))
      return(tidyr::tibble(
        nxhat = nxhat,
        beta  = betas,
        D     = sum(px * d),
        R     = 0
      ))
    }
    parallel_rb_lapply(betas, function(beta) {
      rd <- best_run_rd(nreps, nxhat, beta)
      tidyr::tibble(
        nxhat = nxhat,
        beta  = beta,
        D     = average_distortion(rd),
        R     = mutual_information(rd)
      )
    })
  }

  # repeat results for a range of nxhats
  betas <- c(0, exp(seq(log(1), log(100), length.out = 40)))
  nxhats <- seq_len(nx)
  rds    <- list()
  for (nxhat in nxhats) {
    message("nxhat = ", nxhat)
    rds[[nxhat]] <- scan_best_rds(nxhat, betas, nreps = 50)
  }

  # collect results and make nxhat a factor
  all_rds <- do.call("rbind", rds) |>
    dplyr::mutate(nxhat = factor(nxhat, levels = 9:1))

  # remove dominated results
  nondominated_rds <- list()
  for (i in nxhats) {
    u <- all_rds |> filter(nxhat == i)
    dominated <- sapply(seq_len(nrow(u)), function(i) {
      any(u$D[u$beta < u$beta[i]] < u$D[i])
    })
    nondominated_rds[[i]] <- u[!dominated, ]
  }

  do.call("rbind", nondominated_rds)
}

# plot the example rate-distortiont tradeoff curves
plot_rd_example <- function(rd_example) {
  # re-create the example
  xs <- tidyr::expand_grid(
    x = seq(0, 1, 0.5),
    y = seq(0, 1, 0.5)
  )
  xs <- as.matrix(xs)
  nx <- nrow(xs)
  px <- rep(1 / nx, nx)

  max_nxhat <- max(as.integer(rd_example$nxhat))
  infeasible_region <- rbind(
      rd_example[rd_example$nxhat == max_nxhat, ],
      tibble(nxhat = max_nxhat, beta = 0, D = 0, R = 0)
    ) |>
    mutate(id = 1)

  rd_example |>
    ggplot(aes(x = D, y = R, color = nxhat)) +
      geom_polygon(data = infeasible_region, fill = "gray90", linetype = 0) +
      geom_line(lwd = 1) +
      scale_x_continuous(
        limits = c(0, max(max(rd_example$D), 0.5)), expand = c(0, 0)
      ) +
      scale_y_continuous(limits = c(0, max(log2(nx), 3.5)), expand = c(0, 0)) +
      scale_colour_viridis(discrete = TRUE, direction = -1) +
      geom_hline(aes(yintercept = log2(nx)), lwd = 0.25, linetype = "dashed") +
      annotate(
        "text", y = log2(9) + 0.1, x = 0.48, label = "$H(X)$", vjust = 0
      ) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position    = c(0.95, 0.35),
        legend.spacing.y   = unit(-0.15, "cm")
      ) +
      guides(color = guide_legend(
        title        = "$|\\widehat{X}|$",
        title.vjust  = 4,
        override.aes = list(fill = NA),
        byrow        = TRUE
      )) +
      labs(
        x = "$D$",
        y = "$R$"
      )
}
