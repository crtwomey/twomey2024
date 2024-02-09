#
# vocabulary_generation.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# language-specific rate-distortion object
rd_vocabulary <- function(vocabulary, x, px, beta) {
  nx      <- nrow(x)
  pw_c    <- do.call("rbind", vocabulary$pxhat_x)
  foci    <- do.call("rbind", vocabulary$foci)
  nwords  <- nrow(pw_c)
  pw      <- rep(1 / nwords, nwords)
  list(
    nxhat   = nwords,
    beta    = beta,
    pxhat_x = t(pw_c),
    pxhat   = pw,
    xhat    = foci,
    x       = x,
    px      = px,
    d       = squared_euclidean_distance
  )
}

# de-novo rate-distortion object
denovo_rd_vocabulary <- function(nwords, x, px, beta) {
  init_rd(x,
    nxhat = nwords,
    beta  = beta,
    px    = px,
    fd    = squared_euclidean_distance
  )
}

# rate-distortion optimization
optimize_rd <- function(rd,
  max_iterations = max_ba_iterations,
  min_diff       = min_ba_iteration_diff,
  ...
) {
  avg_d <- Inf
  for (t in seq_len(max_iterations)) {
    rd      <- update_rd(rd, ...)
    avg_d_i <- average_distortion(rd)
    if (abs(avg_d_i - avg_d) <= min_diff) break
  }
  return(rd)
}

# rate-distortion optimal vocabularies
rd_optimize_vocabularies <- function(wcs, wcs_words) {
  # WCS chip positions in CIE Lab space
  x <- wcs$chips$Lab

  # language-specific inferred communicative needs, p(x)
  lang_pxs <- lapply(wcs_languages, function(lnum) {
    lang_px(wcs, lnum)
  })

  # language-specific inferred betas
  lang_betas <- map_dbl(wcs_languages, function(lnum) {
    lang_beta(wcs, lnum)
  })

  parallel_lapply(wcs_languages, function(lnum) {
    wcs_words |>
      dplyr::filter(language == lnum) |>
      rd_vocabulary(x, lang_pxs[[lnum]], lang_betas[lnum]) |>
      optimize_rd(update_xhat = FALSE)
  })
}

# take a list of RD solutions and return a list of unique modes
unique_modes <- function(rds, parallel = FALSE) {
  # differences between modes
  delta <- NULL
  if (parallel) {
    delta <- parallel_rb_lapply(seq_along(rds), function(i) {
      sapply(seq_along(rds), function(j) {
        best_mapped_delta_rd(rds[[i]], rds[[j]])
      })
    })
  } else {
    delta <- lapply(seq_along(rds), function(i) {
      sapply(seq_along(rds), function(j) {
        best_mapped_delta_rd(rds[[i]], rds[[j]])
      })
    })
    delta <- do.call("rbind", delta)
  }

  # group together identical modes
  adjacency <- (delta < same_mode_threshold) + 0
  agraph    <- igraph::graph.adjacency(adjacency)
  clu       <- igraph::components(agraph)
  grp       <- igraph::groups(clu)
  nmem      <- length(grp)

  relabel_by_size <- seq_len(nmem)
  relabel_by_size[order(clu$csize, decreasing = TRUE)] <- seq_len(nmem)
  mem <- relabel_by_size[clu$membership]

  # use first entry as representative of the mode
  lapply(1:nmem, function(k) {
    list(
      mode = rds[[head(which(mem == k), 1)]],
      size = sum(mem == k)
    )
  })
}


# for each language, generate a set of rate-distortion optimal n-1 word
# vocabularies corresponding to the original n word vocabulary with one word
# deleted.
generate_nminus1_vocabularies <- function(rd_optimal_vocabularies) {
  parallel_lapply(rd_optimal_vocabularies, function(vocabulary) {
    rd_nm1s <- vector("list", length = vocabulary$nxhat)
    for (i in seq_len(vocabulary$nxhat)) {
      rd_nm1         <- vocabulary
      rd_nm1$nxhat   <- rd_nm1$nxhat - 1
      rd_nm1$xhat    <- rd_nm1$xhat[-i, ]
      rd_nm1$pxhat_x <- rd_nm1$xhat[-i, ]
      rd_nm1s[[i]]   <- optimize_rd(rd_nm1)
    }
    unique_modes(rd_nm1s)
  })
}

# check what adding a centroid at a given location to an existing solution does
# based on RD dynamics.
add_xhat <- function(rd, xhat0, pxhat0) {
  rd_np1       <- rd
  rd_np1$nxhat <- rd$nxhat + 1
  rd_np1$xhat  <- rbind(rd$xhat, xhat0)
  rd_np1$pxhat <- c(rd$pxhat, pxhat0)
  rd_np1$pxhat <- rd_np1$pxhat / sum(rd_np1$pxhat)
  optimize_rd(rd_np1)
}

# generate all the unique n + 1 vocabularies for a list of vocabulary modes
unique_nplus1_vocabularies <- function(wcs, vocabulary_modes) {
  # WCS chip positions in CIE Lab space
  x  <- wcs$chips$Lab
  nx <- nrow(x)

  # introduce term with very low initial frequency
  pxhat0 <- initial_frequency_of_new_term
  np1s   <- vector("list", length(vocabulary_modes))

  for (i in seq_along(vocabulary_modes)) {
    rd   <- vocabulary_modes[[i]]$mode
    np1i <- parallel_lapply(seq_len(nx), function(j) {
      add_xhat(rd, x[j, ], pxhat0)
    })
    np1s[[i]] <- unique_modes(np1i, parallel = TRUE)
  }

  return(np1s)
}

# generate n + 1 vocabularies for each language, where each language
# is given by a list of one or more vocabulary modes.
generate_nplus1_vocabularies <- function(wcs, vocabularies) {
  pbapply::pblapply(vocabularies, function(vocabulary) {
    unique_nplus1_vocabularies(wcs, vocabulary)
  })
}

# generate n + 1 vocabularies for each language, where each language
# is associated with a single RD vocabulary of size n
generate_nplus1_rd_vocabularies <- function(wcs, vocabularies) {
  # WCS chip positions in CIE Lab space
  x  <- wcs$chips$Lab
  nx <- nrow(x)

  # introduce term with very low initial frequency
  pxhat0 <- initial_frequency_of_new_term

  pbapply::pblapply(vocabularies, function(rd) {
    np1s <- parallel_lapply(seq_len(nx), function(j) {
      add_xhat(rd, x[j, ], pxhat0)
    })
    unique_modes(np1s, parallel = TRUE)
  })
}

# generate de novo vocabularies with communicative needs and rate-distortion
# trade-off (beta) as the given set of empirical languages, but with n +
# delta_n words. (E.g. we use this for ndelta = -1, 0, and +1)
generate_denovo_vocabularies <- function(wcs, wcs_words, delta_n = 0) {
  # WCS chip positions in CIE Lab space
  x <- wcs$chips$Lab

  # language-specific inferred communicative needs, p(x)
  lang_pxs <- lapply(wcs_languages, function(lnum) {
    lang_px(wcs, lnum)
  })

  # language-specific inferred betas
  lang_betas <- map_dbl(wcs_languages, function(lnum) {
    lang_beta(wcs, lnum)
  })

  pbapply::pblapply(wcs_languages, function(lnum) {
    nwords <- (wcs_words |> dplyr::filter(language == lnum))$nwords[1]
    px     <- lang_pxs[[lnum]]
    beta   <- lang_betas[lnum]
    denovo <- parallel_lapply(seq_len(num_gen_denovo_modes), function(i) {
      denovo_rd_vocabulary(nwords + delta_n, x, px, beta) |>
        optimize_rd()
    })
    unique_modes(denovo, parallel = TRUE)
  })
}

compute_nplus1_variability <- function(
  nplus1_vocabularies,
  denovo_nplus1_vocabularies
) {
  historical <- purrr::map_dfr(seq_along(nplus1_vocabularies), function(lnum) {
    modes <- nplus1_vocabularies[[lnum]]
    tidyr::tibble(
      language  = lnum,
      precursor = 1,
      frequency = purrr::map_int(modes, function(m) m$size)
    )
  })

  denovo <- purrr::map_dfr(seq_along(denovo_nplus1_vocabularies),
    function(lnum) {
      precursors <- denovo_nplus1_vocabularies[[lnum]]
      purrr::map_dfr(seq_along(precursors), function(i) {
        modes <- precursors[[i]]
        tidyr::tibble(
          language  = lnum,
          precursor = i,
          frequency = purrr::map_int(modes, function(m) m$size)
        )
      })
    }
  )

  combined <- rbind(historical, denovo)
  combined |>
    dplyr::group_by(language, precursor) |>
    dplyr::filter(n() > 1) |>
    dplyr::summarize(
      sd_f = sd(frequency)
    ) |>
    dplyr::ungroup() |>
    dplyr::summarize(
      mean = mean(sd_f / 330),
      se   = sd(sd_f / 330) / sqrt(n())
    )
}
