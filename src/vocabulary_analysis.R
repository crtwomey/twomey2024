#
# vocabulary_analysis.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# find the minimum distance assignments of all the n-1+1 de novo vocabularies
# to the original empirical vocabulary.
compute_denovo_nminus1plus1_distances <- function(
  wcs_words,
  denovo_nminus1plus1_vocabularies,
  wcs_distance
) {
  pbapply::pblapply(wcs_languages, function(lnum) {
    vocab    <- wcs_words |> dplyr::filter(language == lnum)
    nminus1s <- denovo_nminus1plus1_vocabularies[[lnum]]
    denovo   <- purrr::map_dfr(seq_along(nminus1s), function(nminus1) {
      nplus1s <- nminus1s[[nminus1]]
      purrr::map_dfr(seq_along(nplus1s), function(nplus1) {
        rd     <- nplus1s[[nplus1]]$mode
        pxxhat <- t(rd$pxhat_x) * rd$px
        list(
          nminus1 = nminus1,
          nplus1  = nplus1,
          term    = seq_len(ncol(pxxhat)),
          px_xhat = lapply(seq_len(ncol(pxxhat)), function(i) {
            pxxhat[, i] / sum(pxxhat[, i])
          })
        )
      })
    })

    distances <- parallel_lapply(seq_len(nrow(denovo)), function(i) {
      purrr::map_dbl(seq_len(nrow(vocab)), function(j) {
        wcs_distance(denovo$px_xhat[[i]], vocab$px_xhat[[j]])
      })
    })
    distances <- cbind(denovo[, 1:2], do.call("rbind", distances))

    indexes     <- dplyr::distinct(denovo[, 1:2])
    assignments <- apply(indexes, 1, function(ijs) {
      cost <- as.matrix(distances |>
        dplyr::filter(nminus1 == ijs[1] & nplus1 == ijs[2]) |>
        dplyr::select(-nminus1, -nplus1))
      lp.assign(cost)
    })

    tibble::tibble(indexes, assignments)
  })
}

# find the minimum distance assignments of all the n-1+1 de novo vocabularies
# to each other (pairwise comparisons within a language).
compute_denovo_nminus1plus1_pairwise_distances <- function(
  denovo_nminus1plus1_vocabularies,
  wcs_distance
) {
  pbapply::pblapply(wcs_languages, function(lnum) {
    nminus1s <- denovo_nminus1plus1_vocabularies[[lnum]]
    denovo   <- purrr::map_dfr(seq_along(nminus1s), function(nminus1) {
      nplus1s <- nminus1s[[nminus1]]
      purrr::map_dfr(seq_along(nplus1s), function(nplus1) {
        rd     <- nplus1s[[nplus1]]$mode
        pxxhat <- t(rd$pxhat_x) * rd$px
        list(
          nminus1 = nminus1,
          nplus1  = nplus1,
          term    = seq_len(ncol(pxxhat)),
          px_xhat = lapply(seq_len(ncol(pxxhat)), function(i) {
            pxxhat[, i] / sum(pxxhat[, i])
          })
        )
      })
    })

    # compute upper triangle of distance matrix
    n   <- nrow(denovo)
    ijs <- lapply(1:n, function(i) cbind(i, seq_len(i)))
    ijs <- do.call("rbind", ijs)

    term_distances <- parallel_lapply(seq_len(nrow(ijs)), function(k) {
      i <- ijs[k, 1]
      j <- ijs[k, 2]
      wcs_distance(denovo$px_xhat[[i]], denovo$px_xhat[[j]])
    })

    distances <- do.call("c", term_distances)
    distance  <- matrix(Inf, n, n)
    distance[upper.tri(distance, diag = TRUE)] <- distances
    distance[lower.tri(distance)] <- t(distance)[lower.tri(distance)]

    # compute upper triangle of min distance between modes
    denovo_mode_ids <- denovo |>
      dplyr::select(nminus1, nplus1) |>
      dplyr::distinct()

    # compute upper triangle of min_distance matrix
    n   <- nrow(denovo_mode_ids)
    ijs <- lapply(1:n, function(i) cbind(i, seq_len(i)))
    ijs <- do.call("rbind", ijs)

    mode_distances <- parallel_lapply(seq_len(nrow(ijs)), function(k) {
      i <- ijs[k, 1]
      j <- ijs[k, 2]
      iis <- which(
        denovo$nminus1 == denovo_mode_ids$nminus1[i] &
        denovo$nplus1  == denovo_mode_ids$nplus1[i]
      )
      jjs <- which(
        denovo$nminus1 == denovo_mode_ids$nminus1[j] &
        denovo$nplus1  == denovo_mode_ids$nplus1[j]
      )
      lp.assign(distance[iis, jjs])$objval
    })

    min_distances <- do.call("c", mode_distances)
    min_distance  <- matrix(Inf, n, n)
    min_distance[
      upper.tri(min_distance, diag = TRUE)
    ] <- min_distances
    min_distance[
      lower.tri(min_distance)
    ] <- t(min_distance)[lower.tri(min_distance)]
    min_distance
  })
}

# builds a table of distances and frequencies for the provided vocabularies
build_frequency_and_distance_table <- function(
  denovo_nminus1plus1_distances,
  denovo_nminus1plus1_vocabularies
) {
  purrr::map_dfr(wcs_languages, function(lnum) {
    distances    <- denovo_nminus1plus1_distances[[lnum]]
    vocabularies <- denovo_nminus1plus1_vocabularies[[lnum]]
    frequencies  <- unlist(lapply(vocabularies, function(nminus1) {
      sapply(nminus1, function(m) m$size)
    }))
    distances |>
      tibble::add_column(language = lnum, .before = 1) |>
      tibble::add_column(
        distance  = purrr::map_dbl(distances$assignments, function(a) a$objval),
        frequency = frequencies
      )
  })
}

# computes weighted softmax avoiding numerical under/overflow
weighted_softmax <- function(ws, xs) {
  ys <- xs + log(ws)
  a  <- max(ys)
  exp(ys - log(sum(exp(ys - a))) - a)
}

# compute x * log(x) avoiding numerical issues near x = 0
safe_xlogx <- function(x) ifelse(dplyr::near(x, 0), 0, x * log(x))

# assign likelihoods based on frequency, distance, and a given distance scale
likelihoods_from_soft_threshold <- function(fds, s) {
  # likelihoods based on relative frequency and distance
  fds |>
    dplyr::group_by(language, nminus1) |>
    dplyr::mutate(total = sum(frequency)) |>
    dplyr::group_by(language) |>
    dplyr::mutate(likelihood = weighted_softmax(
      frequency / total, -s * distance
    )) |>
    dplyr::ungroup() |>
    dplyr::select(-total)
}

# determine active range of scale paramter
compute_scale_parameter_prior_bounds <- function(
  fds,
  d_scales         = likelihood_distance_scales,
  change_threshold = likelihood_change_threshold
) {
  # likelihood as a function of the scale parameter
  likelihoods <- t(sapply(d_scales, function(s) {
    (fds |>
      likelihoods_from_soft_threshold(s) |>
      dplyr::group_by(language, nminus1) |>
      dplyr::summarize(likelihood = sum(likelihood)) |>
      dplyr::ungroup()
    )$likelihood
  }))

  # language and n-1 vocabulary information
  llnums <- fds |>
    dplyr::group_by(language, nminus1) |>
    dplyr::summarize(nmodes = n()) |>
    dplyr::ungroup()

  # range based on the absolute value of the change in entropy as the scale
  # parameter increase. estimate bounds for each language based on magnitude
  # of change exceeding a small relative threshold (relative to the maximum
  # change for the language)
  ll_bounds <- t(sapply(seq_len(max(llnums$language)), function(lnum) {
    lang <- llnums$language == lnum
    if (sum(lang) <= 1) return(c(NA, NA))
    lls    <- apply(likelihoods[, lang], 1, function(p) -sum(safe_xlogx(p)))
    deltas <- abs(diff(lls))
    deltas <-  deltas / max(deltas)
    lb     <- range(which(deltas > change_threshold))
    (d_scales[-1])[lb]
  }))

  # return a range between the median lower and upper bounds across languages
  lbs <- apply(ll_bounds, 2, function(x) median(x, na.rm = TRUE))
  10^seq(log10(lbs[1]), log10(lbs[2]), length.out = 1000)
}

# computes weighted softmax avoiding numerical under/overflow
weighted_softmax2 <- function(w, s, d) {
  y <- -s * d + log(w)
  a <- max(y)
  exp(-s * d + log(w * d) - log(sum(exp(y - a))) - a)
}

fisher_information <- function(s, m, w, d) {
  ii <- unique(m)

  # likelihood of each precursor
  qs <- weighted_softmax(w, -s * d)
  ps <- sapply(ii, function(i) sum(qs[m == i]))

  z <- sum(weighted_softmax2(w, s, d))
  a <- sapply(ii, function(i) {
    sum(weighted_softmax2(w[m == i], s, d[m == i]))
  })

  sum(ps * (z - a)^2)
}

# approximate the likelihood of an n-1 vocabulary being the ancestor of a
# vocabulary based on the distance of +1 versions of the n-1 vocabulary to the
# focal n word vocabulary.
compute_precursor_likelihoods <- function(
  fd_table,
  d_scales = likelihood_distance_scales
) {

  # relevant language information
  fds <- fd_table |>
    dplyr::select(language, nminus1, frequency, distance)

  # compute likelihoods for every choice of scale parameter
  likelihoods <- t(sapply(d_scales, function(s) {
    (fds |>
      likelihoods_from_soft_threshold(s)
    )$likelihood
  }))

  # uniform prior over distance scale parameter
  ns <- length(d_scales)
  scale_prior <- rep(1 / ns, ns)

  # marginalize estimates over scale parameter
  fd_table |>
    dplyr::mutate(
      likelihood = colSums(likelihoods * scale_prior)
    ) |>
    dplyr::group_by(language, nminus1) |>
    dplyr::summarize(likelihood = sum(likelihood))
}

# convenience for getting a specific language's precursor likelihoods in
# descending order (most to least likely)
get_language_likelihoods <- function(
  lnum, likelihoods
) {
  likelihoods |>
    dplyr::filter(language == lnum) |>
    dplyr::group_by(nminus1) |>
    dplyr::arrange(desc(likelihood))
}

# compute distances between all pairs of n+1 vocabularies
# (historical and de novo).
compute_nplus1_pairwise_distances <- function(
  nplus1_vocabularies,
  nplus1_denovo_vocabularies,
  wcs_distance
) {
  compute_px_xhats <- function(type, m, rd, size) {
    pxxhat <- t(rd$pxhat_x) * rd$px
    list(
      type    = type,
      mode    = m,
      size    = size,
      term    = seq_len(ncol(pxxhat)),
      px_xhat = lapply(seq_len(ncol(pxxhat)), function(i) {
        pxxhat[, i] / sum(pxxhat[, i])
      })
    )
  }

  pbapply::pblapply(wcs_languages, function(lnum) {
    historical <- nplus1_vocabularies[[lnum]]
    denovo     <- nplus1_denovo_vocabularies[[lnum]]

    historical_terms <- purrr::map_dfr(seq_along(historical), function(i) {
      compute_px_xhats(
        "historical", i, historical[[i]]$mode, historical[[i]]$size
      )
    })

    denovo_terms <- purrr::map_dfr(seq_along(denovo), function(i) {
      compute_px_xhats(
        "denovo", i, denovo[[i]]$mode, denovo[[i]]$size
      )
    })

    all_terms <- rbind(historical_terms, denovo_terms)

    # compute upper triangle of distance matrix
    n   <- nrow(all_terms)
    ijs <- lapply(seq_len(n), function(i) cbind(i, seq_len(i)))
    ijs <- do.call("rbind", ijs)

    term_distances <- parallel_lapply(seq_len(nrow(ijs)), function(k) {
      i <- ijs[k, 1]
      j <- ijs[k, 2]
      wcs_distance(
        all_terms$px_xhat[[i]],
        all_terms$px_xhat[[j]]
      )
    })
    deltas <- do.call("c", term_distances)
    n      <- nrow(all_terms)
    delta  <- matrix(Inf, n, n)
    delta[upper.tri(delta, diag = TRUE)] <- deltas
    delta[lower.tri(delta)] <- t(delta)[lower.tri(delta)]
    list(
      terms     = all_terms,
      distances = delta
    )
  })
}

# compute pairwise distances between n+1 vocabularies based on the minimum
# earth mover's distance distance assignment between them.
compute_nplus1_min_distances <- function(
  nplus1_pairwise_distances
) {
  pbapply::pblapply(wcs_languages, function(lnum) {
    pairwise_distances <- nplus1_pairwise_distances[[lnum]]

    hmodes <- unique(
      (pairwise_distances$terms |>
        dplyr::filter(type == "historical"))$mode
    )
    dmodes <- unique(
      (pairwise_distances$terms |>
        dplyr::filter(type == "denovo"))$mode
    )

    ijs <- rbind(
      tidyr::expand_grid(i = hmodes, j = hmodes) |>
        tibble::add_column(ti = "historical", tj = "historical", .before = 1),
      tidyr::expand_grid(i = hmodes, j = dmodes) |>
        tibble::add_column(ti = "historical", tj = "denovo", .before = 1),
      tidyr::expand_grid(i = dmodes, j = dmodes) |>
        tibble::add_column(ti = "denovo", tj = "denovo", .before = 1)
    )

    min_distances <- parallel_lapply(seq_len(nrow(ijs)), function(k) {
        iis <- which(
          pairwise_distances$term$type == ijs$ti[k] &
          pairwise_distances$term$mode == ijs$i[k]
        )
        jjs <- which(
          pairwise_distances$term$type == ijs$tj[k] &
          pairwise_distances$term$mode == ijs$j[k]
        )
        cost <- pairwise_distances$distances[iis, jjs]
        lpSolve::lp.assign(cost)
    })

    ijs |> tibble::add_column(
      min_distance = purrr::map_dbl(min_distances, function(d) d$objval)
    )
  })
}

# compute absolute and effective number of modes
vocabulary_modes_summary <- function(vocabularies) {
  f_nwords <- function(modes) as.integer(modes[[1]]$mode$nxhat)
  f_nmodes <- function(modes) length(modes)
  f_effnm  <- function(modes) {
    msize <- purrr::map_dbl(modes, function(m) m$size)
    pmode <- msize / sum(msize)
    2^sum(pmode * log2(1 / pmode))
  }

  tibble::tibble(
    nwords     = purrr::map_int(vocabularies, f_nwords),
    nmodes     = purrr::map_int(vocabularies, f_nmodes),
    eff_nmodes = purrr::map_dbl(vocabularies, f_effnm)
  )
}

# compute the absolute and effective number of modes for both historical and
# de novo n+1 vocabularies.
summarize_vocabulary_modes <- function(
  nplus1_vocabularies,
  nplus1_denovo_vocabularies
) {
  # summarize historical modes
  historical <- vocabulary_modes_summary(nplus1_vocabularies) |>
    tibble::add_column(
      language = wcs_languages,
      type     = "historical",
      .before  = 1
    )

  # summarize denovo modes
  denovo <- vocabulary_modes_summary(nplus1_denovo_vocabularies) |>
    tibble::add_column(
      language = wcs_languages,
      type     = "denovo",
      .before  = 1
    )

  rbind(historical, denovo)
}

# compute RD efficiency for a given set of vocabulary modes
vocabulary_modes_efficiency <- function(modes) {
  tibble::tibble(
    nwords    = purrr::map_int(modes, function(m) as.integer(m$mode$nxhat)),
    beta      = purrr::map_dbl(modes, function(m) m$mode$beta),
    frequency = purrr::map_int(modes, function(m) as.integer(m$size)),
    average_distortion = purrr::map_dbl(
      modes, function(m) average_distortion(m$mode)
    ),
    mutual_information = purrr::map_dbl(
      modes, function(m) mutual_information(m$mode)
    )
  ) |>
    dplyr::mutate(weight = frequency / sum(frequency))
}

# compute RD efficiencies for n+1 historical and de novo vocabulary modes
compute_vocabulary_modes_efficiency <- function(
  nplus1_vocabularies,
  nplus1_denovo_vocabularies
) {
  efficiency <- function(vocabularies, vtype, lnum) {
    vocabulary_modes_efficiency(vocabularies[[lnum]]) |>
      tibble::add_column(
        language = lnum,
        type     = vtype,
        .before  = 1
      )
  }

  # summarize historical modes
  historical <- purrr::map_dfr(seq_along(nplus1_vocabularies),
    function(lnum) efficiency(nplus1_vocabularies, "historical", lnum)
  )

  # summarize denovo modes
  denovo <- purrr::map_dfr(seq_along(nplus1_denovo_vocabularies),
    function(lnum) efficiency(nplus1_denovo_vocabularies, "denovo", lnum)
  )

  rbind(historical, denovo)
}
