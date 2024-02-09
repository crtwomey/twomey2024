#
# word_cluster_generation.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# returns a tibble with a row for every term in every language in the WCS
build_wcs_terms <- function(wcs) {
  # language vocabularies filtered by minimum ussage and minimum assignment
  wcs_lang <- lapply(wcs_languages, function(lnum) {
    wcs_language(wcs$languages, lnum)
  })

  # language-specific inferred communicative needs, p(x)
  lang_pxs <- lapply(wcs_languages, function(lnum) {
    lang_px(wcs, lnum)
  })

  # p(x|xhat) propto p(xhat|x) * p(x)
  px_xhats <- lapply(wcs_languages, function(lnum) {
    pxhat_x <- wcs_lang[[lnum]]$pw_c
    pxxhat  <- pxhat_x * lang_pxs[[lnum]]
    sweep(pxxhat, 2, colSums(pxxhat), "/")
  })

  # color terms from every WCS language
  tibble(
    language = do_c_sapply(wcs_languages, function(lnum) {
      rep(lnum, wcs_lang[[lnum]]$nwords)
    }),
    nterms = do_c_sapply(wcs_languages, function(lnum) {
      rep(wcs_lang[[lnum]]$nwords, wcs_lang[[lnum]]$nwords)
    }),
    term = do_c_sapply(wcs_languages, function(lnum) {
      seq_len(wcs_lang[[lnum]]$nwords)
    }),
    pxhat_x = do_c_sapply(wcs_languages, function(lnum) {
      lapply(seq_len(ncol(wcs_lang[[lnum]]$pw_c)), function(i) {
        wcs_lang[[lnum]]$pw_c[, i]
      })
    }),
    px_xhat = do_c_sapply(wcs_languages, function(lnum) {
      lapply(seq_len(ncol(px_xhats[[lnum]])), function(j) {
        px_xhats[[lnum]][, j]
      })
    }),
    avg_color = map_chr(px_xhat, function(p) px_xhat2color(wcs, p))
  ) |>
    # add foci data
    rowwise() |>
      mutate(foci = list(wcs_lang[[language]]$foci[term, ])) |>
    ungroup()
}

# distance between all solution pairs
compute_term_distances <- function(
  wcs_terms,
  wasserstein_distance
) {
  # compute upper triangle of distance matrix
  n   <- nrow(wcs_terms)
  ijs <- lapply(1:n, function(i) cbind(i, seq_len(i)))
  ijs <- do.call("rbind", ijs)

  term_distances <- parallel_lapply(seq_len(nrow(ijs)), function(k) {
    i <- ijs[k, 1]
    j <- ijs[k, 2]
    wasserstein_distance(
      wcs_terms$px_xhat[[i]],
      wcs_terms$px_xhat[[j]]
    )
  })

  term_distances
}

# pairwise adjacency matrix between all wcs terms
compute_term_adjacency <- function(
  wcs_terms,
  term_distances
) {
  deltas <- do.call("c", term_distances)
  n      <- nrow(wcs_terms)
  delta  <- matrix(Inf, n, n)
  delta[upper.tri(delta, diag = TRUE)] <- deltas
  delta[lower.tri(delta)] <- t(delta)[lower.tri(delta)]

  # full adjacency matrix
  adjacency <- exp(-term_adjacency_distance_scale * delta)
  diag(adjacency) <- 0

  adjacency
}

# graph clustering based on modularity (greedy optimization)
compute_word_clusters <- function(
  wcs,
  wcs_terms,
  term_adjacency
) {
  agraph <- graph.adjacency(term_adjacency,
    mode     = "undirected",
    weighted = TRUE
  )
  opt <- cluster_fast_greedy(agraph)
  grp <- igraph::groups(opt)

  # reorder and rename groups by number of languages in that category
  grp_num_langs <- map_int(grp, function(g) {
    length(unique(wcs_terms$language[g]))
  })
  grp_order  <- order(grp_num_langs, decreasing = TRUE)
  grp        <- grp[grp_order]
  grp_map    <- tibble(original = names(grp), renamed = seq_along(grp))
  grp_member <- as_vector(
    tibble(original = as.character(membership(opt))) |>
      left_join(grp_map, by = "original") |>
      select(renamed)
  )
  names(grp_member) <- NULL
  names(grp)        <- seq_along(grp)
  grp_num_langs     <- grp_num_langs[grp_order]

  grp_avg_pxhat_x <- lapply(grp, function(g) {
    colMeans(do.call("rbind", wcs_terms$pxhat_x[g]))
  })

  grp_avg_px_xhat <- lapply(grp, function(g) {
    colMeans(do.call("rbind", wcs_terms$px_xhat[g]))
  })

  grp_avg_colors <- map_chr(grp_avg_px_xhat, function(px_xhat) {
    px_xhat2color(wcs, px_xhat)
  })

  list(
    opt             = opt,
    grp             = grp,
    grp_member      = grp_member,
    grp_num_langs   = grp_num_langs,
    grp_avg_pxhat_x = grp_avg_pxhat_x,
    grp_avg_px_xhat = grp_avg_px_xhat,
    grp_avg_colors  = grp_avg_colors
  )
}

# wcs word assignments based on clustering
compute_wcs_words <- function(
  wcs,
  wcs_terms,
  word_clusters
) {
  # language-specific inferred communicative needs, p(x)
  lang_pxs <- lapply(wcs_languages, function(lnum) {
    lang_px(wcs, lnum)
  })

  p2c <- function(px_xhat) px_xhat2color(wcs, px_xhat)

  # find and merge synonyms
  wcs_terms |>
    mutate(word = word_clusters$grp_member) |>
    arrange(language, word) |>
    select(language, word, pxhat_x, foci) |>
    # average pxhat_x and foci by assigned word cluster
    group_by(language, word) |>
      mutate(
        pxhat_x = list(colMeans(
          matrix(do.call("rbind", pxhat_x), ncol = 330)
        )),
        foci = list(colMeans(
          matrix(do.call("rbind", foci), ncol = 3)
        ))
      ) |>
    ungroup() |>
    # remove synonyms (duplicates)
    distinct(language, word, .keep_all = TRUE) |>
    # add word px_xhat and average color
    rowwise() |>
      mutate(
        px_xhat = list(
          pxhat_x * lang_pxs[[language]] / sum(
            pxhat_x * lang_pxs[[language]]
          )
        ),
        avg_color = p2c(px_xhat)
      ) |>
    ungroup() |>
    # compute nwords
    group_by(language) |>
      mutate(nwords = length(language)) |>
    ungroup()
}

# frequency of each vocabulary size across the WCS languages
compute_nword_frequencies <- function(wcs_words) {
  wcs_words |>
    dplyr::group_by(language) |>
    dplyr::summarize(lang_nwords = head(nwords, 1)) |>
    dplyr::group_by(lang_nwords) |>
    dplyr::summarize(nlangs = length(lang_nwords))
}
