#
# word_cluster_assignment.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# note that the order of the arguments needs to match the order of the column
# names for the data frames this function is applied to.
mode2terms <- function(p2c,
  language, type, parent_mode_id, mode_id, mode_nterms, mode_freq, mode
) {
  nterms <- mode$nxhat
  tibble::tibble(
    language       = rep(language, nterms),
    type           = rep(type, nterms),
    parent_mode_id = rep(parent_mode_id, nterms),
    mode_id        = rep(mode_id, nterms),
    mode_nterms    = rep(mode_nterms, nterms),
    mode_freq      = rep(mode_freq, nterms),
    mode_term      = seq_len(nterms),
    pxhat_x        = lapply(seq_len(nterms), function(i) mode$pxhat_x[i, ]),
    foci           = lapply(seq_len(nterms), function(i) mode$xhat[i, ])
  ) |> dplyr::rowwise() |>
    dplyr::mutate(
      px_xhat   = list(pxhat_x * mode$px / sum(pxhat_x * mode$px)),
      avg_color = p2c(px_xhat)
    )
}

# returns a helper function for constructing term mode tables
term_modes_table_builder <- function(wcs) {
  # language-specific inferred communicative needs, p(x)
  lang_pxs <- lapply(wcs_languages, function(lnum) {
    lang_px(wcs, lnum)
  })

  p2c <- function(px_xhat) px_xhat2color(wcs, px_xhat)
  m2t <- function(...) mode2terms(p2c, ...)

  function(language, type, parent_id, mode_id, nterms, freq, modes) {
    tibble::tibble(
      language       = language,
      type           = type,
      parent_mode_id = parent_id,
      mode_id        = mode_id,
      mode_nterms    = nterms,
      mode_freq      = freq,
      mode           = modes
    ) |> purrr::pmap_df(m2t)
  }
}

# historical n term modes table
build_historical_n_term_modes_table <- function(
  table_builder,
  rd_optimal_vocabularies
) {
  nterms <- purrr::map_int(rd_optimal_vocabularies, function(l) l$nxhat)
  table_builder(
    language  = seq_along(rd_optimal_vocabularies),
    type      = "historical_n",
    parent_id = NA,
    mode_id   = 1,
    nterms    = nterms,
    freq      = 1,
    modes     = rd_optimal_vocabularies
  )
}

# de novo n term modes table
build_denovo_n_term_modes_table <- function(
  table_builder,
  denovo_vocabularies
) {
  dfs <- parallel_lapply(seq_along(denovo_vocabularies), function(lnum) {
    modes  <- denovo_vocabularies[[lnum]]
    nterms <- purrr::map_int(modes, function(m) as.integer(m$mode$nxhat))
    freq   <- purrr::map_int(modes, function(m) as.integer(m$size))
    table_builder(
      language  = lnum,
      type      = "denovo_n",
      parent_id = NA,
      mode_id   = seq_along(modes),
      nterms    = nterms,
      freq      = freq,
      modes     = lapply(modes, function(m) m$mode)
    )
  })
  do.call("rbind", dfs)
}

# historical n+1 term modes table
build_historical_nplus1_term_modes_table <- function(
  table_builder,
  nplus1_vocabularies
) {
  dfs <- parallel_lapply(seq_along(nplus1_vocabularies), function(lnum) {
    modes  <- nplus1_vocabularies[[lnum]]
    nterms <- purrr::map_int(modes, function(m) as.integer(m$mode$nxhat))
    freq   <- purrr::map_int(modes, function(m) as.integer(m$size))
    table_builder(
      language  = lnum,
      type      = "historical_nplus1",
      parent_id = 1,
      mode_id   = seq_along(modes),
      nterms    = nterms,
      freq      = freq,
      modes     = lapply(modes, function(m) m$mode)
    )
  })
  do.call("rbind", dfs)
}

# de novo n+1 term modes table
build_denovo_nplus1_term_modes_table <- function(
  table_builder,
  denovo_nplus1_vocabularies
) {
  dfs <- parallel_lapply(seq_along(denovo_nplus1_vocabularies), function(lnum) {
    parent_modes <- denovo_nplus1_vocabularies[[lnum]]
    purrr::map_df(seq_along(parent_modes), function(parent_id) {
      child_modes <- denovo_nplus1_vocabularies[[lnum]][[parent_id]]
      nterms <- purrr::map_int(child_modes, function(m) {
        as.integer(m$mode$nxhat)
      })
      freq <- purrr::map_int(child_modes, function(m) as.integer(m$size))
      table_builder(
        language  = lnum,
        type      = "denovo_nplus1",
        parent_id = parent_id,
        mode_id   = seq_along(child_modes),
        nterms    = nterms,
        freq      = freq,
        modes     = lapply(child_modes, function(m) m$mode)
      )
    })
  })
  do.call("rbind", dfs)
}

# combine all term tables
build_term_modes_table <- function(
  historical_n_term_modes_table,
  historical_nplus1_term_modes_table,
  denovo_n_term_modes_table,
  denovo_nplus1_term_modes_table
) {
  rbind(
    historical_n_term_modes_table,
    historical_nplus1_term_modes_table,
    denovo_n_term_modes_table,
    denovo_nplus1_term_modes_table
  )
}

# compute assignments of terms to words
compute_modes_table_word_assignments <- function(
  assign_word,
  term_modes_table
) {
  ws <- parallel_lapply(seq_len(nrow(term_modes_table)), function(i) {
    assign_word(term_modes_table$px_xhat[[i]])
  })
  word_assignments <- do.call("c", ws)

  term_modes_table |>
    tibble::add_column(word = word_assignments) |>
    dplyr::group_by(type, language, parent_mode_id, mode_id) |>
    dplyr::mutate(nwords = length(unique(word))) |>
    dplyr::ungroup()
}
