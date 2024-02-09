#
# word_transitions.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# builds a table of connections between parent and child terms
build_n_to_nplus1_table <- function(
  word_modes_table
) {
  # return historical or denovo depending on the passed type string
  get_type <- function(s) stringr::str_split_fixed(s, "_", 2)[, 1]

  n_df <- word_modes_table |>
    dplyr::filter(grepl("n$", type)) |>
    dplyr::select(-parent_mode_id) |>
    dplyr::mutate(
      type           = get_type(type),
      parent_mode_id = mode_id
    ) |>
    dplyr::select(-mode_id)

  nplus1_df <- word_modes_table |>
    dplyr::filter(grepl("nplus1$", type)) |>
    dplyr::mutate(
      type        = get_type(type),
      nplus1words = nwords
    ) |>
    dplyr::select(-nwords)

  nplus1_df |>
    dplyr::left_join(
      n_df |> dplyr::select(type, language, parent_mode_id, mode_term, word),
      by     = c("type", "language", "parent_mode_id", "mode_term"),
      suffix = c("_n1", "_n0")
    ) |>
    dplyr::group_by(type, language, parent_mode_id, mode_id) |>
    dplyr::mutate(
      nwords = sum(!is.na(word_n0))
    ) |>
    dplyr::ungroup()
}

# reformats the n -> n+1 table as a table of transitions
build_n_to_nplus1_transitions <- function(
  n_to_nplus1_table
) {
  word_levels <- max(
    n_to_nplus1_table[c("word_n0", "word_n1")], na.rm = TRUE
  )

  n_to_nplus1_table |>
    dplyr::filter(nplus1words == nwords + 1) |>
    dplyr::select(type, nplus1words, word_n0, word_n1, mode_freq) |>
    dplyr::mutate(
      word_n0 = factor(word_n0, levels = seq_len(word_levels)),
      word_n1 = factor(word_n1, levels = seq_len(word_levels))
    ) |>
    dplyr::group_by(type, nplus1words) |>
      tidyr::pivot_wider(
        names_from  = word_n0,
        values_from = mode_freq,
        values_fn   = sum,
        names_sort  = TRUE
      ) %>%
    dplyr::ungroup() |>
    dplyr::rename(added = "NA") |>
    dplyr::group_by(type) |>
      dplyr::arrange(nplus1words, word_n1) |>
      tidyr::complete(nplus1words, word_n1) |>
    dplyr::ungroup()
}

# reformats the transition table as a set of transition matricies. Returns a
# list of lists, at the top level one for each type (historical and de novo),
# and for each type a list of transition matricies from the available range of
# n -> n+1 transitions (e.g. 3 -> 4, 4 -> 5, etc.)
build_joint_transition_matricies <- function(
  n_to_nplus1_transitions
) {
  nplus1_min <- min(n_to_nplus1_transitions$nplus1words)
  nplus1_max <- max(n_to_nplus1_transitions$nplus1words)

  lapply(unique(n_to_nplus1_transitions$type), function(tp) {
    lapply(nplus1_min:nplus1_max, function(np1) {
      np1_transitions <- as.matrix(n_to_nplus1_transitions |>
        dplyr::filter(type == tp, nplus1words == np1) |>
        dplyr::select(-type, -nplus1words, -word_n1) |>
        dplyr::mutate(
          dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0))
        )
      )
      np1_transitions / sum(np1_transitions)
    })
  })
}

# keeps the same structure as the joint transition tables (see
# build_joint_transition_matricies above), but normalizes each column so we
# have the conditional probability that a word in a vocabulary of size n
# transitions to a word in a vocabulary of size n+1.
build_conditional_transition_matricies <- function(
  n_to_nplus1_transitions
) {
  nplus1_min <- min(n_to_nplus1_transitions$nplus1words)
  nplus1_max <- max(n_to_nplus1_transitions$nplus1words)

  lapply(unique(n_to_nplus1_transitions$type), function(tp) {
    lapply(nplus1_min:nplus1_max, function(np1) {

      np1_transitions <- as.matrix(n_to_nplus1_transitions |>
        dplyr::filter(type == tp, nplus1words == np1) |>
        dplyr::select(-type, -nplus1words, -word_n1) |>
        dplyr::mutate(
          dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0))
        )
      )
      sweep(np1_transitions, 2, colSums(np1_transitions), "/")
    })
  })
}
