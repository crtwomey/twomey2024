#
# wcs_distance.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# Generate a function that computes the earth mover's distance between two
# distributions using the WCS color chips' CIE Lab space coordinages.
build_wcs_distance_function <- function(wcs) {
  x <- wcs$chips$Lab
  function(px, qx) {
    emdw(x, px, x, qx,
      dist     = "euclidean",
      max.iter = max_emdw_iterations
    )
  }
}

# returns a function that takes a matrix of terms (one term per row) and
# returns an integer assignment to the minimum cost word.
build_word_assignment_function <- function(
  wcs_distance,
  word_clusters
) {
  word_px_xhat <- do.call("rbind", word_clusters$grp_avg_px_xhat)
  nwords       <- nrow(word_px_xhat)
  nstimuli     <- ncol(word_px_xhat)

  function(term_px_xhat) {
    term_px_xhat <- matrix(term_px_xhat, ncol = nstimuli)
    nterms       <- nrow(term_px_xhat)
    purrr::map_int(seq_len(nterms), function(i) {
      # compute the cost of assigning the term to each candidate word
      cs <- purrr::map_dbl(seq_len(nwords), function(j) {
        wcs_distance(term_px_xhat[i, ], word_px_xhat[j, ])
      })

      # assign the term to the minimum cost word
      head(which(cs == min(cs)), 1)
    })
  }
}
