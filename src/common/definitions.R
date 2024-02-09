#
# common.R
#
# Parameters and common utilities used throughout the project.
# All code is open source under the GPL v3 (see LICENSE)
#

# indexes of languages from each dataset
wcs_languages <- 1:110
bk_languages  <- 111:130

# total number of languages in each dataset
n_wcs <- length(wcs_languages)
n_bk  <- length(bk_languages)

# check which dataset a language comes from
is_wcs <- function(language) language <= 110
is_bk  <- function(language) language >  110

# minimum number of times a term was used to include it in a language's
# vocabulary. BK data only has one participant per language, so we can't
# require more than one use. In the WCS data we require at least two
# participants to have used the term to include it.
min_usage_threshold_bk  <- 1
min_usage_threshold_wcs <- 2

# return dataset dependent usage threshold
min_usage_threshold <- function(language) {
  ifelse(is_wcs(language),
    min_usage_threshold_wcs,
    min_usage_threshold_bk
  )
}

# minimum number of times a term is the best choice for a color
min_assignment_threshold <- 1

# load WCS language vocabulary, returning only words that are
# the best choice for at least min.assignments colors chips.
wcs_language <- function(wcs, language,
  min_assignments = min_assignment_threshold
) {
  foci   <- wcs[[language]]$foci
  pw_c   <- wcs[[language]]$pw_c
  nwords <- wcs[[language]]$nwords
  best   <- wcs[[language]]$best

  # determine which words best represent at least one color
  keep <- which(best >= min_assignments)

  return(list(
    foci   = foci[keep, ],
    pw_c   = pw_c[, keep],
    nwords = length(keep)
  ))
}

# determine which words in the term map pass the assignment threshold
word_assignment_threshold <- function(pw_c,
  min_assignments = min_assignment_threshold
) {
  assignment <- apply(pw_c, 1, function(p) head(which(p == max(p)), 1))
  total      <- sapply(seq_len(ncol(pw_c)), function(i) sum(assignment == i))
  return(total > min_assignments)
}

# all-pairs Euclidean distances
squared_euclidean_distance <- function(x, y) {
  s <- matrix(0, nrow(x), nrow(y))
  for (i in seq_len(ncol(x))) {
    s <- s + outer(x[, i], y[, i], "-")^2
  }
  return(s)
}

# returns the closest (in CIE Lab space) WCS chip to each row of x,
# where x is an n by 3 matrix of coordinates in CIE Lab.
get_nearest_chips <- function(wcs_chips, x) {
  d       <- squared_euclidean_distance(x, wcs_chips$Lab)
  nearest <- apply(d, 1, function(x) head(which(x == min(x)), 1))
  return(as.vector(nearest))
}

# tolerance threshold chosen based on diminishing returns in improvement
# to RMSE balanced with keeping entropy high (see Twomey et al. 2021 SI
# Fig. C2a).
rmse_threshold_tolerance <- 0.25

# language-specific inferences of px
lang_px <- function(wcs, language) {
  tolerances <- sapply(wcs$inferred[[language]]$invs, function(i) i$tol)
  rmse       <- sapply(wcs$inferred[[language]]$fits, function(f) f$RMSE)
  tol        <- tolerances <= rmse_threshold_tolerance
  opt        <- which(rmse == min(rmse[tol]) & tol)
  wcs$inferred[[language]]$invs[[opt]]$qx
}

# language-specific beta (need for precision) from previous inferences
lang_beta <- function(wcs, language) {
  tolerances <- sapply(wcs$inferred[[language]]$invs, function(i) i$tol)
  rmse       <- sapply(wcs$inferred[[language]]$fits, function(f) f$RMSE)
  tol        <- tolerances <= rmse_threshold_tolerance
  opt        <- which(rmse == min(rmse[tol]) & tol)
  wcs$inferred[[language]]$fits[[opt]]$opt$par
}

# convert p(x|xhat) to a displayable RGB color using the colorscience package
px_xhat2color <- function(wcs, px_xhat) {
  lab_colors <- colSums(wcs$chips$Lab * px_xhat)
  xyz_colors <- colorscience::Lab2XYZ(lab_colors, illuminant = "C")
  rgb_colors <- pmax(colorscience::XYZ2RGB(xyz_colors, illuminant = "C"), 0)
  rgb_colors <- pmin(rgb_colors, 1)
  rgb(rgb_colors[1], rgb_colors[2], rgb_colors[3], maxColorValue = 1)
}

# sapply a function f to each element of v, and c-bind the results
do_c_sapply <- function(v, f) do.call("c", sapply(v, f))

# parameter controlling exponential decay of distances for conversion from
# distances to adjacency matrix (compute_term_adjacency_matrix)
term_adjacency_distance_scale <- 1

# minimum mean squared distance between modes to consider them distinct
same_mode_threshold <- 100

# initial range for scale parameter converting distances to likelihoods
likelihood_distance_scales <- 10^seq(-5, 5, 0.01)

# cutoff for estimating range of variation in likelihoods
likelihood_change_threshold <- 1E-05

# initial frequency of new terms (for n + 1 analysis)
initial_frequency_of_new_term <- 1E-9

# number of de novo modes to generate per language
num_gen_denovo_modes <- 1000

# maximum number of emdw (earth mover distance solver) iterations
max_emdw_iterations <- 1000

# stopping criteria for Blahut-Arimoto (max iterations or min difference)
max_ba_iterations     <- 1000
min_ba_iteration_diff <- 1E-16

# example languages to use in Fig. 5
fig5_language_examples <- c(23, 9)
