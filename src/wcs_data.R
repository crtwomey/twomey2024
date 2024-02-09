#
# wcs_data.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

# load and aggregate WCS data
load_wcs_data <- function(
  wcs_chips_file,
  wcs_vocabularies_file,
  wcs_info_file,
  wcs_infs_file
) {
  # load WCS p(word|color) maps and average Lab foci coordinates
  wcs        <- readRDS(wcs_vocabularies_file)
  nlanguages <- 110

  # load WCS chips to get their Lab coordinates
  wcs_chips <- readRDS(wcs_chips_file)
  nchips    <- 330

  # load WCS language information
  wcs_info <- read.table(wcs_info_file, sep = "\t", header = TRUE)

  # load best-fit WCS inferred distributions
  infs <- readRDS(wcs_infs_file)

  list(
    languages  = wcs,
    chips      = wcs_chips,
    info       = wcs_info,
    inferred   = infs,
    nlanguages = nlanguages,
    nchips     = nchips
  )
}
