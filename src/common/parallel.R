#
# parallel.R
#
# Helper functions for parallel computation.
# All code is open source under the GPL v3 (see LICENSE)
#

# run parallel
parallel_lapply <- function(v, f,
  max_ncores = parallelly::availableCores(omit = 1)
) {
  parallel::mclapply(v, f,
    mc.cores = min(length(v), max_ncores)
  )
}

parallel_rb_lapply <- function(...) {
  r <- parallel_lapply(...)
  do.call("rbind", r)
}

parallel_cb_lapply <- function(...) {
  r <- parallel_lapply(...)
  do.call("cbind", r)
}

parallel_c_lapply <- function(...) {
  r <- parallel_lapply(...)
  do.call("c", r)
}
