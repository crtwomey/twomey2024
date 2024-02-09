#
# wcs_plot.R
#
# All code is open source under the GPL v3 (see LICENSE)
#

plot_rd <- function(wcs, rd,
  threshold = 1 / rd$nxhat,
  ...
) {
  px_xhat <- NULL
  if (!exists("px_xhat", rd)) {
    px_xhat <- t(rd$pxhat_x) * rd$px
    px_xhat <- t(px_xhat) / colSums(px_xhat)
    rownames(px_xhat) <- NULL
    colnames(px_xhat) <- NULL
  } else {
    px_xhat <- rd$px_xhat
  }

  xhat_colors <- map_chr(seq_len(nrow(px_xhat)), function(i) {
    px_xhat2color(wcs, px_xhat[i, ])
  })

  wnums <- seq_len(nrow(rd$pxhat_x))
  polys <- lapply(wnums, function(wnum) {
    pxhat_x <- rd$pxhat_x[wnum, ]
    above_threshold <- map_lgl(pxhat_x, function(p) {
      p > threshold
    })
    boundary_polygons(wcs, above_threshold) |>
      add_column(color = xhat_colors[wnum])
  })

  for (i in seq_along(polys)[-1]) {
    polys[[i]]$group <- polys[[i]]$group + max(polys[[i - 1]]$group)
    polys[[i]]$path  <- polys[[i]]$path  + max(polys[[i - 1]]$path)
  }
  polys <- do.call("rbind", polys)

  plot_wcs_polys(wcs, polys, ...)
}

boundary_polygons <- function(wcs, chips_in_shape,
  offset0 = 0,
  inset   = 0.1
) {
  neighbors <- wcs_chip_neighbors(wcs)

  # neighborhood information within the shape
  bounds <- neighbors[chips_in_shape, ]

  coords <- wcs$chips$coordinates
  xs <- list(
    "l" = c(-1, -1),
    "t" = c(-1,  1),
    "r" = c(1,  1),
    "b" = c(1, -1)
  )
  ys <- list(
    "l" = c(1, -1),
    "t" = c(-1, -1),
    "r" = c(-1,  1),
    "b" = c(1,  1)
  )

  to_edge <- function(cnum, edge, type) {
    x0 <- ifelse(coords[cnum, 2] == 1, offset0, coords[cnum, 2])
    y0 <- coords[cnum, 1]
    tibble(cnum, vertex = c("start", "end"),
      x = xs[[edge]] * 0.5 + x0,
      y = ys[[edge]] * 0.5 + y0,
      edge, type
    )
  }

  lvls_ltrb <- c("l", "t", "r", "b")
  lvls_oc   <- c("open", "closed")
  nodes <- bounds |>
    mutate(across(!cnum,
      ~ ifelse(is.na(.), "open", ifelse(. %in% cnum, "none", "closed"))
    )) |>
    pivot_longer(!cnum, names_to = "edge", values_to = "type") |>
    filter(type != "none") |>
    mutate(
      edge = factor(edge, levels = lvls_ltrb),
      type = factor(type, levels = lvls_oc)
    ) |>
    pmap_dfr(to_edge) |>
    mutate(z = zorder_coord(x + 0.5, y + 0.5))

  #--------------------------------------------------------------------------
  # edge lookup table for matching ends of segments together
  starts <- nodes |> filter(vertex == "start") |> select(-vertex)
  ends   <- nodes |> filter(vertex == "end")   |> select(-vertex)

  nz     <- max(nodes$z)
  lookup <- vector(mode = "integer", length = 4 * nz)
  dim(lookup) <- c(nz, 4)
  lookup[cbind(starts$z, starts$edge)] <- seq_len(nrow(starts))

  #--------------------------------------------------------------------------
  # polygon paths

  # valid directions with precedence
  #   r <-/-> l
  #   t <-/-> b
  ltrb <- list( #     1 2 3 4
    c(1, 2, 4), # l : l t _ b
    c(2, 3, 1), # t : l t r _
    c(3, 4, 2), # r : _ t r b
    c(4, 1, 3)  # b : l _ r b
  )

  # find all closed polygons
  next_node <- rep(NA, nrow(starts))
  paths     <- list()
  i         <- 1
  repeat {
    k       <- 1
    path    <- rep(NA, nrow(starts))
    path[k] <- i
    repeat {
      js <- lookup[ends$z[i], ltrb[[ends$edge[i]]]]
      j  <- head(js[js != 0], 1)
      if (length(j) != 1) {
        warning("path following failed")
        return(NULL)
      }
      next_node[i] <- j
      path[k]      <- j
      if (is.na(next_node[j])) {
        i <- j
        k <- k + 1
      } else {
        break
      }
    }
    paths  <- c(paths, list(c(path[seq_len(k)], path[1])))
    unused <- which(is.na(next_node))
    if (length(unused) > 0) {
      i <- unused[1]
    } else {
      break
    }
  }

  # record groupings
  polys <- map_dfr(seq_along(paths), function(i) {
    starts[paths[[i]], ] |>
      select(-z) |>
      add_column(group = i, .before = 1)
  })

  #--------------------------------------------------------------------------
  # identify path segments and path tips (changes from one type of segment to
  # another within the same polygon grouping)
  k       <- 1
  tips    <- c()
  path    <- rep(NA, nrow(polys))
  path[1] <- k
  for (i in 2:length(path)) {
    same_group <- identical(polys$group[i - 1], polys$group[i])
    same_type  <- identical(polys$type[i - 1], polys$type[i])
    if (!same_group || !same_type) {
      k <- k + 1
    }
    if (same_group && !same_type) {
      tips <- c(tips, i)
    }
    path[i] <- k
  }
  polys <- polys |> add_column(path)

  #--------------------------------------------------------------------------
  # start and finish indexes for each group
  group_lims <- lapply(unique(polys$group), function(g) {
    gs <- which(polys$group == g)
    c(head(gs, 1), tail(gs, 1))
  })

  prev_edge <- map_int(seq_len(nrow(polys)), function(i) {
    if (i == group_lims[[polys$group[i]]][1]) {
      group_lims[[polys$group[i]]][2]
    } else {
      i - 1L
    }
  })

  lookup_inset_x <- matrix(c(
    1, 1, NA, 1,
    1, 0, -1, NA,
    NA, -1, -1, -1,
    1, NA, -1, 0
  ), nrow = 4, ncol = 4, byrow = TRUE)

  lookup_inset_y <- matrix(c(
    0, 1, NA, -1,
    1, 1, 1, NA,
    NA, 1, 0, -1,
    -1, NA, -1, -1
  ), nrow = 4, ncol = 4, byrow = TRUE)

  # open / closed inset correction
  lookup_inset_oc_x <- matrix(c(
    0, 1, 0, 1,
    1, 1, 1, 1
  ), nrow = 2, ncol = 4, byrow = TRUE)

  lookup_inset_oc_y <- matrix(c(
    1, 0, 1, 0,
    1, 1, 1, 1
  ), nrow = 2, ncol = 4, byrow = TRUE)

  inset_xs <- map_dbl(seq_len(nrow(polys)), function(i) {
    ep <- polys$edge[prev_edge[i]]
    ei <- polys$edge[i]
    tp <- polys$type[prev_edge[i]]
    ti <- polys$type[i]
    lookup_inset_x[ep, ei] *
      lookup_inset_oc_x[ti, ei] *
      lookup_inset_oc_x[tp, ep] *
      inset
  })

  inset_ys <- map_dbl(seq_len(nrow(polys)), function(i) {
    ep <- polys$edge[prev_edge[i]]
    ei <- polys$edge[i]
    tp <- polys$type[prev_edge[i]]
    ti <- polys$type[i]
    lookup_inset_y[ep, ei] *
      lookup_inset_oc_y[ti, ei] *
      lookup_inset_oc_y[tp, ep] *
      inset
  })

  # inset edges
  polys <- polys |> mutate(
    x = x + inset_xs,
    y = y + inset_ys
  )

  #--------------------------------------------------------------------------
  # fix repeated path tips
  path_tips       <- polys[tips, ]
  path_tips$path  <- polys$path[tips - 1]
  path_tips$type  <- polys$type[tips - 1]
  path_tips$group <- polys$group[tips - 1]

  # place tips at end of their respective paths
  polys <- rbind(polys, path_tips)
  polys <- polys[order(polys$path), ]

  polys
}

wcs_chip_neighbors <- function(wcs) {
  # WCS grid : 10 rows x 41 columns
  grd <- matrix(NA, nrow = 10, ncol = 41)
  grd[wcs$chips$coordinates] <- seq_len(330)

  # padding between achromatic central column (col 1) and main color grid
  # (columns 1 -- 41)
  grd <- cbind(grd[, 1], rep(NA, 10), grd[, -1])

  # padding on all boundaries
  grd <- rbind(
    rep(NA, 44),
    cbind(rep(NA, 10), grd, rep(NA, 10)),
    rep(NA, 44)
  )

  ijs <- expand.grid(2:(nrow(grd) - 1), 2:(ncol(grd) - 1))

  map_dfr(seq_len(nrow(ijs)), function(ij) {
    i <- ijs[ij, 1]
    j <- ijs[ij, 2]
    cnum <- grd[i, j]
    if (is.na(cnum)) return(NULL)
    tibble(cnum = cnum,
      l  = grd[i, j - 1],
      t  = grd[i - 1, j],
      r  = grd[i, j + 1],
      b  = grd[i + 1, j]
    )
  }) |> arrange(cnum)
}

zorder_coord <- function(x, y) {
  masks  <- c(0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF)
  shifts <- c(1, 2, 4, 8)

  for (j in rev(seq_along(shifts))) {
    x <- bitwAnd(bitwOr(x, bitwShiftL(x, shifts[j])), masks[j])
    y <- bitwAnd(bitwOr(y, bitwShiftL(y, shifts[j])), masks[j])
  }

  bitwOr(x, bitwShiftL(y, 1))
}

plot_wcs_polys <- function(wcs, datapoly,
  aspect_ratio   = 1 / 3.6,
  point_size     = 5,
  axis_text_size = 10,
  offset0        = 0,
  poly_alpha     = 0.75,
  boundry_lwd    = 0.25,
  label          = NULL
) {
  y_lims <- seq(1, 10)
  y_max  <- 10
  x_max  <- 41.5
  ggplot(datapoly, aes(x = x, y = y)) +
    geom_polygon(
      aes(group = group, fill = color),
      alpha = poly_alpha
    ) +
    geom_path(aes(
      group = path,
      color = color,
      alpha = ifelse(type == "closed", 1, 0)
    ), linetype = "solid", lwd = boundry_lwd) +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_y_reverse(NULL,
      breaks = y_lims,
      labels = LETTERS[seq_along(y_lims)],
      expand = c(0.1, 0.1),
      limits = c(10.5, 0.5)
    ) +
    scale_x_continuous(NULL,
      breaks = c(offset0, seq(1, 40) + 1),
      labels = c(0, seq(1, 40)),
      expand = c(0.01, 0.01),
      limits = c(-0.5, x_max)
    ) +
    annotate("rect",
      xmin = -0.5, xmax = 0.5, ymin = 0.5, ymax = y_max + 0.5,
      color = "gray", fill = NA, lwd = 0.1
    ) +
    annotate("rect",
      xmin = 1.5, xmax = 41.5, ymin = 1.5, ymax = y_max - 0.5,
      color = "gray", fill = NA, lwd = 0.1
    ) +
    annotate(
      geom  = "text", x = 41.5, y = 0.5, hjust = 1, size = 3,
      label = label
    ) +
    theme_minimal() +
    theme(
      aspect.ratio     = aspect_ratio,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        size   = axis_text_size, angle = 90, hjust = 1, vjust = 0.5,
        margin = margin(t = -5)
      ),
      axis.text.y = element_text(
        size   = axis_text_size, hjust = 0.5, vjust = 0.5,
        margin = margin(r = -2)
      )
    ) +
    theme(legend.position = "none") +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
}
