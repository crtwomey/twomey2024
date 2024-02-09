#
# plotting.R
#
# Common plotting functions and defintions used throughout this codebase.
# All code is open source under the GPL v3 (see LICENSE)
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(paletteer)
  library(tikzDevice)
})

# set ggplot2 theme
theme_set(theme_classic() + theme(
  panel.grid.minor.y = element_line(),
  panel.grid.major.y = element_line(),
  panel.grid.minor.x = element_line(),
  panel.grid.major.x = element_line(),
) + theme(
  plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
))

scale_fill_continuous <- function(...) {
  scale_fill_viridis_c(..., option = "rocket")
}

scale_colour_continuous <- function(...) {
  scale_colour_viridis_c(..., option = "viridis")
}

# discrete color scheme options from
# https://github.com/EmilHvitfeldt/r-color-palettes
palette_choice <- c(
  "awtools::a_palette",
  "ggsci::dark_uchicago",
  "ggsci::uniform_startrek",
  "ggthemes::wsj_colors6",
  "ochRe::lorikeet",
  "rcartocolor::Bold",
  "RColorBrewer::Set1",
  "wesanderson::Darjeeling1",
  "yarrr::basel"
)[4]

scale_fill_discrete <- function(...) {
  scale_fill_paletteer_d(..., palette_choice)
}

scale_colour_discrete <- function(...) {
  scale_color_paletteer_d(..., palette_choice)
}

# output figure
fig <- function(figure, filename, width, height, use_tikz = FALSE) {
  ext <- ifelse(use_tikz, ".tex", ".pdf")
  f   <- ifelse(use_tikz, tikz, pdf)
  output_file <- file.path(fig_dir, paste0(filename, ext))
  f(output_file, width = width, height = height)
  print(figure)
  dev.off()
  output_file
}

# wrapper for executing base plot code
print.base_plot <- function(pl) pl()

# assign colors to WCS chip positions (defaults to WCS color stimuli)
wcs_color_points <- function(wcs, chip_colors = wcs$chips$RGB) {
  # WCS color chip coordinates (row and column)
  coords        <- as.data.frame(wcs$chips$coordinates)
  names(coords) <- c("y", "x")
  tibble(
    x     = coords$x,
    y     = coords$y,
    color = chip_colors
  )
}

# assign colors to WCS chip rects (defaults to WCS color stimuli)
wcs_color_chips <- function(wcs,
  chip_colors = wcs$chips$RGB,
  offset0     = 0
) {
  ids    <- factor(seq_along(wcs$chips$RGB))
  values <- tibble(
    id    = ids,
    color = chip_colors
  )
  coords    <- wcs$chips$coordinates
  make_chip <- function(coord) {
    x0 <- ifelse(coord[2] == 1, offset0, coord[2])
    y0 <- coord[1]
    tibble(
      x  = c(-1, -1, 1, 1) * 0.5 + x0,
      y  = c(-1, 1, 1, -1) * 0.5 + y0
    )
  }
  chips <- map_dfr(ids, function(i) {
    make_chip(coords[i, ]) |> add_column(id = i)
  })
  merge(values, chips, by = c("id"))
}

chip_colors_from_mapping <- function(wcs, px_xhat, pxhat_x) {
  color <- px_xhat2color(wcs, px_xhat)
  map_chr(pxhat_x, function(p) {
    adjustcolor(color, alpha.f = p)
  })
}

# plot_wcs_chip_colors
plot_wcs_color_chips <- function(chips,
  chip_colors    = chips$RGB,
  point_size     = 5,
  axis_text_size = 10,
  offset0        = 0,
  outline        = "gray",
  cell_outline   = "white",
  bg             = NA,
  label          = NULL
) {
  y_lims <- seq(1, 10)
  y_max  <- max(chips$y)
  x_max  <- max(chips$x)
  ggplot(chips, aes(x = x, y = y)) +
    annotate("rect",
      xmin = -0.5, xmax = 0.5, ymin = 0.5, ymax = y_max,
      color = outline, fill = bg, lwd = 0.1
    ) +
    annotate("rect",
      xmin = 1.5, xmax = 41.5, ymin = 1.5, ymax = y_max - 1.0,
      color = outline, fill = bg, lwd = 0.1
    ) +
    geom_polygon(aes(fill = color, group = id),
      color = cell_outline, lwd = 0.05
    ) +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_y_reverse(NULL,
      breaks = y_lims,
      labels = LETTERS[seq_along(y_lims)],
      expand = c(0.1, 0.1)
    ) +
    scale_x_continuous(NULL,
      breaks = c(offset0, seq(1, x_max - 1) + 1),
      labels = c(0, seq(1, x_max - 1)),
      expand = c(0.01, 0.01)
    ) +
    annotate(
      geom  = "text", x = 41.5, y = 0.5, hjust = 1, size = 4,
      label = label
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        size   = axis_text_size, angle = 90, hjust = 1, vjust = 0.5,
        margin = margin(t = -5)
      ),
      axis.text.y = element_text(
        size   = axis_text_size, hjust = 0.5, vjust = 0.5,
        margin = margin(r = -2)
      ),
      plot.margin = margin(t = 0, b = 0, l = 0, r = 0)
    ) +
    theme(legend.position = "none") +
    coord_fixed()
}

# generate a reference plot of the WCS color stimuli
plot_wcs_color_stimuli <- function(wcs, ...) {
  wcs_color_chips(wcs) |>
    plot_wcs_color_chips(
      label = "WCS color stimuli",
      ...
    )
}

# plot a term
plot_term <- function(wcs, px_xhat, pxhat_x, ...) {
  wcs_color_chips(wcs,
    chip_colors_from_mapping(wcs, px_xhat, pxhat_x)
  ) |>
    plot_wcs_color_chips(...)
}

# plot a word (same as plotting a term)
plot_word <- plot_term

# standardized order for displaying color words
plotting_color_order <- c(
  8, 5, 11, 2, 1, 14, 13, 3, 15, 6, 12, 4, 7, 10, 9
)
