# Slide versions of figures for presentation/proposal-presentation.qmd.
#
# save_slide() re-themes a ggplot for the slides: larger text, no plot title or
# subtitle (the slide title does that job) and the legend at the bottom in two rows.
# Long facet labels are wrapped in the notebooks with labeller = label_wrap_gen(16).
# It uses the default "sans" family (Helvetica / Arial), the same font as the slides.
# The notebooks' own HTML output is unchanged; only the saved PNGs are restyled.

library(ggplot2)

theme_slides <- function(base_size = 18, legend = "bottom") {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      legend.position  = legend,
      legend.text      = element_text(hjust = 0),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      strip.text       = element_text(size = rel(0.9), face = "bold"),
      panel.grid.minor = element_blank()
    )
}

save_slide <- function(p, file, width = 11, height = 5.5, base_size = 18,
                       legend = "bottom", legend_rows = if (legend == "bottom") 2 else NULL,
                       dpi = 150) {
  p <- p + theme_slides(base_size, legend)
  if (!is.null(legend_rows)) {
    p <- p + guides(colour = guide_legend(nrow = legend_rows),
                    fill   = guide_legend(nrow = legend_rows))
  }
  ggsave(file, p, width = width, height = height, dpi = dpi, bg = "white")
}
