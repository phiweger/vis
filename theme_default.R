# https://github.com/stan-dev/bayesplot/blob/master/R/bayesplot-colors.R

scheme_level_names <- function() {
  c("light",
    "light_highlight",
    "mid",
    "mid_highlight",
    "dark",
    "dark_highlight")
}

red =
list("#DCBCBC", "#C79999", "#B97C7C", "#A25050", "#8F2727", "#7C0000")

library(ggplot2)

df <- read.table('conserved.tsv', sep='\t', stringsAsFactors=F, header=T)
names(df) = c('index', 'entropy', 'sum')

theme_default <- function(base_size = getOption("bayesplot.base_size", 12),
                          base_family = getOption("bayesplot.base_family", "serif")) {
  theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 0.4),
      axis.ticks = element_line(size = 0.3),
      strip.background = element_blank(),
      strip.text = element_text(size = rel(0.9)),
      strip.placement = "outside",
      # strip.background = element_rect(fill = "gray95", color = NA),
      panel.spacing = unit(1.5, "lines"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.text = element_text(size = 13),
      legend.text.align = 0,
      legend.key = element_blank()
    )
}

ggplot(df[df$entropy>0,], aes(x=index, y=entropy, color=entropy)) +
    geom_hline(yintercept=1, color='#8F2727') +
    geom_point(size=0.5) +
    scale_color_gradient(high='#8F2727', low='#DCBCBC') +
    theme_default()