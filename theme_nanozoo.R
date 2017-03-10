# http://docs.ggplot2.org/dev/vignettes/themes.html
# data
# https://cran.r-project.org/web/packages/HistData/HistData.pdf
# thread
# http://stackoverflow.com/questions/6736378/how-do-i-change-the-background-color-of-a-plot-made-with-ggplot2
# custom palette
# https://www.r-bloggers.com/creating-color-palettes-in-r/
# color cheat sheet
# https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf

theme_nanozoo = 
    theme_bw() +
    theme(
        axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "navy"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(colour = "grey40"),
        panel.grid.minor = element_blank()
        )

# colors:

# like 
# [bayesplot](https://cran.r-project.org/web/packages/bayesplot/vignettes/MCMC.html)
# package have 4 shades of each color, giving 12 in total
# also: how to get tones of one color in continuous cases

# grün: #9bd754
# blau: #2689ef
# pink: #ff1d6c
# grey62

# shades pink (dark to light):
# #f4abc5
# #ff1d6c # original
# #ff5993
# #f486ad
# #f4abc5


# blue, Economist-like background
# #e0efff

# google color and use google color picker to find shades
# have RColorBrewer palettes, i.e. continuous and discrete ones

# continuous white to all 3 colours
# continuous all 3 pairs of colours
# discrete all 3
# discrete steps (same as continuous from white to color)

data(iris)
df = iris

cols <- c(
    "virginica" = "#ff1d6c", 
    "versicolor" = "grey62", 
    "setosa" = "#f486ad"
    )

cols_shade_pink <- c(
    "virginica" = "#9bd754", 
    "versicolor" = "#2689ef", 
    "setosa" = "#ff1d6c"
    )

ggplot(df, aes(x=Sepal.Length, y=Petal.Width, color=Species)) +
    geom_jitter() +
    scale_colour_manual(values=cols) +
    theme_bw() +
    theme(panel.background = element_rect(fill = '#e0efff'))





th <- theme(
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    #panel.border=element_rect(color="grey50", fill=NA),
    axis.ticks=element_blank(),
    axis.text.y=element_blank(),
    axis.title=element_blank(),
    axis.ticks=element_blank()
    )

scale_colour_gradient(low = "#9bd754", high = "#ff1d6c", guide=F)



data(iris)
df = iris



df = as.data.frame(cbind(1:10, rbinom(10, 100, 0.2)))
names(df) = c('x', 'y')

ggplot(df, aes(x=x, y=y)) + geom_line()