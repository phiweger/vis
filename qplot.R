source('~/code/vis/theme_minimal.R')

library(ggplot2)

# Some functions for EDA of 1-dimensional data, e.g. (marginal) 
# distributions.

# example data.
# samples = rnorm(1e4, mean=0, sd=1)

# q .. quick
qhist <- function(data=NULL, label='foo', nbins=100, ratioscale=.5){
	df <- data.frame(label=data)
	binwidth <- dist(range(data)) / nbins # divide in 100 bins
	p <- ggplot(df, aes(x=label)) +
			geom_histogram(binwidth=binwidth) +
			theme_minimal() +
			theme(
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank())

	plotinfo <- ggplot_build(p)
	x <- as.numeric(dist(as.numeric(unlist(
		plotinfo$panel$ranges[[1]][1]
		))))
	y <- as.numeric(dist(as.numeric(unlist(
		plotinfo$panel$ranges[[1]][7]
		))))
	print(c(x, y))
	r = x / y * ratioscale
	p <- p + coord_fixed(ratio=r) # ratio = y / x
	return(p)
}

qfreq <- function(data=NULL, label='foo', nbins=100){
	df <- data.frame(label=data)
	binwidth <- dist(range(data)) / nbins # divide in 100 bins
	p <- ggplot(df, aes(x=label)) +
			geom_freqpoly(binwidth=binwidth) +
			theme_minimal()
	return(p)
}

qdens <- function(data=NULL, label='foo'){
	df <- data.frame(label=data)
	
	p <- ggplot(df, aes(x=label)) +
			geom_density() +
			theme_minimal()
	return(p)
}

# TODO
qscatter <- function(data=NULL, label1, label2){
	p <- ggplot(data, aes_string(x=label1, y=label2)) +
			geom_point() +
			theme_minimal()
	return(p)
}

# TODO:

# For timeseries.
# * calendar plot
# * sparklines

# Small multiples.
# Circular plot.
