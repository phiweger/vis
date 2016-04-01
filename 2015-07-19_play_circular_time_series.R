library(magrittr)
library(ggplot2)
library(dplyr)

time = c(1:360)
value = abs(rnorm(360))

df = cbind(time, value) %>% as.data.frame()

ggplot(df, aes(x=time, y=value, col=time)) + geom_path() + coord_polar()


library(Quandl)
sp <- Quandl('YAHOO/INDEX_GSPC') #, collapse='monthly')

library(lubridate)
# > lubridate::year(sp$Date[1])
# [1] 2015
sp$Date = as.Date(sp$Date)
sp$year = sapply(sp$Date, year)


# http://stackoverflow.com/questions/7958298/how-do-you-convert-posix-date-to-day-of-year-in-r
# doy <- strftime(today, format = "%j")
sp$doj = as.numeric(strftime(sp$Date, format = "%j"))

library(data.table)
setnames(sp, 'Adjusted Close', 'adjclose')


sp$year = as.factor(sp$year)

# https://danieljhocking.wordpress.com/2014/12/03/lags-and-moving-means-in-dplyr/
sp %<>% dplyr::mutate(
	mean = rollapply(
		data=adjclose,
		FUN=mean,
		width=80)
	)

data = sp[(year(sp$Date) > 2005),]

ggplot(data, 
	aes(x=doj, y=adjclose)) + 
	geom_line(aes(col=year)) +
	scale_colour_hue(h=c(0, 360)) +
	geom_point(x=data$doj[1], y=data$adjclose[1]) +
	geom_point(x=data$doj[nrow(data)], y=data$adjclose[nrow(data)]) +
	coord_polar()
# http://docs.ggplot2.org/current/scale_hue.html

