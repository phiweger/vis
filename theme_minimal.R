# http://stackoverflow.com/questions/4094094/modifying-fonts-in-ggplot2

theme_min = function (size=10, font=NA, face='plain', 
    panelColor=backgroundColor, axisColor='#999999', 
    gridColor=gridLinesColor, textColor='black') 
{
    theme_text = function(...)
        ggplot2::theme_text(family=font, face=face, colour=textColor, 
            size=size, ...)

opts(
    axis.text.x = theme_text(),
    axis.text.y = theme_text(),
    axis.line = theme_blank(),
    axis.ticks = theme_segment(colour=axisColor, size=0.25),
    panel.border = theme_rect(colour=backgroundColor),
    legend.background = theme_blank(),
    legend.key = theme_blank(),
    legend.key.size = unit(1.5, 'lines'),
    legend.text = theme_text(hjust=0),
    legend.title = theme_text(hjust=0),
    panel.background = theme_rect(fill=panelColor, colour=NA),
    panel.grid.major = theme_line(colour=gridColor, size=0.33),
    panel.grid.minor = theme_blank(),
    strip.background = theme_rect(fill=NA, colour=NA),
    strip.text.x = theme_text(hjust=0),
    strip.text.y = theme_text(angle=-90),
    plot.title = theme_text(hjust=0),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'lines'))
}

##Create a custom font type. Could be 'F', 'TEST', whatever
# windowsFonts(F = windowsFont('Wide Latin'))

##and insert this line of code into the original code I list above: 
# + theme_min(font='F', size=10) 