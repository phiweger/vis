library(ggplot2)
library(magrittr)
library(dplyr)
library(ggthemes)

# fp = '/Users/owl/Desktop/data_molpath/tmp/y.csv'
fp = '/Users/owl/Desktop/data_molpath/tmp/...'
df = read.table(fp, sep=',', header=F, stringsAsFactors=F)

names(df) = c(
     'start',
     'type',
     'dist_rel',
     'coverage',
     'len',
     'len_rel',
     'num_switch',
     'num_tfbs',
     'num_dmr',
     'has_dmr',
     'num_cpg',
     'methyl_base_mean',
     'methyl_base_var',
     'methyl_diff_mean',
     'methyl_diff_var',
     'gene_name',
     'gene_id',
     'gene_len',
     'gene_strand',
     'gene_chrom',
     'gene_start',
     'gene_end',
     'expression',
     'padj',
     'compartment'
     )


df = df[df$padj < 0.05,]
# introduces NA if padj == nan
df = df[!is.na(df$gene_name),]
# removes NA


smry = df %>% 
     select(
          -gene_id, 
          -gene_name, 
          -gene_len,
          -gene_strand,
          -gene_chrom,
          -gene_start,
          -gene_end,
          -compartment, 
          -has_dmr,
          -len_rel,
          -start,
          -padj) %>% 
     group_by(type) %>% summarise_each(funs(mean))
     # group_by(type) %>% summarise_each(funs(mean, sd))
     # make smry_sd

# smry$bin = cut(smry$expression, breaks=seq(-2, 2, 0.5))
# for the figure display

# http://stackoverflow.com/questions/15215457/standardize-data-columns-in-r
t = scale(smry[, 2:ncol(smry)])

mock = matrix(ncol=3, nrow=1) %>% as.data.frame()
names(mock) = c('index', 'value', 'type')
index = c(1:ncol(t))

for (i in seq_along(1:nrow(t))){
     # print(i)
     name = smry$type[i]
     s = t[i,]

     u = cbind(index, s) %>% as.data.frame()
     names(u) = c('index', 'value')
     u %<>% mutate(type = name)
     mock = rbind(mock, u)
}

# remove 1st mock row
mock = mock[2:nrow(mock),]



# rename switches
cs_rename = function(cs_string){
     cs_string %<>% strsplit('_') %>% 
          unlist() %>% 
          '['(c(3,5)) %>% 
          paste(., collapse=' > ') 
     return(cs_string)
}



cs_from = function(cs_string){
     cs_string %<>% strsplit(' > ') %>% 
          unlist() %>% 
          '['(1)
     return(paste0('from ', cs_string))
}



cs_to = function(cs_string){
     cs_string %<>% strsplit(' > ') %>% 
          unlist() %>% 
          '['(2)
     return(paste0('to ', cs_string))
}



mock$type = sapply(mock$type, cs_rename)
smry$type = sapply(smry$type, cs_rename)

j = inner_join(mock, smry)


# t2 = t[3,]
# u = cbind(index, t2) %>% as.data.frame()
# names(u) = c('index', 'value')
# u %<>% mutate(type = 'a')

# t3 = t[4,]
# v = cbind(index, t3) %>% as.data.frame()
# names(v) = c('index', 'value')
# v %<>% mutate(type = 'b')

# w = rbind(u, v)

p = 
ggplot(j, aes(x=index, y=value)) + 
     geom_path() + 
     coord_polar() +
     scale_x_continuous(breaks=1:14) +
     scale_y_continuous(limits = c(-5, 5)) +
     theme(text = element_text(size=10)) + 
     facet_wrap(~ type)

fp_img_001 = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/radar_plot_switch_001.pdf'

pdf(fp_img_001, width=7, height=7)
p
dev.off()

# superimpose all: not too informative
# q = 
# ggplot(mock, aes(x=index, y=value)) + 
#      geom_path(aes(color=expression)) + 
#      coord_polar() +
#      scale_x_continuous(breaks=1:14) +
#      scale_y_continuous(limits = c(-5, 5)) +
#      theme(text = element_text(size=10)) +
#      facet_wrap(~ type)



fp_img_002 = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015-07-22_radar_plot_switch_002.pdf'

# > brewer.pal(name='RdYlGn', 11)
#  [1] "#A50026" "#D73027" "#F46D43" "#FDAE61" "#FEE08B" "#FFFFBF" "#D9EF8B"
#  [8] "#A6D96A" "#66BD63" "#1A9850" "#006837"

r = 
ggplot(j, aes(x=index, y=value)) + 
     geom_path(aes(color=expression), size=1) + 
     coord_polar() +
     scale_x_continuous(breaks=1:14) +
     scale_y_continuous(limits = c(-5, 5)) +
     scale_color_continuous(limits = c(-2, 2), low='#1A9850', high='#D73027') + #
     theme(
          text = element_text(size=10),  
          strip.background = element_blank()) +
     xlab('feature') + 
     ylab('value') +
     facet_wrap(~ type)

pdf(fp_img_002, width=7, height=7)
r
dev.off()





fp_img_003 = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015-07-24_radar_plot_switch_tss_002.pdf'

j$bin = cut(j$expression, seq(-2.5, 2.5, 0.5), right=F)
# > levels(j$bin)
# [1] "(-2,-1.5]" "(-1.5,-1]" "(-1,-0.5]" "(-0.5,0]"  "(0,0.5]"   "(0.5,1]"
# [7] "(1,1.5]"   "(1.5,2]"
# > length(levels(j$bin))
# [1] 8

cols = colorRampPalette(c('red', 'grey50', 'green'))(10) 
# cols = colorRampPalette(c('red', 'black', 'green'))(10) 
# include more colours if needed

# automate this step in a function!
# cols2 = c(
#      '[-2.5,-2)' = '#00FF00',
#      '[-2,-1.5)' = '#00C600',
#      '[-1.5,-1)' = '#008D00',
#      '[-1,-0.5)' = '#005400',
#      '[-0.5,0)' = '#001C00',
#      '[0,0.5)' = '#1C0000',
#      '[0.5,1)' = '#550000',
#      '[1,1.5)' = '#8D0000',
#      '[1.5,2)' = '#C60000',
#      '[2,2.5)' = '#FF0000'
#      )

cols2 = c(
     '[-2.5,-2)' = '#00FF00',
     '[-2,-1.5)' = '#1CE21C',
     '[-1.5,-1)' = '#38C638',
     '[-1,-0.5)' = '#54A954',
     '[-0.5,0)' = '#708D70',
     '[0,0.5)' = '#8D7070',
     '[0.5,1)' = '#A95454',
     '[1,1.5)' = '#C63838',
     '[1.5,2)' = '#E21C1C',
     '[2,2.5)' = '#FF0000'
     )

r = 
ggplot(
     j,
     #transform(j, discrete=cut(expression, seq(-2, 2, 0.5))), 
     aes(x=index, y=value)
     ) + 
     geom_path(aes(colour=factor(bin)), size=1) + 
     coord_polar() +
     theme_bw() +
     scale_x_continuous(breaks=1:14) +
     scale_y_continuous(limits=c(-5, 5)) +
     scale_color_manual(
          name='expression (log2fc)',
          limits=names(cols2),
          breaks=rev(names(cols2)),
          #values=brewer.pal(name='Set1', n=8)
          values=cols2
          ) +
     theme(
          text = element_text(size=10),  
          strip.background = element_blank(),
          strip.text = element_text(colour='black'),
          axis.text = element_text(colour='grey50'),
          axis.ticks = element_line(colour='grey50'),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.key = element_rect(colour=NA),
          axis.title = element_text(size=12)
          ) +
     xlab('feature') + 
     ylab('value') +
     facet_wrap(~ type)
     
pdf(fp_img_003, width=7, height=7)
r
dev.off()



#                                                   --------------------
# more plots

data = df
data$type = sapply(data$type, cs_rename)
data$from = sapply(data$type, cs_from)
data$to = sapply(data$type, cs_to)

colour_map = c(
     'from Apro' = '#4169e1',
     'from Rpro' = '#ffa500',
     'from RegE' = '#008000',
     'from TranR' = '#ffd700',
     'from RHet' = '#db7093'
     )

s = 
ggplot(data, aes(y=expression, x=to, colour=from)) + 
     geom_hline(yintersect=0, color='grey50') +
     geom_boxplot(outlier.shape=NA) + 
     facet_wrap(~ from, ncol=5) +
     scale_y_continuous(limits=c(-5, 5)) +
     theme_bw() +
     theme(
          axis.text.x  = element_text(
               size = 10,
               angle = 45,
               colour = 'grey50',
               vjust = 1,
               hjust = 1),
          axis.text.y = element_text(size=10, colour = 'grey50'),
          text = element_text(size=10), 
          axis.title = element_text(size=12), 
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_line(colour='grey50'),
          legend.key = element_rect(colour=NA),
          strip.text = element_text(size=10)
     ) +
     # scale_colour_discrete(name='baseline state') +
     scale_color_manual(
          name='baseline state',
          values=colour_map
          ) +
     
     ylab('expression (log2fc)') + 
     xlab('target state') 
     #ggtitle('gene expression by chromatin state transition')

fp_img_004 = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015-07-24_xswitchtype_yexpression.pdf'   
pdf(fp_img_004, width=8, height=3)
s
dev.off()


t = 
ggplot(data, aes(y=expression, x=dist_rel)) + 
     geom_point(size=0.6, alpha=0.2) + 
     facet_wrap(~ type) +
          scale_y_continuous(limits=c(-5, 5)) +
          scale_x_continuous(breaks=c(0.25, 0.5, 0.75)) +
     theme_bw() +
     theme(
          axis.text.x  = element_text(
               size = 10,
               angle = 45,
               colour = 'grey50',
               vjust = 1,
               hjust = 1),
          axis.text.y = element_text(size=10, colour = 'grey50'),
          text = element_text(size=10), 
          axis.title = element_text(size=12), 
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_line(colour='grey50'),
          legend.key = element_rect(colour=NA),
          strip.text = element_text(size=10)
     ) +
     #scale_colour_discrete(name='baseline state') +
     ylab('expression (log2fc)') + 
     xlab('relative distance to start of gene') 

fp_img_005 = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015-07-24_xreldist_yexpression.pdf'   
pdf(fp_img_005, width=7, height=7)
t
dev.off()


