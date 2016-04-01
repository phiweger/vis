library(ggplot2)
library(magrittr)
library(dplyr)

# fp = '/Users/owl/Desktop/data_molpath/tmp/y.csv'
fp = '/Users/owl/Desktop/data_molpath/tmp/harvest_by_gene_2015-07-20.csv'
df = read.table(fp, sep=',', header=F, stringsAsFactors=F)

names(df) = c(
    'gene_name',
    'gene_id',
    'expression',
    'padj',
    'compartment',
    'num_switch',
    'num_switch_rel',
    'num_dmr',
    'num_dmr_rel',
    'num_tfbs',
    'num_tfbs_rel',
    'majority_count',
    'majority_coverage'
     )

# clean = df$majority_count

# for (i in 1:length(clean)){
#     clean[i] = gsub('\\(', '', clean[i])
#     clean[i] = gsub('\\)', '', clean[i])
#     clean[i] = gsub("\\'", '', clean[i])
#     clean[i] = gsub("\\,", '', clean[i])
# }
# df$majority_count = clean

data = df[df$padj < 0.05,]
data = data[!is.na(data$gene_name),]

p = 
ggplot(data, aes(
    x=compartment,
    y=expression
    )) +
    
    geom_boxplot()

fp_out = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015_07_21_xexpression_ycompartment_genebody.pdf'
pdf(fp_out, width=3, height=5)
p
dev.off()


# owls-MacBook-Air:subcompartments owl$ cut -f2,3,4 GSE63525_GM12878_subcompartments.bed | awk '{print $2-$1"\t"$3}' | sort -k2,2 | groupBy -g 2 -c 1 -o sum
# A1  400600000
# A2  581400000
# B1  349400000
# B2  435700000
# B3  855500000
# B4  11000000
# NA  248700000



cs_rename = function(cs_string){
     cs_string %<>% strsplit('_') %>% 
          unlist() %>% 
          '['(c(3,5)) %>% 
          paste(., collapse=' > ') 
     return(cs_string)
}


data = df[df$padj < 0.05,] 
data = data[!is.na(data$gene_name),]

data$majority_coverage = sapply(data$majority_coverage, cs_rename)
data$majority_count = sapply(data$majority_count, cs_rename)

q = 
ggplot(data, aes(
    x=majority_count,
    y=expression
    )) +
    
    geom_boxplot() +
    theme(axis.text.x  = element_text(
    size = 10,
    angle = 45,
    colour = "black",
    vjust = 1,
    hjust = 1)
    )


fp_out = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015_07_21_xexpression_ymajorvotecount_genebody.pdf'
pdf(fp_out, width=7, height=5)
q
dev.off()

fp_out = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015_07_21_xexpression_ymajorvotecoverage_genebody.pdf'
pdf(fp_out, width=7, height=5)
r
dev.off()



s = 
ggplot(data, aes(
    x=compartment,
    y=expression,
    )) +
    
    geom_boxplot() +
    facet_wrap(~ majority_count) 

    # theme(axis.text.x  = element_text(
    # size = 10,
    # angle = 45,
    # colour = "black",
    # vjust = 1,
    # hjust = 1)
    # )

fp_out = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015_07_21_xexpression_ycompartment_facetmajorvotecount_genebody.pdf'
pdf(fp_out, width=9, height=7)
s
dev.off()



tmp = data[!is.na(data$compartment != 'NA'),]
comp = tmp$compartment

for (i in 1:length(comp)){
    if (grepl('A', comp[i])){
        comp[i] = 'A'
    } else {
        comp[i] = 'B'
    }
}

tmp$comp = comp

t = 
ggplot(tmp, aes(
    x=comp,
    y=expression,
    )) +
    
    geom_boxplot() +
    facet_wrap(~ majority_count) 

fp_out = '/Users/owl/Desktop/on_genomics/lab_journal_img/img_2015-07_july/2015_07_21_xexpression_ycompartment_facetmajorvotecount_pureAB_genebody.pdf'
pdf(fp_out, width=5, height=7)
t
dev.off()

