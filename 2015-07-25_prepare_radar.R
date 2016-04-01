prepare_radar = function(df, padj){

df = df[df$padj < padj,]
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



mock$type = sapply(mock$type, cs_rename)
smry$type = sapply(smry$type, cs_rename)

j = inner_join(mock, smry)

return(j)
}