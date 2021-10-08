library(dplyr)

# read in all the full GO terms
a <- read.table('../deseq2_1.5FC/aus25_vs_pcr1/agrigo_nopadjust_aus25_pcr1_upreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(a)
dim(a)
b <- read.table('../deseq2_1.5FC/aus25_vs_pcr2/agrigo_nopadjust_aus25_pcr2_upreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(b)
dim(b)
c <- read.table('../deseq2_1.5FC/aus30_vs_pcr1/agrigo_nopadjust_aus30_pcr1_upreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(c)
dim(c)
d <- read.table('../deseq2_1.5FC/aus30_vs_pcr2/agrigo_nopadjust_aus30_pcr2_upreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(d)
dim(d)
e <- read.table('../deseq2_1.5FC/aus25_vs_pcr1/agrigo_nopadjust_aus25_pcr1_downreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(e)
dim(e)
f <- read.table('../deseq2_1.5FC/aus25_vs_pcr2/agrigo_nopadjust_aus25_pcr2_downreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(f)
dim(f)
g <- read.table('../deseq2_1.5FC/aus30_vs_pcr1/aaagrigo_nopadjust_aus30_pcr1_downreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(g)
dim(g)
h <- read.table('../deseq2_1.5FC/aus30_vs_pcr2/agrigo_nopadjust_aus30_pcr2_downreg.txt', header = T, sep = '\t') %>% 
  select(c(GO_acc,term_type,Term))
head(h)
dim(h)

# row bind all the terms
df <- bind_rows(a,b,c,d,e,f,g,h)
df
dim(df)

# keep only unique terms
un_df <- df %>% 
  unique()
head(un_df)
dim(un_df)

write.table(un_df,'test3.txt', col.names = T, row.names = F, quote = F, sep = '\t')

# bring in the lists of intersecting GO terms
up <- read.table('upreg_GO_venn_intersections.txt', header = F, sep = '\t')
head(up)
dim(up)

# bind the DFs on the GO term IDs
up_df <- left_join(up, un_df, by=c('V2'='GO_acc')) %>% 
  rename(intersection=V1) %>% 
  rename(GO_id=V2)
head(up_df)
dim(up_df)

dn <- read.table('downreg_GO_venn_intersections.txt', header = F, sep = '\t')
head(dn)
dim(dn)

dn_df <- left_join(dn, un_df, by=c('V2'='GO_acc')) %>% 
  rename(intersection=V1) %>% 
  rename(GO_id=V2)
head(dn_df)
dim(dn_df)

write.table(up_df, 'upreg_GO_venn_intersections_terms_annotated.txt', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(dn_df, 'downreg_GO_venn_intersections_terms_annotated.txt', col.names = T, row.names = F, quote = F, sep = '\t')
