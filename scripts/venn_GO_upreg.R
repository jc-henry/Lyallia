library(UpSetR)
library(plyr)
library(dplyr)

### First use venn to get the lists of intersections (UpSetR won't do this)
# this will give the venn/upset intersects. 
# however, if just wanting direct overlaps then the intersect(aus25_pcr1_df,aus25_pcr2_df) would be better
# or Reduce(intersect, list(aus25_pcr1_df,aus25_pcr2_df,pcr1_pcr2)) for more than 2

# import lists of significant GO terms from comparisons
# note: not all datasets had significant GO terms
aus25_pcr1_df <- read.table('../deseq2_1.5FC/aus25_vs_pcr1/agrigo_nopadjust_aus25_pcr1_upreg.txt', header = T, sep = '\t') %>% 
  filter(pvalue < 0.05)
aus25_pcr2_df <- read.table('../deseq2_1.5FC/aus25_vs_pcr2/agrigo_nopadjust_aus25_pcr2_upreg.txt', header = T, sep = '\t') %>% 
  filter(pvalue < 0.05)
aus30_pcr1_df <- read.table('../deseq2_1.5FC/aus30_vs_pcr1/agrigo_nopadjust_aus30_pcr1_upreg.txt', header = T, sep = '\t') %>% 
  filter(pvalue < 0.05)
aus30_pcr2_df <- read.table('../deseq2_1.5FC/aus30_vs_pcr2/agrigo_nopadjust_aus30_pcr2_upreg.txt', header = T, sep = '\t') %>% 
  filter(pvalue < 0.05)

# get only the GO terms
aus25_pcr1_col <- aus25_pcr1_df %>% 
  select(GO_acc)
aus25_pcr2_col <- aus25_pcr2_df %>% 
  select(GO_acc)
aus30_pcr1_col <- aus30_pcr1_df %>% 
  select(GO_acc)
aus30_pcr2_col <- aus30_pcr2_df %>% 
  select(GO_acc)

# make vectors of the sets
aus25_pcr1 <- aus25_pcr1_col$GO_acc
aus25_pcr2 <-aus25_pcr2_col$GO_acc
aus30_pcr1 <- aus30_pcr1_col$GO_acc
aus30_pcr2 <- aus30_pcr2_col$GO_acc


# need to set the list up like this to display the set names
library(venn)
input <- list(AUS25_PCR1=aus25_pcr1,AUS25_PCR2=aus25_pcr2,AUS30_PCR1=aus30_pcr1,AUS30_PCR2=aus30_pcr2)

tmp <- venn(input)                  # run venn
isect <- attr(tmp, "intersection")  # get the intersections as a list of lists

# convert the lists to a dataframe and write out
lst_df <- ldply (isect, data.frame)
write.table(lst_df, 'upreg_GO_venn_intersections.txt', col.names = F, row.names = F, quote = F, sep = '\t')


### make a good-looking venn diagram
library("ggvenn")
pdf('upreg_GO_venn.pdf')
ggvenn(input, fill_color = c("#7D9D33", "#DCC949", "#CD8862", "#775B24"),
  stroke_size = 0.5, set_name_size = 4)
dev.off()

### run UpSetR just in case


# add DEG comparisons as colnames for lists
colnames(aus25_pcr1_col)<-"AUS25_PCR1"
colnames(aus25_pcr2_col)<-"AUS25_PCR2"
colnames(aus30_pcr1_col)<-"AUS30_PCR1"
colnames(aus30_pcr2_col)<-"AUS30_PCR2"

listInput <- c(aus30_pcr1_col,aus30_pcr2_col,aus25_pcr1_col,aus25_pcr2_col)

pdf('upreg_GO_upset.pdf')
upset(fromList(listInput), order.by = "freq", matrix.color = "#2F3D70", main.bar.color = "#2F3D70", sets.bar.color = "#7C7189")
dev.off()

