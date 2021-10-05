library(UpSetR)
library(venn)
library(plyr)

### First use venn to get the lists of intersections (UpSetR won't do this)
# this will give the venn/upset intersects. 
# however, if just wanting direct overlaps then the intersect(aus25_pcr1,aus25_pcr2) would be better
# or Reduce(intersect, list(aus25_pcr1,aus25_pcr2,aus30_aus25,pcr1_pcr2)) for more than 2

# import lists of DEGs from comparisons
aus25_pcr1 <- read.table('../deseq2_1.5FC/aus25_vs_pcr1/geneIDonly_diffex_aus25_pcr1_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus25_pcr2 <- read.table('../deseq2_1.5FC/aus25_vs_pcr2/geneIDonly_diffex_aus25_pcr2_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus30_aus25 <- read.table('../deseq2_1.5FC/aus30_vs_aus25/geneIDonly_diffex_aus30_aus25_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus30_pcr1 <- read.table('../deseq2_1.5FC/aus30_vs_pcr1/geneIDonly_diffex_aus30_pcr1_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus30_pcr2 <- read.table('../deseq2_1.5FC/aus30_vs_pcr2/geneIDonly_diffex_aus30_pcr2_padj0.5_1.5fc.txt', header = F, sep = '\t')
pcr1_pcr2 <- read.table('../deseq2_1.5FC/pcr1_vs_pcr2/geneIDonly_diffex_pcr1_pcr2_padj0.5_1.5fc.txt', header = F, sep = '\t')

# make vectors of the sets
aus25_pcr1 <- aus25_pcr1$V1
aus25_pcr2 <-aus25_pcr2$V1
aus30_aus25 <- aus30_aus25$V1
aus30_pcr1 <- aus30_pcr1$V1
aus30_pcr2 <- aus30_pcr2$V1
pcr1_pcr2 <- pcr1_pcr2$V1

# need to set the list up like this to display the set names
input <- list(aus25_pcr1=aus25_pcr1,aus25_pcr2=aus25_pcr2,aus30_aus25=aus30_aus25,
                aus30_pcr1=aus30_pcr1,aus30_pcr2=aus30_pcr2,pcr1_pcr2=pcr1_pcr2)

tmp <- venn(input)                  # run venn
isect <- attr(tmp, "intersection")  # get the intersections as a list of lists

# convert the lists to a dataframe and write out
lst_df <- ldply (isect, data.frame)
write.table(lst_df, 'venn_intersections.txt', col.names = F, row.names = F, quote = F, sep = '\t')


### run UpSetR for the visuals (6 sets too many for venn diagram)

# re-import lists of DEGs from comparisons (been altered from the venn intersections)
aus25_pcr1 <- read.table('../deseq2_1.5FC/aus25_vs_pcr1/geneIDonly_diffex_aus25_pcr1_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus25_pcr2 <- read.table('../deseq2_1.5FC/aus25_vs_pcr2/geneIDonly_diffex_aus25_pcr2_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus30_aus25 <- read.table('../deseq2_1.5FC/aus30_vs_aus25/geneIDonly_diffex_aus30_aus25_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus30_pcr1 <- read.table('../deseq2_1.5FC/aus30_vs_pcr1/geneIDonly_diffex_aus30_pcr1_padj0.5_1.5fc.txt', header = F, sep = '\t')
aus30_pcr2 <- read.table('../deseq2_1.5FC/aus30_vs_pcr2/geneIDonly_diffex_aus30_pcr2_padj0.5_1.5fc.txt', header = F, sep = '\t')
pcr1_pcr2 <- read.table('../deseq2_1.5FC/pcr1_vs_pcr2/geneIDonly_diffex_pcr1_pcr2_padj0.5_1.5fc.txt', header = F, sep = '\t')

# add DEG comparisons as colnames for lists
colnames(aus25_pcr1)<-"AUS25_PCR1"
colnames(aus25_pcr2)<-"AUS25_PCR2"
colnames(aus30_aus25)<-"AUS30_AUS25"
colnames(aus30_pcr1)<-"AUS30_PCR1"
colnames(aus30_pcr2)<-"AUS30_PCR2"
colnames(pcr1_pcr2)<-"PCR1_PCR2"

listInput <- c(aus25_pcr1,aus25_pcr2,aus30_aus25,aus30_pcr1,aus30_pcr2,pcr1_pcr2)

pdf('upset_gene_1.5FC.pdf')
upset(fromList(listInput), nsets=6, order.by = "freq", matrix.color = "#2F3D70", main.bar.color = "#2F3D70", sets.bar.color = "#7C7189")
dev.off()

