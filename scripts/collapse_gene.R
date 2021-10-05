library(tidyverse)
library(gplots)
library(matrixTests)
library(strex)

# read in the BLASTx output
bla <- read.table('../blastx/blastx_pacbio_cdhit_araport.outfmt6', sep = "\t", header = F)
head(bla)
dim(bla)
 
# make a dataframe of transcript ID to Araport gene ID to araport transcript ID
Name <- bla$V1
araport_gene <- str_before_first(bla$V2, '\\.')

df_ara <- data.frame(Name, araport_gene)
head(df_ara)
dim(df_ara)
  
# make objects of the TPM values from quant.sf files
# keep only transcript ID and TPM column renamed as sample ID
aus25_1 <- read.table('../salmon/aus25_1.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus25_1 = TPM) %>% 
  as.data.frame()
head(aus25_1)  
 
aus25_3 <- read.table('../salmon/aus25_3.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus25_3 = TPM) %>% 
  as.data.frame()
aus25_4 <- read.table('../salmon/aus25_4.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus25_4 = TPM) %>% 
  as.data.frame()
aus25_5 <- read.table('../salmon/aus25_5.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus25_5 = TPM) %>% 
  as.data.frame()
aus25_7 <- read.table('../salmon/aus25_7.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus25_7 = TPM) %>% 
  as.data.frame()
aus25_9 <- read.table('../salmon/aus25_9.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus25_9 = TPM) %>% 
  as.data.frame()
aus30_0 <- read.table('../salmon/aus30_0.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus30_0 = TPM) %>% 
  as.data.frame()
aus30_2 <- read.table('../salmon/aus30_2.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus30_2 = TPM) %>% 
  as.data.frame()
aus30_3 <- read.table('../salmon/aus30_3.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus30_3 = TPM) %>% 
  as.data.frame()
aus30_5 <- read.table('../salmon/aus30_5.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(aus30_5 = TPM) %>% 
  as.data.frame()
mac1_6 <- read.table('../salmon/mac1_6.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac1_6 = TPM) %>% 
  as.data.frame()
mac1_7 <- read.table('../salmon/mac1_7.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac1_7 = TPM) %>% 
  as.data.frame()
mac3_2 <- read.table('../salmon/mac3_2.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac3_2 = TPM) %>% 
  as.data.frame()
mac3_5 <- read.table('../salmon/mac3_5.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac3_5 = TPM) %>% 
  as.data.frame()
mac3_7 <- read.table('../salmon/mac3_7.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac3_7 = TPM) %>% 
  as.data.frame()
mac3_9 <- read.table('../salmon/mac3_9.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac3_9 = TPM) %>% 
  as.data.frame()
mac3_10 <- read.table('../salmon/mac3_10.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(mac3_10 = TPM) %>% 
  as.data.frame()
pcr1_4 <- read.table('../salmon/pcr1_4.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr1_4 = TPM) %>% 
  as.data.frame()
pcr1_5 <- read.table('../salmon/pcr1_5.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr1_5 = TPM) %>% 
  as.data.frame()
pcr1_7 <- read.table('../salmon/pcr1_7.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr1_7 = TPM) %>% 
  as.data.frame()
pcr1_9 <- read.table('../salmon/pcr1_9.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr1_9 = TPM) %>% 
  as.data.frame()
pcr1_10 <- read.table('../salmon/pcr1_10.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr1_10 = TPM) %>% 
  as.data.frame()
pcr2_3 <- read.table('../salmon/pcr2_3.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr2_3 = TPM) %>% 
  as.data.frame()
pcr2_4 <- read.table('../salmon/pcr2_4.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr2_4 = TPM) %>% 
  as.data.frame()
pcr2_8 <- read.table('../salmon/pcr2_8.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr2_8 = TPM) %>% 
  as.data.frame()
pcr2_10 <- read.table('../salmon/pcr2_10.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr2_10 = TPM) %>% 
  as.data.frame()
pcr2_18 <- read.table('../salmon/pcr2_18.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pcr2_18 = TPM) %>% 
  as.data.frame()
pdc1_2 <- read.table('../salmon/pdc1_2.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc1_2 = TPM) %>% 
  as.data.frame()
pdc1_4 <- read.table('../salmon/pdc1_4.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc1_4 = TPM) %>% 
  as.data.frame()
pdc1_5 <- read.table('../salmon/pdc1_5.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc1_5 = TPM) %>% 
  as.data.frame()
pdc1_7 <- read.table('../salmon/pdc1_7.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc1_7 = TPM) %>% 
  as.data.frame()
pdc1_8 <- read.table('../salmon/pdc1_8.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc1_8 = TPM) %>% 
  as.data.frame()
pdc3_1 <- read.table('../salmon/pdc3_1.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc3_1 = TPM) %>% 
  as.data.frame()
pdc3_4 <- read.table('../salmon/pdc3_4.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc3_4 = TPM) %>% 
  as.data.frame()
pdc3_8 <- read.table('../salmon/pdc3_8.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc3_8 = TPM) %>% 
  as.data.frame()
pdc3_10 <- read.table('../salmon/pdc3_10.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pdc3_10 = TPM) %>% 
  as.data.frame()
pjda6_1 <- read.table('../salmon/pjda6_1.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda6_1 = TPM) %>% 
  as.data.frame()
pjda6_2 <- read.table('../salmon/pjda6_2.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda6_2 = TPM) %>% 
  as.data.frame()
pjda6_3 <- read.table('../salmon/pjda6_3.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda6_3 = TPM) %>% 
  as.data.frame()
pjda6_6 <- read.table('../salmon/pjda6_6.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda6_6 = TPM) %>% 
  as.data.frame()
pjda6_7 <- read.table('../salmon/pjda6_7.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda6_7 = TPM) %>% 
  as.data.frame()
pjda6_9 <- read.table('../salmon/pjda6_9.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda6_9 = TPM) %>% 
  as.data.frame()
pjda10_1 <- read.table('../salmon/pjda10_1.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda10_1 = TPM) %>% 
  as.data.frame()
pjda10_5 <- read.table('../salmon/pjda10_5.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda10_5 = TPM) %>% 
  as.data.frame()
pjda10_8 <- read.table('../salmon/pjda10_8.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda10_8 = TPM) %>% 
  as.data.frame()
pjda10_21 <- read.table('../salmon/pjda10_21.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(pjda10_21 = TPM) %>% 
  as.data.frame()
rba5_2 <- read.table('../salmon/rba5_2.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(rba5_2 = TPM) %>% 
  as.data.frame()
rba5_3 <- read.table('../salmon/rba5_3.out/quant.sf', sep = "\t", header = TRUE) %>% 
  subset(select = c(Name, TPM)) %>% 
  left_join(df_ara, by='Name') %>% 
  group_by(araport_gene) %>% 
  summarise('TPM'=sum(TPM)) %>% 
  rename(rba5_3 = TPM) %>% 
  as.data.frame()

# make the dataframe with all TPM values against the renamed genes
# also remove the NA for araport ID (failed blast matches)
df_tpm <- Reduce(merge,list(aus25_1,aus25_3,aus25_4,aus25_5,aus25_7,aus25_9,
                            aus30_0,aus30_2,aus30_3,aus30_5,
                            mac1_6,mac1_7,
                            mac3_2,mac3_5,mac3_7,mac3_9,mac3_10,
                            pcr1_4,pcr1_5,pcr1_7,pcr1_9,pcr1_10,
                            pcr2_3,pcr2_4,pcr2_8,pcr2_10,pcr2_18,
                            pdc1_2,pdc1_4,pdc1_5,pdc1_7,pdc1_8,
                            pdc3_1,pdc3_4,pdc3_8,pdc3_10,
                            pjda6_1,pjda6_2,pjda6_3,pjda6_6,pjda6_7,pjda6_9,
                            pjda10_1,pjda10_5,pjda10_8,pjda10_21,
                            rba5_2,rba5_3)) %>% 
  drop_na()

head(df_tpm)
dim(df_tpm)

# export the file
write.table(df_tpm, 'all_transcript__including_zeros_TPM.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

# remove any transcript with no reads mapping
df_nozero <- df_tpm[rowSums(df_tpm[, -(1)]) > 0, ]

# export the file
write.table(df_nozero, 'all_transcript_no_zero_TPM.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

# create a vector of gene(cluster) names
gene_id <- as.data.frame(df_nozero$araport_gene)

# convert dataframe to a numeric matrix
mat_no_id <- data.matrix(df_nozero[,2:49])

# perform a Kruskal-Wallis rank sum test on each row of the matrix (group=taxa)
kw_group <- row_kruskalwallis(mat_no_id[,1:48], c('aus25','aus25','aus25','aus25','aus25','aus25',
                                                  'aus30','aus30','aus30','aus30',
                                                  'mac1','mac1',
                                                  'mac3','mac3','mac3','mac3','mac3',
                                                  'pcr1','pcr1','pcr1','pcr1','pcr1',
                                                  'pcr2','pcr2','pcr2','pcr2','pcr2',
                                                  'pdc1','pdc1','pdc1','pdc1','pdc1',
                                                  'pdc3','pdc3','pdc3','pdc3',
                                                  'pjda6','pjda6','pjda6','pjda6','pjda6','pjda6',
                                                  'pjda10','pjda10','pjda10','pjda10',
                                                  'rba5','rba5'))

# merge the two dataframes
mer_df <- cbind(gene_id,  kw_group)

# rename the gene id column
mer_df_rename <- mer_df %>%
  rename(Araport = 'df_nozero$araport_gene')

# extract genes with pvalue < 0.05
less_05 <- subset(mer_df_rename, mer_df_rename[,6] < 0.05)

# remove all cols except Gene id
sig <- subset(less_05, select = -c(obs.tot, obs.groups, df, statistic, pvalue)) 

# export the file
write.table(sig, 'sig_diff_araports.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')





# make a TPM object with the "_i" attached (needed for grep to be specific since no -x flag in R) 
sig_p <- sig %>%
  mutate(Hold = paste(Araport, '_i', sep = ''))

# need to do the same as above for the input df_nozero
df_nozero_p <- df_nozero %>%
  mutate(hold = paste(araport_gene, '_i', sep = ''))

# grep and merge the files
clus_idx2 <- sapply(sig_p$Hold, grep, df_nozero_p$hold) # create an index
clus_idx1 <- sapply(seq_along(clus_idx2), function(i) rep(i, length(clus_idx2[[i]]))) # duplicate the original indices so the results align
clus_merge <- cbind(sig[unlist(clus_idx1),,drop=F], df_nozero_p[unlist(clus_idx2),,drop=F]) # merge the datasets with cbind aligned on the new indices


# remove unnecessary "gene" column
# remove all cols except Gene id
fin_df <- subset(clus_merge, select = -c(araport_gene, hold))
head(fin_df)

# write out the final df if case needed later
write.table(fin_df, 'cleaned_final_sigdiff.txt', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

### heatmap of all sigdiff genes
# setup with log scaling (necessary for making heatmap).
copy_df <- fin_df
rownames(copy_df) <- NULL
name_df <- column_to_rownames(copy_df, var = "Araport")
fin_mat <- data.matrix(name_df)
Log_LC <- log10(fin_mat+1)
SC_LC <- scale(Log_LC,scale = TRUE)
y <- data.matrix(SC_LC)

## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete") # row clustering (t transposes the data since cor correlates on columns)
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") # column clustering


pdf("sigdiff_gene_heatmap.pdf", width = 8, height = 6)
mycol <- colorpanel(100, "green", "black", "red")

heatmap.2(
  y,
  Rowv = as.dendrogram(hr),
  Colv = as.dendrogram(hc),
  col = mycol,
  density.info = "none",
  trace = "none",
  dendrogram = "both",
  scale = "row",
  labRow = FALSE,
  labCol = NULL,
  srtCol=45,  adjCol = c(1,1),
  cexCol = 0.5,
  margins = c(10,7)
)
dev.off()

