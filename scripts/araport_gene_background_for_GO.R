library(tidyverse)
library(strex)

# read in the BLASTx output
bla <- read.table('../blastx/blastx_pacbio_cdhit_araport.outfmt6', sep = "\t", header = F)
head(bla)
dim(bla)
 
# make a dataframe of Araport gene ID 
araport_gene <- as.data.frame(str_before_first(bla$V2, '\\.'))
head(araport_gene)
dim(araport_gene)

df <- distinct(araport_gene)
head(df)
dim(df)

write.table(df, 'araport_gene_background_for_GO.txt', col.names = F, row.names = F, quote = F, sep = '\t')  


