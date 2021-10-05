library(tidyverse)
library(strex)

# read in the BLASTx output
bla <- read.table('../blastx/blastx_pacbio_cdhit_araport.outfmt6', sep = "\t", header = F)
head(bla)
dim(bla)
 
# make a dataframe of transcript ID to Araport gene ID to araport transcript ID
Name <- bla$V1
gene <- str_before_first(bla$V2, '\\.')

df_ara <- data.frame(Name, gene)
head(df_ara)
tail(df_ara)
dim(df_ara)

write.table(df_ara, 'tx2gene.txt', col.names = F, row.names = F, quote = F, sep = '\t')

