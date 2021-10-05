library(tidyverse)
library(strex)

# read in the BLASTx output
bla <- read.table('../blastx/blastx_pacbio_cdhit_araport.outfmt6', sep = "\t", header = F)
head(bla)
dim(bla)
 
# make a dataframe of transcript ID to Araport gene ID (to araport transcript ID (with bitscore)
transcript <- as.data.frame(bla$V1)
gene_id <- as.data.frame(str_before_first(bla$V2, '\\.'))
bit <- as.data.frame(bla$V12)
df <- cbind(transcript, gene_id, bit)
colnames(df)[1]<-'transcript'
colnames(df)[2]<-'gene_id'
colnames(df)[3]<-'bit'
head(df)

# keep only genes with highest bitscores
high <- df %>% 
  group_by(gene_id) %>% 
  slice_max(bit, with_ties = F) %>% 
  select(-bit) %>% 
  ungroup() %>% 
  arrange(as.numeric(transcript))
head(high)
dim(high)


# write out a list of transcript and gene mappings
write.table(as.data.frame(high), 'transcript_araport_mapping_nonredundant.txt', col.names = F, row.names = F, quote = F, sep = '\t')

# write out separate transcript and gene IDs
tran <- high %>% 
  select(transcript)
write.table(as.data.frame(tran), 'transcripts_from_araport_mapping_nonredundant.txt', col.names = F, row.names = F, quote = F, sep = '\t')
gen <- high %>% 
  select(gene_id)
write.table(as.data.frame(gen), 'genes_from_araport_mapping_nonredundant.txt', col.names = F, row.names = F, quote = F, sep = '\t')




