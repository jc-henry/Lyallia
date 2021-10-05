library(tidyverse)
library(readxl)
library(reshape)

# reading in GO terms from uniprot
# used the background gene list (all unique Araport genes)
tair <- read.csv('uniprot-yourlist_20210922.txt', header = TRUE, sep = '\t')
dim(tair)

# keep only the Arabidopsis IDs with the most GO terms (in case of duplicates)
tair_len <- tair %>%
  subset(select = c(yourlist.M20210921F248CABF64506F29A91F8037F07B67D11F9A24E, Gene.ontology.IDs)) %>% 
  mutate(golen=str_length(Gene.ontology.IDs)) %>% # make a new column of string (go term) lengths
  group_by(yourlist.M20210921F248CABF64506F29A91F8037F07B67D11F9A24E) %>%              # group by gene ID
  slice_max(golen) %>%                # keep only the ID with the most GO terms
  filter(golen != 0) %>%              # drop any row with no GO terms
  distinct() %>%                      # retain only unique rows
  select(-golen)                      # drop the length column
dim(tair_len)

# write this out and separate the go_id in excel!
write.table(tair_len, file = 'GO_linear_out.txt', col.names = F, row.names = F, quote = F,  sep = '\t')

# bring back in as tab separated
in_tair <- read.csv('GO_linear_in.txt', header = F, sep = '\t')
head(in_tair)
dim(in_tair)

# use the Melt function to move all the column values to rows.
# retains the information as to which column the GO terms came from.
go_df_fil <- melt(in_tair, id.vars = c('V1')) %>%
  select(-variable) %>%
  mutate(golen=str_length(value)) %>% # make a new column of string (go term) lengths (zero when no mapping)
  filter(golen != 0) %>%             # drop any row with no GO terms
  select(-golen)
head(go_df_fil)
dim(go_df_fil)

# write out the background with GO annotations
write.table(go_df_fil, 'background_genes_GO.txt', col.names = F, row.names = F, quote = F,  sep = '\t')


