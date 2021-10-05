library(tidyverse)

# read in the background GO annotations
go <- read.table('../../blastx/background_genes_GO.txt', header = F, sep = '\t')
head(go)
dim(go)

# read in all DEGs
deg <- read.table('geneIDonly_diffex_aus25_pcr1_padj0.5_1.5fc.txt', header = FALSE, sep = '\t')
head(deg)
dim(deg)

# join the background Go annotations to the DEG names 
# remove any gene name without GO terms
merg_deg <- left_join(deg, go, by = 'V1') %>%
  drop_na()
head(merg_deg)
dim(merg_deg) 

write.table(merg_deg,'GO_terms_degs_aus25_pcr1.txt', col.names = F, row.names = F, quote = F, sep = '\t')

# read in upregulated genes
up <- read.table('upreg_geneIDonly_aus25_pcr1_padj0.5_1.5fc.txt', header = FALSE, sep = '\t')
head(up)
dim(up)

merg_up <- left_join(up, go, by = 'V1') %>%
  drop_na()
head(merg_up)
dim(merg_up) 

write.table(merg_up,'GO_terms_upreg_aus25_pcr1.txt', col.names = F, row.names = F, quote = F, sep = '\t')

# read in downregulated genes
down <- read.table('downreg_geneIDonly_aus25_pcr1_padj0.5_1.5fc.txt', header = FALSE, sep = '\t')
head(down)
dim(down)

merg_down <- left_join(down, go, by = 'V1') %>%
  drop_na()
head(merg_down)
dim(merg_down) 

write.table(merg_down,'GO_terms_downreg_aus25_pcr1.txt', col.names = F, row.names = F, quote = F, sep = '\t')