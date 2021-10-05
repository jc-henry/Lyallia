# takes in the starting fasta and renames according to a mapping file/dataframe
# currently using phylotools v0.2.2

# bring in the mapping file as a dataframe
df <- read.table('transcript_araport_mapping_nonredundant.txt', header = F, sep = '\t')
head(df)
class(df)
dim(df)

# run rename.fasta for renaming the fasta
rename.fasta(infile = 'trimmed_header_nonredundant_gene.fasta', df, 'araport_name_cdhit.fasta')
