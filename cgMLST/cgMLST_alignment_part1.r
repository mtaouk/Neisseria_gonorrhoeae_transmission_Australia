#!/usr/bin/env Rscript

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/NG_transmission/cgMLST_alignment")

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Read in cgMLST file
mydata <- fread("cgMLST.tsv", header=T, sep="\t")

# Remove suffix from sample names 
mydata$FILE <- gsub(".contigs.fa", "", mydata$FILE)

# Convert to long format
mydata2 <- mydata %>% pivot_longer(starts_with("NEI"), values_to="group", names_to="allele")

# Remove suffix from allele names
mydata2$allele <- gsub(".fasta", "", mydata2$allele)

# Add new column by merging gene name with group- for pulling out sequences by header name
# Then arrange by allele for faster splitting later
finaldf <- mydata2 %>% 
		mutate(newname=paste(allele,"_", group, sep="")) %>% 
		select(FILE, allele, newname) %>% 
		arrange(allele, FILE)

# Write to new file (large file)
fwrite(finaldf, file="cgMLST_transposed.tsv", col.names=F, row.names=F, sep="\t", quote=F)