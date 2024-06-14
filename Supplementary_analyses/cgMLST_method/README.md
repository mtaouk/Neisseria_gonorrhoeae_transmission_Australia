# Applying our cgMLST clustering method to another dataset

490 *N. gonorrhoeae* genomes were made available through [De Silva et
al. 2016](https://10.1016/S1473-3099(16)30157-8) (accessions can be
found in
[SraRunTable.csv](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/cgMLST_method/SraRunTable.csv)).

These genomes have similar level metadata to our study
(within-individual, within-site). We applied the same cgMLST scheme,
allele calling and refining methods as described in our study to assess
the generalizability of this method. After performing allele calling on
the De Silva genomes, we refined the cgMLST schema to 1,506 genes were
present in 95% of isolates and used these genes as the schema. This was
comparable to our refined scheme of 1,495 genes. 1,476 of the genes
retained in the De Silva 95% schema were included in our refined scheme.
All cgMLST result files including the refined schema and alleles called
for each gene at each isolate can be found in
[cgMLST_results](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/tree/main/Supplementary_analyses/cgMLST_method/cgMLST_results).

Here we show the overall distribution of pairwise allelic differences
across the De Silva genomes:

## Load libraries

```         
# Libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(scales)
library(patchwork)
library(ggbreak)
```

## Read inputs and transform data

The input matrix can be downloaded from
[pairwise_distances_desilva_cgMLST.txt](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/cgMLST_method/pairwise_distances_desilva_cgMLST.txt).

The input metadata can be downloaded from
[patient_ids.csv](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/cgMLST_method/patient_ids.csv).

```         
cgMLST_matrix = as.matrix(read.csv("pairwise_distances_desilva_cgMLST.txt", row.names =1, sep = "\t"))

# melt matrix keeping only top of pyramid
cgMLST_melt = data.frame(ID1=rownames(cgMLST_matrix)[row(cgMLST_matrix)[upper.tri(cgMLST_matrix)]], ID2=colnames(cgMLST_matrix)[col(cgMLST_matrix)[upper.tri(cgMLST_matrix)]], dist=cgMLST_matrix[upper.tri(cgMLST_matrix)])

# Read in data that tells us which isolates are the intra-patient ones
cat = read.csv("patient_ids.csv")
cat = subset(cat, select = c("SRR", "pt_id", "new"))

# Merge cgMLST_melt with cat based on ID1 and ID2 columns
cgMLST_melt = cgMLST_melt %>% merge(cat, by.x = "ID1", by.y = "SRR", all.x = F, all.y = F) %>% merge(cat, by.x = "ID2", by.y = "SRR", all.x = F, all.y = F)

# Check for NAs
sum(is.na(cgMLST_melt))
  
```

## Plot the distribution of pairwise allelic differences

```         
data_merged_density = data.frame(table(cgMLST_melt$dist))
data_merged_density$Var1 = as.numeric(as.character(data_merged_density$Var1))

ggplot(data_merged_density, aes(x=Var1, y=Freq)) +
  geom_area(alpha = 0.7, fill = "#2b8cbe", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8,angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab(" ") +
  xlab("Pairwise allelic difference") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE, big.mark = " "))
```
![00001f](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/assets/90819350/9646d68f-4e4b-4a69-afde-ac4299a9895e)

## Plot the distribution of pairwise allelic differences for the calibration isolates (within individual, within site)

```         

# Create a new column called "check" that indicates whether the Patient_no.x and Patient_no.y columns are the same or different and filter based on intrapatient isolates only

cgMLST_melt$check = ifelse(cgMLST_melt$pt_id.x == cgMLST_melt$pt_id.y, "same", "diff")
cgMLST_melt <- filter(cgMLST_melt, cgMLST_melt$check == "same")
cgMLST_melt <- cgMLST_melt %>%
  mutate(type = case_when(new.x == "intrapatient" ~ "Within individuals",
                          TRUE ~ "Within site"))

cgMLST_melt_sum = data.frame(table(cgMLST_melt$dist, cgMLST_melt$type)) %>%
  filter(Freq > 0)
cgMLST_melt_sum$Var1 = as.numeric(as.character(cgMLST_melt_sum$Var1))

ggplot(cgMLST_melt_sum, aes(x=Var1, y=Freq, fill = factor(Var2, level = c("Within site", "Within individuals")), color = factor(Var2, level = c("Within site", "Within individuals")))) +
  geom_area(alpha = 0.7, position = "identity", color = "black") +
  geom_vline(xintercept = 7, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size=0),
        legend.key.size = unit(2.5, 'mm'),
        legend.text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 3)) +
  ylab(" ") +
  xlab("Pairwise allelic difference") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE, big.mark = " ")) +
  scale_fill_manual(values = c("#fdbb84", "#ef6548"), na.translate=FALSE)  +
  scale_x_continuous(limits = c(0, 50))
```
![000039](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/assets/90819350/2f7d8726-b90a-4065-b5ce-e29e2a7eed5a)

## Generate clusters

```         
dist <- as.dist(cgMLST_matrix) 

hc <- hclust(dist, method = "single")

hc_cut <- cutree(hc, h=7)

#convert to dataframe
hc_cut_df <- as.data.frame(hc_cut)
```

The clustering results can be downloaded from
[hc_cut.csv](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/cgMLST_method/hc_cut.csv).
