---
title: "desilva"
author: "Mona Taouk"
date: "2023-07-10"
output: html_document
---

```{r message=FALSE, warning=FALSE}
# Libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(scales)
library(patchwork)
library(ggbreak)

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/NG_transmission/deSilva/")
```

```{r message=FALSE, warning=FALSE}
# Read in pairwise_distances_0.95_cgMLST.txt and convert to a matrix
cgMLST_matrix = as.matrix(read.csv("pairwise_distances_desilva_cgMLST.txt", row.names =1, sep = "\t"))

# melt matrix keeping only top of pyramid
cgMLST_melt = data.frame(ID1=rownames(cgMLST_matrix)[row(cgMLST_matrix)[upper.tri(cgMLST_matrix)]], ID2=colnames(cgMLST_matrix)[col(cgMLST_matrix)[upper.tri(cgMLST_matrix)]], dist=cgMLST_matrix[upper.tri(cgMLST_matrix)])

# Read in intrapatient.csv
cat = read.csv("patient_ids.csv")
cat = subset(cat, select = c("SRR", "pt_id", "new"))

# Merge cgMLST_melt with cat based on ID1 and ID2 columns
cgMLST_melt = cgMLST_melt %>% merge(cat, by.x = "ID1", by.y = "SRR", all.x = F, all.y = F) %>% merge(cat, by.x = "ID2", by.y = "SRR", all.x = F, all.y = F)

# Check for NAs
sum(is.na(cgMLST_melt))

# Create a new column called "check" that indicates whether the Patient_no.x and Patient_no.y columns are the same or different
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

```{r}
cgMLST_melt_all = data.frame(ID1=rownames(cgMLST_matrix)[row(cgMLST_matrix)[upper.tri(cgMLST_matrix)]], ID2=colnames(cgMLST_matrix)[col(cgMLST_matrix)[upper.tri(cgMLST_matrix)]], dist=cgMLST_matrix[upper.tri(cgMLST_matrix)])

data_merged_density = data.frame(table(cgMLST_melt_all$dist))
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

```{r}
# threshold

#without-threshold
dist <- as.dist(cgMLST_matrix) 

hc <- hclust(dist, method = "single")

hc_cut <- cutree(hc, h=7)

#convert to dataframe
hc_cut_df <- as.data.frame(hc_cut)
```
