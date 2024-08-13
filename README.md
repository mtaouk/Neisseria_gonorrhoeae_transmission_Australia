# Longitudinal Genomic Analysis of *Neisseria gonorrhoeae* Transmission Dynamics in Australia

Mona L. Taouk, George Taiaroa, Sebastian Duchene, Soo Jen Low, Charlie
K. Higgs, Darren Y. J. Lee, Shivani Pasricha, Nasra Higgins, Danielle J.
Ingle, Benjamin P. Howden, Marcus Y. Chen, Christopher K. Fairley, Eric
P. F. Chow, Deborah A. Williamson

This GitHub contains all code used in this study. All reads used in this
study are available in the NCBI database under BioProjects PRJNA520805
and PRJNA1033374. Accessions are available in the Supplementary Dataset
published with this study. All analyses performed for this study can be
replicated using the code below, however input names, directory paths
and file names may need to be edited accordingly.

## Quality check reads

### 1. Kraken

Trimmed paired end reads for each genome were input into Kraken2
(v2.1.1) to investigate for species contamination:

```         
for i in $(cat IDs.txt); do kraken2 --db /home/linuxbrew/db/kraken2/pluspf --gzip-compressed --paired /home/taouk/NGtransmission/reads/${i}_1.fq.gz /home/taouk/NGtransmission/reads/${i}_2.fq.gz --report ${i}.kraken.txt; done
```

The top species match and percentage were extracted from each result
file and summarised:

```         
grep -m1 -P "\tS\t" *.txt | sort -k2nr | sed 's/.txt: /\t/' > kraken_results.txt
```

### 2. Sequencing depth

Trimmed paired end reads for each genome were aligned to the NCCP11945
reference genome (GenBank accession NC_011035.1) using Minimap2
(v2.17-r941) and Samtools (v1.10, using HTSlib v1.10.2) to calculate the
number of reads aligning to the reference genome:

```         
minimap2 -t10 -ax sr NCCP11945.fa AUSMDU00008753_1.fq.gz AUSMDU00008753_2.fq.gz | samtools fastq -f 2 - | grep -c "^@" > AUSMDU00008753_count.txt
```

The outputs were combined:

```         
grep "" *_count.txt | sed 's/_count.txt:/\t/' > test.txt
```

The average read length for each genome was calculated using fq
(v0.11.0) and extracted from the results:

```         
fq --ref NCCP11945.fa AUSMDU00008753_1.fq.gz AUSMDU00008753_2.fq.gz > AUSMDU00008753.yield.tab
```

```         
grep -m1 -P "AvgLen" *.yield.tab | sort -k2nr > yield.all3.txt
```

The average read depth was calculated as follows: (number of reads
mapping to reference \* average read length)/reference genome length

## Assemblies

### 1. Shovill

*De novo* genome assemblies were generated using Shovill (v1.1.0):

```         
shovill --outdir assemblies/AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz
```

The names of each genome's contigs were changed from the default:

```         
for i in \$(cat IDs.txt); do mv \${i}/contigs.fa \${i}/\${i}.contigs.fa;
done
```

### 2. Number of contigs

The number of contigs for each genome was calculated:

```         
grep -c "^>" */*contigs.fa > contigs_all.txt
```

## Genotyping

### 1. MLST

MLSTs were identified using mlst (v2.23.0; PubMLST database accessed May
5th 2023):

```         
mlst --legacy --scheme neisseria assemblies/*/*.contigs.fa > mlst_results.txt
```

### 2. NG-MAST

NG‐MASTs were assigned using NGmaster (1.0.0; NG-MAST v2.0 PubMLST
database accessed May 5th, 2023):

```         
ngmaster assemblies/*/*.contigs.fa > ngmast_results.txt
```

### 3. NG-STAR

The NG-STAR typing scheme was used to determine the phenotypic profile
of seven resistance genes (*penA*, *mtrR*, *porB*, *ponA*, *gyrA*,
*parC* and *23S* rRNA), using pyngSTar(database version 2.0). NG-STAR
alleles and profiles obtained from the database hosted by the
<a href="https://ngstar.canada.ca/alleles/loci_selection?lang=en" title="Public Health Agency of Canada">Public
Health Agency of Canada</a>:

```         
python3 /home/taouk/pyngSTar/pyngSTar.py -f -a -i assemblies/*/*.contigs.fa -p /home/taouk/pyngSTar/pyngSTarDB/ > ngstar_results.txt
```

## cgMLST

### 1. Prepare schema

It is recommended to prepare the schema for each new dataset.

A prodigal training file for *N. gonorrhoeae* was made using the
NCCP11945 reference genome using prodigal (v.2.6.3):

```         
prodigal -i NCCP11945.fa -t neisseria_gonorrhoeae.trn -p single
```

Training file:
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/neisseria_gonorrhoeae.trn">cgMLST/neisseria_gonorrhoeae.trn</a>.

The *N. gonorrhoeae* cgMLST schema developed by
<a href="https://doi: 10.1093/infdis/jiaa002" title="Harrison et al.">Harrison
et al.</a> was downloaded from
<a href="https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=schemeInfo&scheme_id=62" title="PubMLST">PubMLST</a>
and prepared for use using chewBBACA (v2.8.5):

```         
chewBBACA.py PrepExternalSchema -i scheme_directory -o scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn --cpu 50 --st 0.3
```

### 2. Allele calling

Allele calling was performed on all *N. gonorrhoeae* genomes, using the
assemblies as input:

```         
chewBBACA.py AlleleCall -i paths_to_assemblies.txt -g scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn -o results --cpu 50
```

### 3. Refining schema

The schema was refined to only include genes present in 95% of genomes
in the dataset:

```         
chewBBACA.py ExtractCgMLST -i results/results_20211011T232621/results_alleles.tsv --r results/results_20220602T014232/RepeatedLoci.txt -o evaluate95 --t 0.95
```

cgMLST allele calling final results:

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLST.tsv">cgMLST/Results/cgMLST.tsv</a>:
Table of alleles for each gene in the scheme for all isolates.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLSTschema.txt">cgMLST/Results/cgMLSTschema.txt</a>:
1495 genes included in final schema for this study.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/mdata_stats.tsv">cgMLST/Results/mdata_stats.tsv</a>:
Number of and percentage of uncalled genes per genome in the dataset.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/presence_absence.csv">cgMLST/Results/presence_absence.csv</a>:
Gene presence/absence table. (1 present, 0 absent).

### 4. cgMLST allelic differences

The cgMLST results table was transformed to a symmetrical distance
matrix using cgmlst-dists:

```         
cgmlst-dists results/evaluate/cgMLST.tsv > cgMLST_matrix.txt
```

Resulting matrix:
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLST_matrix.csv.zip">cgMLST/Results/cgMLST_matrix.csv.zip</a>.

The values from this matrix were adjusted in
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia?tab=readme-ov-file#12-determine-threshold-for-clustering">step
12</a>.

### 5. cgMLST alignment (only including 5,881 genomes)

R Script to preprocess the data:
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/cgMLST_alignment_part1.r">cgMLST/cgMLST_alignment_part1.r</a>
This script takes the cgMLST.tsv result file and converts it into a long
transposed version that lists the allele called at each gene for each
isolate.

Shell script to make the alignments:
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/cgMLST_alignment_part2.r">cgMLST/cgMLST_alignment_part2.r</a>
This script takes the transposed output file from the above R script and
uses this as input. This script also requires a second input
(all_alleles_cat.fasta), which is a fasta file that contains all the
alleles for all the genes in the cgMLST database. This can be generated
by concatenating the fasta files containing all alleles for each gene
together. This file will be \~1.7GB.

### 6. 95% core alignment (only including 5,881 genomes)

The alignment of 5,881 genomes representing a unique infection each, was
processed to create a 95% soft core using trimAl (v1.4.rev15):

```         
trimal -in unique_pass.aln -out stripped_alignment.fasta -gt 0.05 -threads 40
```

### 10. Maximum likelihood phylogeny

A maximum likelihood phylogenetic tree was generated using IQ-tree
(v2.0.3):

```         
iqtree -s stripped_alignment.fasta -B 1000 -T 40
```

The best-fitting nucleotide substitution model based on the lowest
Bayesian Information Criterion (BIC) was GTR+I+R10.

### 11. Timed phylogeny

Molecular dating of ancestral events was performed using the
least-squares dating (LSD) software (v0.3), with the maximum likelihood
phylogeny generated in the above step used as input:

```         
/home/taouk/lsd-0.3beta-master/src/lsd -d dates.tsv -i 95gapsMP.tree -c -r a
```

### 12. Determine threshold for clustering

The following R code details steps to determine a threshold for
clustering, adjusting the threshold for time between sample collection
dates and performing clustering:
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Threshold.Rmd">cgMLST/Threshold.Rmd</a>.

The data required for the code can be downloaded from
<a href="https://doi.org/10.26188/25989001">doi.org/10.26188/25989001</a>.

## Odds ratios

To assess variables associated with cluster persistence, we applied a
multivariable logistic regression model via generalised estimating
equation (GEE) to calculate adjusted odds ratios using an independence
model with geepack (v1.3.9). The following specified variables were
included in the models: age group, sex, size of transmission cluster,
phenotypic resistance to penicillin, phenotypic resistance to
tetracycline, phenotypic resistance to ciprofloxacin, phenotypic
resistance or decreased susceptibility to ceftriaxone and phenotypic
resistance to azithromycin. For sex, an ‘unknown’ category was included
to accommodate missing data, with all other categories having complete
data. For antimicrobial resistance phenotype, isolates were grouped
binarily as either resistant or susceptible/less susceptible/decreased
susceptibility except for ceftriaxone where there were no resistant
persistent isolates and isolates were grouped as either decreased
susceptibility or susceptible.

The following R code was used with
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Odds_ratio/odds_metadata.csv">Odds_ratio/odds_metadata.csv</a>
as input.

The input spreadsheet lists persistence in a binary system where 1 is
persistent and 0 is non persistent.

```         
library(tidyverse)
library(geepack)

Metadata = read.table("odds.csv", header = TRUE, sep= ",")

fit <- geeglm(formula = Persistence ~ Sex + AgeGroupSum + cluster_size + PEN + TET + CTRIX + CIPRO + AZITH, 
        data = Metadata, 
        id = cluster, 
        family = binomial, 
        corstr = "independence")

broom::tidy(x = fit, exp=T, conf.int = TRUE)
```

Results:

```         
  term                 estimate std.error statistic   p.value   conf.low conf.high
   <chr>                   <dbl>     <dbl>     <dbl>     <dbl>      <dbl>     <dbl>
 1 "(Intercept)"         0.00386   2.18        6.47  0.0110     0.0000534     0.279
 2 "SexFemale"           2.91      0.448       5.70  0.0170     1.21          7.00 
 3 "SexOther"            0.749     0.704       0.169 0.681      0.189         2.97 
 4 "SexUnknown"          0.111     1.86        1.39  0.239      0.00289       4.29 
 5 "AgeGroupSum20-29"    0.682     0.392       0.957 0.328      0.316         1.47 
 6 "AgeGroupSum30-39"    0.616     0.337       2.07  0.151      0.318         1.19 
 7 "AgeGroupSum40-49"    0.483     0.458       2.52  0.113      0.197         1.19 
 8 "AgeGroupSum≥50"      0.458     0.504       2.41  0.121      0.171         1.23 
 9 "cluster_size"        1.02      0.00801     4.82  0.0281     1.00          1.03 
10 "PENSUS"              2.71      1.22        0.669 0.413      0.249        29.5  
11 "TETSUS"              1.80      0.928       0.403 0.525      0.292        11.1  
12 "CTRIXSUS"            0.130     1.07        3.62  0.0570     0.0159        1.06 
13 "CIPROSUS"            2.62      1.07        0.811 0.368      0.323        21.2  
14 "AZITHSUS"          126.        1.16       17.3   0.0000326 12.9        1238.
```

## Whole genome alignments

### 1. Call variants

The trimmed paired end reads were aligned to the NCCP11945 reference
genome using Snippy (v4.3.5), requiring a minimum of ten supporting
reads and a variant frequency of 0.9 or greater:

```         
snippy --cpus 8 --minfrac 0.9 --mincov 10 --ref NCCP11945.fa --cleanup --outdir snippy/AUSMDU00008753 --prefix AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz
```

A pseudoalignment of all genomes was generated:

```         
snippy-core snippy/* --ref NCCP11945.fa
```

### 2. Whole genome ML phylogeny

Recombination filtering was performed using Gubbins (v2.4.1) with
default settings and the full Snippy psuedoalignments as input:

```         
run_gubbins.py --threads 10 core.full.aln
```

Following Gubbins, a SNP alignment was generated using snp-sites (v1)
and the Gubbins filtered alignment (_ploymorphic_sites.fasta) as input:

```         
snp-sites -c -o core.full.Gubbins.SNPs.aln core.full.Gubbins.aln
```

The number of constant sites was also calculated using snp-sites (v1):

```         
snp-sites -C core.full.aln
```

A ML phylogenetic tree was inferred using IQ-tree (v2.0.3), with the
best-fitting nucleotide substitution model chosen based on the lowest
BIC and the number of constant sites specified:

```         
iqtree -s core.full.Gubbins.SNPs.aln -B 1000 -T 60 -fconst 484258,580068,533403,495503
```

### 3. Subset alignments

Cluster whole genome pseudoalignments were subset from the entire
dataset pseudoalignments for each of the previously defined cgMLST
clusters with at least 5 isolates collected after 1st July 2019. A list
of these isolates and their clusters can be found in
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/clustersover5.txt" title="clustersover5.txt">Timed_trees/clustersover5.txt</a>.

Note: Isolates were renamed to append the decimal date to each ID for
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia?tab=readme-ov-file#6-bayesian-hierarchical-model">step
6</a>.

### 4. Recombination filtering

Recombination filtering was performed using Gubbins (v2.4.1) with
default settings for all clusters with the full whole genome
psuedoalignments as input:

```         
run_gubbins.py --threads 10 full_alignments/group_14_snippy.fasta
```

For each cluster gubbins filtered alignment (_ploymorphic_sites.fasta), a SNP alignment was
generated:

```         
snp-sites -c -o group_14_snippy_gubbins_SNPsites.fasta group_14_snippy_gubbins.fasta
```

### 5. Cluster ML and timed phylogenies

For each cluster whole genome pseudoalignment, the number of constant
sites was calculated using snp-sites (v1):

```         
snp-sites -C group_14_snippy.fasta
```

For each cluster gubbins filtered SNP alignment, ML phylogenetic trees
were inferred using IQ-tree (v2.0.3), with the best-fitting nucleotide
substitution model chosen based on the lowest BIC and the number of
constant sites specified:

```         
iqtree -s group_14_snippy_gubbins_SNPsites.fasta -B 1000 -T 60 -fconst 484048,535722,549723,494991
```

All ML phylogenies can be found in
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/ML_trees.zip">Timed_trees/ML_trees.zip</a>.

Molecular dating of ancestral events was performed on the resulting ML
trees, using the least-squares dating (LSD) software (v0.3) with a rate
of 4.5x10-6 substitutions per site as previously defined:

```         
/home/taouk/lsd-0.3beta-master/src/lsd -d dates.txt -i group_14_ML.tree -c -r a -w rate.txt
```

The dates.txt input is a tab separated file that has the IDs in one
column and the date of collection in decimal format in another for each
isolate. A separate file was generated for each cluster. The decimal
date for each isolate can be found in
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/dates.txt" title="dates.txt">Timed_trees/dates.txt</a>.

All LSD timed phylogenies can be found in
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/LSD_trees.zip">Timed_trees/LSD_trees.zip</a>.

### 6. Bayesian hierarchical model

The timed phylogenetic trees for each cluster were input into a Bayesian hierarchical model. We applied a birth-death skyline process to all trees. Under this model, the epidemiological process follows a birth-death process, where branching events in a phylogenetic tree are informative about transmission, and the sampling process is directly modelled. The epidemiological parameters can change in a piecewise fashion over two intervals, with the interval time estimated as part of the model. Each transmission cluster had a fixed average infection duration of three months and shared a sampling proportion parameter, meaning that sampling intensity is a free parameter that is estimated here. However, each cluster was allowed to have independent effective reproductive numbers (Re) and epidemic origin times. The Re values were permitted to vary at a specific point in time, referred to as the "time slice," which represented a significant point of change in Re values. Although Re values were independent for each cluster, they were governed by a single Gamma prior distribution with two hyperparameters: shape and mean. The shape of this distribution influenced the skew, with smaller values indicating more variation in Re among clusters, and larger values suggesting similarity in Re values across clusters.

In our model, the Re values for clusters followed a Gamma distribution, where the shape parameter reflected heterogeneity in the spread among clusters. We maintained a fixed infection duration of three months in our hierarchical model, reflecting what we considered was the most plausible scenario. The birth-death model design resulted in different magnitudes of Re values, although the trends remained consistent.  We assumed a uniform prior distribution for all origin times, ranging from zero to six months. Before the first genome of each cluster was sampled, we set the sampling proportion to 0.0, and thereafter, it was modelled with a Beta(1, 30) prior distribution to capture our assumption that the sampling proportion was at most 10%. The hyperparameters of the Gamma distribution for Re values were assigned Gamma(10, 1) priors. Additionally, the time slice allowing Re value changes was assigned a uniform prior distribution between March 20 and March 30, 2020.

We sampled the posterior distribution using Markov chain Monte Carlo, implemented in BEAST2.6, with a chain length of 109 steps and sampling every 105 steps. As the phylogenetic trees were fixed, there were no calculations of phylogenetic likelihood, making this analysis computationally more efficient compared to those involving both phylogenetic and phylodynamic likelihoods.

The XML: 
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/hierarchical_with_feast_D90days.xml">Timed_trees/hierarchical_with_feast_D90days.xml</a>

# Supplementary Analyses

In addition to methods used in the main manuscript and analyses,
supplementary analyses are included in this GitHub. All supplementary
analyses can be found in <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/tree/main/Supplementary_analyses">Supplementary_analyses</a>. This folder includes two analyses:

-   <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/tree/main/Supplementary_analyses/cgMLST_method">Supplementary_analyses/cgMLST_method</a>:
    Using cgMLST to define transmission groups in another publicly
    available international *N. gonorrheoae* dataset.
-   <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/tree/main/Supplementary_analyses/SNP_alignment">Supplementary_analyses/SNP_alignment</a>:
    Comparing the use of a traditional strict core alignment for
    phylogeny and clustering compared to cgMLST.
