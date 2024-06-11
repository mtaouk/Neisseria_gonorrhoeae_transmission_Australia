# Longitudinal Genomic Analysis of *Neisseria gonorrhoeae* Transmission Dynamics in Australia

Mona L. Taouk, George Taiaroa, Sebastian Duchene, Soo Jen Low, Charlie K. Higgs, Darren Y. J. Lee, Shivani Pasricha, Nasra Higgins, Danielle J. Ingle, Benjamin P. Howden, Marcus Y. Chen, Christopher K. Fairley, Eric P. F. Chow, Deborah A. Williamson


## METHODS

## QC reads 

### 1. Kraken
Trimmed paired end reads for each genome were input into Kraken2 (v2.1.1) to investigate for species contamination:

`for i in $(cat IDs.txt); do kraken2 --db /home/linuxbrew/db/kraken2/pluspf --gzip-compressed --paired /home/taouk/NGtransmission/reads/${i}_1.fq.gz /home/taouk/NGtransmission/reads/${i}_2.fq.gz --report ${i}.kraken.txt; done`

The top species match and percentage were extracted from each result file and summarised:

`grep -m1 -P "\tS\t" *.txt | sort -k2nr | sed 's/.txt: /\t/' > kraken_results.txt`

### 2. Sequencing depth
Trimmed paired end reads for each genome were aligned to the NCCP11945 reference genome (GenBank accession NC_011035.1) using Minimap2 (v2.17-r941) and Samtools (v1.10, using HTSlib v1.10.2) to calculate the number of reads aligning to the reference genome:

`minimap2 -t10 -ax sr NCCP11945.fa AUSMDU00008753_1.fq.gz AUSMDU00008753_2.fq.gz | samtools fastq -f 2 - | grep -c "^@" > AUSMDU00008753_count.txt`

The outputs were combined:

`grep "" *_count.txt | sed 's/_count.txt:/\t/' > test.txt`

The average read length for each genome was calculated using fq (v0.11.0) and extracted from the results:

`fq --ref NCCP11945.fa AUSMDU00008753_1.fq.gz AUSMDU00008753_2.fq.gz > AUSMDU00008753.yield.tab`

`grep -m1 -P "AvgLen" *.yield.tab | sort -k2nr > yield.all3.txt`

The average read depth was calulated as follows: (number of reads mapping to reference * average read length)/reference genome length

## Assemblies

### 1. Shovill

De novo genome assemblies were generated using Shovill (v1.1.0):

`shovill --outdir assemblies/AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz`

The names of each genome's contigs were changed from the default:

`for i in $(cat IDs.txt); do mv ${i}/contigs.fa ${i}/${i}.contigs.fa; done`

### 2. Number of contigs
The number of contigs for each genome was calculated:

`grep -c "^>" */*contigs.fa > contigs_all.txt`

## Genotyping

### 1. MLST

MLSTs were identified using mlst (v2.23.0; PubMLST database accessed May 5th 2023):

`mlst --legacy --scheme neisseria assemblies/*/*.contigs.fa > mlst_results.txt`

### 2. NG-MAST

NG‐MASTs were assigned using NGmaster (1.0.0; NG-MAST v2.0 PubMLST database accessed May 5th, 2023):

`ngmaster assemblies/*/*.contigs.fa > ngmast_results.txt` 

### 3. NG-STAR

The NG-STAR typing scheme was used to determine the phenotypic profile of seven resistance genes (*penA*, *mtrR*, *porB*, *ponA*, *gyrA*, *parC* and *23S* rRNA), using pyngSTar(database version 2.0). NG-STAR alleles and profiles obtained from the database hosted by the <a href="https://ngstar.canada.ca/alleles/loci_selection?lang=en" title="Public Health Agency of Canada">Public Health Agency of Canada</a>:

`python3 /home/taouk/pyngSTar/pyngSTar.py -f -a -i assemblies/*/*.contigs.fa -p /home/taouk/pyngSTar/pyngSTarDB/ > ngstar_results.txt`

## cgMLST

### 1. Prepare schema
A prodigal training file for *N. gonorrhoeae* was made using the NCCP11945 reference genome using prodigal (v.2.6.3):

`prodigal -i NCCP11945.fa -t neisseria_gonorrhoeae.trn -p single`

Training file: <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/neisseria_gonorrhoeae.trn" title="neisseria_gonorrhoeae.trn">neisseria_gonorrhoeae.trn</a>.

The *N. gonorrhoeae* cgMLST schema developed by <a href="https://doi: 10.1093/infdis/jiaa002" title="Harrison et al.">Harrison et al.</a> was downloaded from <a href="https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=schemeInfo&scheme_id=62" title="PubMLST">PubMLST</a> and prepared for use using chewBBACA (v2.8.5):

`chewBBACA.py PrepExternalSchema -i scheme_directory -o scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn --cpu 50 --st 0.3`

### 2. Allele calling
Allele calling was performed on all *N. gonorrhoeae* genomes, using the assemblies as input:

`chewBBACA.py AlleleCall -i paths_to_assemblies.txt -g scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn -o results --cpu 50`

### 3. Refining schema 
The schema was refined to only include genes present in 95% of genomes in the dataset:

`chewBBACA.py ExtractCgMLST -i results/results_20211011T232621/results_alleles.tsv --r results/results_20220602T014232/RepeatedLoci.txt -o evaluate95 --t 0.95` 

cgMLST allele calling final results:

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLST.tsv">cgMLST.tsv</a>: Table of alleles for each gene in the scheme for all isolates.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLSTschema.txt">cgMLSTschema.txt</a>: 1495 genes included in final schema for this study.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/mdata_stats.tsv">mdata_stats.tsv</a>: Number of and percentage of uncalled genes per genome in the dataset. 

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/presence_absence.csv">presence_absence.csv</a>: Gene presence/absence table. (1 present, 0 absent).

### 4. cgMLST alleleic differences
The cgMLST results table was transformed to a symmetrical distance matrix using cgmlst-dists:

`cgmlst-dists results/evaluate/cgMLST.tsv > cgMLST_matrix.txt` 

Resulting matrix: <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLST_matrix.csv.zip">cgMLST_matrix.csv.zip</a>.

The values from this matrix were adjusted in step 12.  


### 5. cgMLST alginment (only including 5,881 genomes)
R Script to preprocess the data: <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/cgMLST_alignment_part1.r">cgMLST_alignment_part1.r</a>
This script takes the cgMLST.tsv result file and converts it into a long transposed version that lists the allele called at each gene for each isolate.

Shell script to make the alignments: <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/cgMLST_alignment_part2.r">cgMLST_alignment_part2.r</a>
This script takes the transposed output file from the above R script and uses this as input. This script also requires a second input (all_alleles_cat.fasta), which is a fasta file that contains all the alleles for all the genes in the cgMLST database. This can be generated by concatenating the fasta files contiaing all alleles for each gene together. This file will be ~1.7GB.

### 6. 95% core alignment (only including 5,881 genomes)
The alignment of 5,881 genomes representing a unique infection each, was processed to create a 95% soft core using trimal (v1.4.rev15):

`trimal -in unique_pass.aln -out stripped_alignment.fasta -gt 0.05 -threads 40`

### 10. Maximum likelihood phylogeny
A maximum liklihood phylogenetic tree was generated using IQ-tree (v2.0.3):

`iqtree -s stripped_alignment.fasta -B 1000 -T 40`

The best-fitting nucleotide substitution model based on the lowest Bayesian Information Criterion (BIC) was GTR+I+R10.

### 11. Timed phylogeny
Molecular dating of ancestral events was performed using the least-squares dating (LSD) software (v0.3), with the maximum liklihood phylogeny generated in the above step used as input:

`/home/taouk/lsd-0.3beta-master/src/lsd -d dates.tsv -i 95gapsMP.tree -c`

### 12. Determine threshold for clustering
The following R code details steps to determine a threshold for clustering, ajusting the threshold for time between sample collection dates and performing clustering: <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Threshold.Rmd" title="Threshold.Rmd">Threshold.Rmd</a>.

The data required for the code can be downloaded from <a href="https://doi.org/10.26188/25989001">doi.org/10.26188/25989001</a>.

## Odds ratios

Add in the R code and input files when finalised.

## Whole genome alignments

### 1. Call variants

The trimmed paried end reads were aligned to the NCCP11945 reference genome using Snippy (v4.3.5), requiring a minimum of ten supporting reads and a variant frequency of 0.9 or greater:

`snippy --cpus 8 --minfrac 0.9 --mincov 10 --ref NCCP11945.fa --cleanup --outdir snippy/AUSMDU00008753 --prefix AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz`

A pseudoalignment of all genomes was generated:

`snippy-core snippy/* --ref NCCP11945.fa`

### 2. Subset alignments

Cluster whole genome pseudoalignments were subset from the entire dataset pseudoalignments for each of the previously defined cgMLST clusters with at least 5 isolates collected after 1st July 2019. A list of these isolates and their clusters can be found in <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/clustersover5.txt" title="clustersover5.txt">clustersover5.txt</a>. Isolates have been renamed to append the decimal date to each ID for Step 4. 

### 3. ML and timed phylogenies

For each cluster alignment, ML phylogenetic trees were inferred using IQ-tree (v2.0.3), with the best-fitting nucleotide substitution model chosen based on the lowest BIC. Molecular dating of ancestral events was performed on the remaining ML trees, using the least-squares dating (LSD) software (0.3) with -r -c parameters and a rate of 4.5x10-6 substitutions per site as previously defined. The subsequent timed trees were used as input for a Bayesian hierarchical model. 

In short, ML phylogenetic trees were made for each subset alignments:

`iqtree -s snippy.fasta -B 1000 -T 60`

and timed phylogenetic trees were made from each ML phylogeny:

`/home/taouk/lsd-0.3beta-master/src/lsd -d dates.txt -i cluster.tree -c -r a -w rate.txt` 

The rate in the rate file was 4.54E-6.
The dates.txt is a tab seperated file that has the IDs in one column and the date of collection in decimal format for each isolate. A sepearte file was generated for each cluster. The decimal date for each isolate can be found in <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/dates.txt" title="dates.txt">dates.txt</a>.

### 4. Bayesian hierarchical model 

The timed phylogenetic trees for each cluster were input into a Bayesian hierarchical model using this script. TO ADD.

