# Longitudinal Genomic Analysis of *Neisseria gonorrhoeae* Transmission Dynamics in Australia

Mona L. Taouk, George Taiaroa, Sebastian Duchene, Soo Jen Low, Charlie K. Higgs, Darren Y. J. Lee, Shivani Pasricha, Nasra Higgins, Danielle J. Ingle, Benjamin P. Howden, Marcus Y. Chen, Christopher K. Fairley, Eric P. F. Chow, Deborah A. Williamson



## QC reads 

### 1. Kraken
`for i in $(cat IDs.txt); do kraken2 --db /home/linuxbrew/db/kraken2/pluspf --gzip-compressed --paired /home/taouk/NGtransmission/reads/${i}_1.fq.gz /home/taouk/NGtransmission/reads/${i}_2.fq.gz --report ${i}.kraken.txt; done`

`grep -m1 -P "\tS\t" *.txt | sort -k2nr | sed 's/.txt: /\t/' > kraken_results.txt`

### 2. NG reads
`minimap2 -t10 -ax sr NCCP11945.fa AUSMDU00008753_1.fq.gz AUSMDU00008753_2.fq.gz | samtools fastq -f 2 - | grep -c "^@" > AUSMDU00008753_count.txt`

`grep "" *_count.txt | sed 's/_count.txt:/\t/' > test.txt`

### 3. Average read length
`fq --ref NCCP11945.fa AUSMDU00008753_1.fq.gz AUSMDU00008753_2.fq.gz > AUSMDU00008753.yield.tab`

`grep -m1 -P "AvgLen" *.yield.tab | sort -k2nr > yield.all3.txt`

## Assemblies

### 1. Shovill
`shovill --outdir assemblies/AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz`

`for i in $(cat IDs.txt); do mv ${i}/contigs.fa ${i}/${i}.contigs.fa; done`

### 2. Number of contigs
`grep -c "^>" */*contigs.fa > contigs_all.txt`

## Genotyping

### 1. MLST
`mlst --legacy --scheme neisseria assemblies/*/*.contigs.fa > mlst_results.txt`

### 2. NG-MAST
`ngmaster assemblies/*/*.contigs.fa > ngmast_results.txt` 

### 3. NG-STAR
`python3 /home/taouk/pyngSTar/pyngSTar.py -f -a -i assemblies/*/*.contigs.fa -p /home/taouk/pyngSTar/pyngSTarDB/ > ngstar_results.txt`

## cgMLST

### 1. Prepare schema
Make training file:

`prodigal -i NCCP11945.fa -t neisseria_gonorrhoeae.trn -p single`

Prepare the downloaded schema:

`chewBBACA.py PrepExternalSchema -i scheme_directory -o scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn --cpu 50 --st 0.3`

The training file I used can be found <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/neisseria_gonorrhoeae.trn" title="Here">here</a>.

The scheme can be found <a href="https://doi: 10.1093/infdis/jiaa002" title="Harrison et al.">Harrison et al.</a> and can be downloaded from <a href="https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=schemeInfo&scheme_id=62" title="PubMLST">PubMLST</a>.

### 2. Allele calling
`chewBBACA.py AlleleCall -i paths_to_assemblies.txt -g scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn -o results --cpu 50`

### 3. Refining schema 
`chewBBACA.py ExtractCgMLST -i results/results_20211011T232621/results_alleles.tsv --r results/results_20220602T014232/RepeatedLoci.txt -o evaluate95 --t 0.95` 

cgMLST allele calling final results:

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLST.tsv">cgMLST.tsv</a>: Table of alleles for each gene in the scheme for all isolates.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/cgMLSTschema.txt">cgMLSTschema.txt</a>: 1495 genes included in final schema for this study.

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/mdata_stats.tsv">mdata_stats.tsv</a>: Number of and percentage of uncalled genes per genome in the dataset. 

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Results/presence_absence.csv">presence_absence.csv</a>: Gene presence/absence table. (1 present, 0 absent).


### 4. Transform allele calling spreadsheet to distance matrix
`cgmlst-dists results/evaluate/cgMLST.tsv > pairwise_distances_0.95_cgMLST.txt` 

### 5. cgMLST alginment (only including 5,881 genomes)
<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/cgMLST_alignment_part1.r" title="R Script to preprocess the data">R Script to preprocess the data</a>

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/cgMLST_alignment_part2.r" title="Bash script to make the alignments">Bash script to make the alignments</a>

### 6. 95% core alignment (only including 5,881 genomes)
`trimal -in unique_pass.aln -out stripped_alignment.fasta -gt 0.05 -threads 40`

### 10. Maximum likelihood phylogeny
`iqtree -s stripped_alignment.fasta -B 1000 -T 40`

### 11. Timed phylogeny
`/home/taouk/lsd-0.3beta-master/src/lsd -d dates.tsv -i 95gapsMP.tree -c`

### 12. Determine threshold and adjust for timed between samples

R code on methods to determine threshold, adjust thresholds and call clusters: <a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/cgMLST/Threshold.Rmd" title="Threshold.Rmd">Threshold.Rmd</a>.

Data can be downloaded from <a href="https://doi.org/10.26188/25989001" title="here">here</a>.

## Snippy

### 1. Run Snippy on whole dataset
`snippy --cpus 8 --minfrac 0.9 --mincov 10 --ref NCCP11945.fa --cleanup --outdir snippy/AUSMDU00008753 --prefix AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz`

`snippy-core snippy/* --ref NCCP11945.fa`

### 2. Gubbins on 5881 unique genomes
`run_gubbins.py 5881.full.aln`

### 3. Whole dataset phylogeny
`coresnpfilter -e -c 1.0 gubbins.output > 5881.core.aln`

`iqtree -s 5881.core.aln -B 1000 -T 40`

### 4. Snp distance matrix including callibration isolates
`run_gubbins.py 6215.full.aln`

`snp-dists gubbins.output.FULL > 6215_snp_distance_matrix`

## BEAST

### 1. Subset alignments
Make alignments for each cluster:

`bash split.sh`

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/split.sh" title="split.sh">split.sh</a>

### 2. Make Maximum likelihood phylogenies and timed phylogenetic trees from alignments
Make iqtree trees for each cluster alignments and then input into LSD to make timed trees:

`bash trees.sh`

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/trees.sh" title="trees.sh">trees.sh</a>

<a href="https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Timed_trees/dates.sh" title="here">dates.sh</a>

### 3. Run BEAST 
Script from Sebastian to be added here

