# NG-transmission-Australia
Code used in the study: Longitudinal Genomic Analysis of *Neisseria gonorrhoeae* Transmission Dynamics in Australia

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

### 1. Prepare Schema
`chewBBACA.py PrepExternalSchema -i scheme_directory -o scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn --cpu 50 --st 0.3`

### 2. Allele Calling
`chewBBACA.py AlleleCall -i paths_to_assemblies.txt -g scheme_prepped --ptf /home/taouk/NGtransmission/cgMLST/neisseria_gonorrhoeae.trn -o results --cpu 50`

### 3. Refining Schema 
`chewBBACA.py ExtractCgMLST -i results/results_20211011T232621/results_alleles.tsv --r results/results_20220602T014232/RepeatedLoci.txt -o evaluate95 --t 0.95` 

### 4. Transform allele calling spreadsheet to distance matrix
`cgmlst-dists results/evaluate/cgMLST.tsv > pairwise_distances_0.95_cgMLST.txt` 

### 5. cgMLST alginment (only including 5,881 genomes)
R Script to preprocess the data

`R cgMLST_alignment_part1.r`

Bash script to make the alignments

`bash cgMLST_alignmen_part2.sh`

### 6. 95% core alignment (only including 5,881 genomes)
`trimal -in unique_pass.aln -out stripped_alignment.fasta -gt 0.05 -threads 40`

### 10. Tree
`iqtree -s stripped_alignment.fasta -B 1000 -T 40`

### 11. Timed tree
`/home/taouk/lsd-0.3beta-master/src/lsd -d dates.tsv -i 95gapsMP.tree -c`

### 12. Determine threshold and adjust for timed between samples
see cgMLST_threshold/Threshold.Rmd
Data can be downloaded from https://doi.org/10.26188/25989001 

## Snippy

### 1. Run Snippy on whole dataset
`snippy --cpus 8 --minfrac 0.9 --mincov 10 --ref NCCP11945.fa --cleanup --outdir snippy/AUSMDU00008753 --prefix AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz`

`snippy-core snippy/* --ref NCCP11945.fa`

### 2. Gubbins on 5881 unique genomes
`run_gubbins.py 5881.full.aln`

### 3. Whole dataset tree
`coresnpfilter -e -c 1.0 gubbins.output > 5881.core.aln`

`iqtree -s 5881.core.aln -B 1000 -T 40`

### 4. Snp distance matrix including callibration isolates
`run_gubbins.py 6215.full.aln`

`snp-dists gubbins.output.FULL > 6215_snp_distance_matrix`

## BEAST

### 1. subset alignments from core
`bash split.sh`

### 2. Make iqtree trees from alignemnts
`bash trees.sh`

### 3. Run BEAST 
script from Sebastain

