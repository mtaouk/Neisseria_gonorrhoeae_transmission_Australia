#!/bin/bash


#1. Rename headers in alleles fasta file to simplified format (eg. NEIS0001_1)
## Replace first and last underscore with tab in headers, then merge first and last field to get new clean headers
grep ">" all_alleles_cat.fasta | sed 's/>//g' > oldheaders
cat oldheaders | sed 's/_/\t/1' | sed 's/\(.*\)_/\1\t/' | awk -F "\t" '{ print $1, $NF }' OFS="_" | sed 's/>//g' > cleanheaders
# Make alias file for new headers
paste oldheaders cleanheaders > alias.txt
rm oldheaders cleanheaders
# Rename sequences using seqkit
seqkit replace -p '^(\S+)' -r '{kv}' -k alias.txt all_alleles_cat.fasta > all_alleles_cat.fasta.renamed
rm alias.txt


#2. Convert renamed alleles fasta to tab-delimited format
seqkit fx2tab -j 4 -Q all_alleles_cat.fasta.renamed > all_alleles_cat.fasta.renamed.tab


#3. Split melted cgMLST data frame by geneid (into 1495 files)
# Create directory to store split files
mkdir lookup
awk -F "\t" '{ print >"lookup/"$2".txt" }' cgMLST_transposed.tsv


#4. Split alleles tab-delimited file by geneid (into 1495 files). We will use this as a lookup file to get sequences for each genome
awk '{ print $1, $1, $2 }' OFS="\t" all_alleles_cat.fasta.renamed.tab | sed 's/_/\t/2' | cut -f1,2,4 > all_alleles_cat.fasta.renamed.tab.tosplit
awk -F "\t" '{ print >"lookup/"$2".lookup" }' all_alleles_cat.fasta.renamed.tab.tosplit


#5. Create list of genes, then loop through each gene, join genomes to sequences
ls lookup/*.txt | awk -F "/" '{ print $2 }' | sed 's/.txt//g' > genes
mkdir individual_genes
for i in `cat genes`; do sort -k1,1g -o lookup/"$i".lookup lookup/"$i".lookup && sort -k3,3g -o lookup/"$i".txt lookup/"$i".txt; done
for i in `cat genes`; do csvtk join --left-join --na "-" -H -T -t -f "3;1" lookup/"$i".txt lookup/"$i".lookup | cut -f1,5 | sed 's/^/>/g' | sed 's/\t/\n/g' > individual_genes/"$i".fasta; done
## check length
seqkit stats individual_genes/*.fasta > length_check.txt


#6. Align the individual gene fastas by submitting a job array on spartan or equivalent- I noticed quite a few genes have uniform lengths throughout, so you may or may not want to 'align' those
cd individual_genes
for i in $(cat ../genes); do echo "mafft --thread 40 ${i}.fasta > ${i}.aln" >> input_aligning.txt; done
dos2unix input_aligning.txt
parallel -j 20 -k --bar {} :::: input_aligning.txt


#7. Concatenate all individual gene alignments
cd ../
seqkit concat individual_genes/NEIS*.aln --full > concat_1495.aln
dos2unix exclude.txt
seqkit grep -v -f exclude.txt concat_1495.aln -o unique_pass.aln
python3 filter_fasta_by_list_of_headers.py --id_file exclude.txt --input concat_1495.aln --output unique_pass.aln
seqtk comp unique_pass.aln > filter_alignment_stats.txt