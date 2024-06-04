#!/bin/bash

mkdir inputs

awk -F "\t" '{print >"inputs/"$2".txt"}' clusters.txt

ls inputs/*.txt | awk -F "/" '{ print $2 }' | sed 's/.txt//g' > IDs.txt

for i in `cat IDs.txt`; do awk -F "\t" '{print $1}' inputs/"$i".txt > inputs/"$i".filtered.txt; done

mkdir snippy

for i in `cat IDs.txt`; do seqkit grep -f inputs/"$i".filtered.txt snippy_decimalnames_unique_pass.aln -o snippy/group_"$i"_snippy.fasta; done
