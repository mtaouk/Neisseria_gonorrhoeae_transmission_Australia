#!/bin/bash
for i in `cat clustersover5.txt`; do mkdir "$i" && cp /home/taouk/NGtransmission/Snippy_subset_2/alignments/group_"$i"_snippy.fasta "$i"/group_"$i"_snippy.fasta; done
for i in `cat clustersover5.txt`; do iqtree -s "$i"/group_"$i"_snippy.fasta -B 1000 -T 60; done
for i in `cat clustersover5.txt`; do cp "$i"/group_"$i"_snippy.fasta.treefile group_"$i".tree; done

###LSD 

bash dates.sh

for i in `cat clustersover5.txt`; do cp /home/taouk/NGtransmission/Snippy_subset_2/full_trees/group_"$i".tree "$i"/group_"$i".tree; done
for i in `cat clustersover5.txt`; do cd "$i" && /home/taouk/lsd-0.3beta-master/src/lsd -d dates_"$i" -i group_"$i".tree -c -r a -w /home/taouk/NGtransmission/LSD/rate.txt && cd ../; done
mkdir timed_trees
for i in `cat clustersover5.txt`; do cp "$i"/group_"$i".tree.result.date.nexus timed_trees/group_"$i"_timed.tree; donecluster_trees