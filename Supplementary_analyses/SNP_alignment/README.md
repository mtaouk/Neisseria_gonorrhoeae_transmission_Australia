# Comparing cgMLST to SNP based approaches

The trimmed paired end reads were aligned to the NCCP11945 reference
genome using Snippy (v4.3.5), requiring a minimum of ten supporting
reads and a variant frequency of 0.9 or greater. Recombination filtering
was performed using Gubbins (v2.4.1) with default settings and the full
Snippy pseudoalignments as input. Following Gubbins, a core SNP
alignment was generated using snp-sites (v1) and the Gubbins filtered
alignment as input with the -c flag. The number of constant sites from
the whole genome pseudoalignment was also calculated using snp-sites
with the -C flag (v1). A ML phylogenetic tree was inferred using IQ-tree
(v2.0.3), with the best-fitting nucleotide substitution model chosen
based on the lowest BIC and the number of constant sites specified.
Molecular dating of ancestral events was performed using the
least-squares dating (LSD) software (v0.3), with the whole dataset
maximum likelihood phylogeny generated here used an input. All code is
here:

```         
snippy --cpus 8 --minfrac 0.9 --mincov 10 --ref NCCP11945.fa --cleanup --outdir snippy/AUSMDU00008753 --prefix AUSMDU00008753 --R1 /home/taouk/NGtransmission/reads/AUSMDU00008753_1.fq.gz --R2 /home/taouk/NGtransmission/reads/AUSMDU00008753_2.fq.gz

snippy-core snippy/* --ref NCCP11945.fa

run_gubbins.py --threads 10 core.full.aln

snp-sites -c -o core.full.Gubbins.SNPs.aln core.full.Gubbins.aln

snp-sites -C core.full.aln

iqtree -s core.full.Gubbins.SNPs.aln -B 1000 -T 60 -fconst 484258,580068,533403,495503

/home/taouk/lsd-0.3beta-master/src/lsd -d dates.tsv -i core_SNP.tree -c -r a
```

Using this method, the core SNP alignment consisted of 8,842 core sites
(polymorphic/variants that are present in all samples). While using a
core SNP alignment can be used to build a high-resolution phylogeny in
many cases and has been the gold standard approach in bacterial
phylogenetics for the past decade, when applied to a large and diverse
dataset as in this case, it results in a shrinking of informative sites,
and a less resolute phylogeny. N. gonorrhoeae is a very diverse species
with much recombination, therefor the number of sites conserved across
various samples is smaller. For example, across the full whole genome
pseudoalignment, there is a minimum of 6% N sites (134,866 bp) for any
genome. As these N sites will be dispersed mostly randomly across the
genome, the chances of any site having at least one N in at least one
sample is high and means that site will be excluded from the core SNP
alignment, even if it is informative to the phylogeny. One way to
increase the number of core SNP sites is to remove genomes with a high
proportion of N sites from the alignment and analysis. In this case,
1,184 genomes have more than 10% N sites. Excluding these would remove
20% of isolates from the phylogeny and clustering analysis – decreasing
our sampling proportion and introducing a high level of uncertainty into
our clustering. As a result, we opted to use a cgMLST method, a much
more permissive way of comparing relatedness across very diverse genomes
in a large dataset.

The concern that using a strict core results in isolates potentially
being classified as more closely related than they would be, stems from
the principle of using a static SNP threshold to define transmission.
For example, a common SNP threshold of 10 SNPs would result in much more
permissive clustering using a core SNP alignment of 8,842 compared to
using the same threshold on a larger core, as the threshold represents a
fraction of the total sites.

Regardless, we have generated a recombination filtered core genome SNP
ML phylogeny. We can see that the overall population structure is mostly
conserved, however the resolution at recent evolutionary events is much
reduced when compared to the concatenated cgMLST phylogeny. We also
generated a timed phylogeny from the recombination masked core genome
SNP ML phylogeny which produced similar results to the timed tree
generated from the concatenated cgMLST alignment. Both phylogenies are
rooted at the midpoint.

## SNP tree

![](images/SNP_tree.pdf){width="188.8cm" height="20cm"}
