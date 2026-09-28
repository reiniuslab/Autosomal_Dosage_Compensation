#!/bin/bash

# summarise mappability scores and nucleotide frequencies per bin
mkdir ucsc

wget -P ucsc https://bismap.hoffmanlab.org/raw/mm10/k50.umap.bedgraph.gz
zcat ucsc/k50.umap.bedgraph.gz | sed -e "1d; s/^chr//" | sort -k1,1 -k2,2n | gzip > ucsc/k50.umap.nochr.bedgraph.gz
bedtools map -a cnv/clonal_DNAseq/genomebins_100kb.bed -b ucsc/k50.umap.nochr.bedgraph.gz -c 4 -o mean | gzip > ucsc/mappability_k50_100kb.bedgraph.gz

gunzip -k GRCm38_68_CAST_EiJ_N-masked.fa.gz
bedtools nuc -fi GRCm38_68_CAST_EiJ_N-masked.fa -bed cnv/clonal_DNAseq/genomebins_100kb.bed | gzip > ucsc/nucfreq_100kb.tsv.gz
rm GRCm38_68_CAST_EiJ_N-masked.fa

wget -P ucsc https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz # assembly gaps
zcat ucsc/gap.txt.gz | cut -f2- | awk '$8 == "no" {OFS="\t"; gsub("^chr","",$1);print}' | sort -k1,1 -k2,2n | gzip > ucsc/gap.nochr.bed.gz
bedtools window -w 500000 -c -a cnv/clonal_DNAseq/genomebins_100kb.bed -b ucsc/gap.nochr.bed.gz | gzip > ucsc/gap_100kb.bedgraph.gz

wget -P ucsc https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
zcat ucsc/rmsk.txt.gz | cut -f6,7,8,11,12,13 | sed "s/^chr//" | sort -k1,1 -k2,2n | gzip > ucsc/rmsk.nochr.bed.gz
bedtools intersect -wao -a cnv/clonal_DNAseq/genomebins_100kb.bed -b ucsc/rmsk.nochr.bed.gz | bedtools map -c 10 -a cnv/clonal_DNAseq/genomebins_100kb.bed -b - | awk '{print $1,$2,$3,$4/($3-$2)}' | gzip > ucsc/rmsk_100kb.bedgraph.gz