#!/bin/bash

##

trim_data(){
	# $1 in1
	# $2 in2
	# $3 folder clean
	# $4 basename

	mkdir ${3}/fastp_reports/${4}
	
	fastp --in1 $1 --in2 $2 --out1 ${3}/${1##*/} --out2 ${3}/${2##*/} --json ${3}/fastp_reports/${4}/fastp.json --html ${3}/fastp_reports/${4}/fastp.html -R ${4}
}
export -f trim_data

align_data(){
	# $1 index
	# $2 Read1
	# $3 Read2
	# $4 output
	minimap2 -ax sr $1 $2 $3 | samtools sort -O bam -o ${4}
	samtools index ${4}
}
export -f align_data

align_stats(){
	# $1 input
	# $2 output
	samtools stats $1 | grep "^SN" | cut -f 2- > $2
}
export -f align_stats

downsample_bam(){
	# $1 in
	# $2 out

	frac=$( samtools idxstats $1 | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=100000/total; if (frac > 1) {print 1} else {print frac}}' )
	samtools view -bs $frac $1 > $2
	samtools index $2
}
export -f downsample_bam

##

bcl2fastq --sample-sheet SampleSheet_clonal_DNAseq.csv --no-lane-splitting -R /huxley.hdd/raw/230802_NB552166_0110_AHJJG3BGXT/ -o ./fastq/clonal_DNAseq

## Quality trim data

mkdir clean clean/clonal_DNAseq clean/clonal_DNAseq/fastp_reports

parallel -j 12 --link --rpl '{%} s:.*/::; s:_S\d{1,3}.*::;' \
trim_data {1} {2} clean/clonal_DNAseq {1%} \
::: fastq/clonal_DNAseq/*_R1_001.fastq.gz ::: fastq/clonal_DNAseq/*_R2_001.fastq.gz

multiqc -d -n qc_report_clonal_DNAseq clean/clonal_DNAseq/fastp_reports

## Align data

mkdir bam bam/clonal_DNAseq

# cat /huxley.hdd/reference_genomes/mouse/GRCm38.mgp/CAST_EiJ_N-masked/chr[123456789XYM]*.fa | bgzip > GRCm38_68_CAST_EiJ_N-masked.fa.gz

parallel -j 12 --link --rpl '{%} s:.*/::; s:_S\d{1,3}.*::;' \
align_data GRCm38_68_CAST_EiJ_N-masked.fa.gz {1} {2} bam/clonal_DNAseq/{1%}.bam \
::: clean/clonal_DNAseq/*_R1_001.fastq.gz ::: clean/clonal_DNAseq/*_R2_001.fastq.gz

## Alignment statistics

mkdir alignment_stats alignment_stats/clonal_DNAseq

parallel -j 12 samtools flagstat {} ">" alignment_stats/clonal_DNAseq/flagstat_{/.}.txt ::: bam/clonal_DNAseq/*.bam
parallel -j 12 align_stats {} alignment_stats/clonal_DNAseq/stats_{/.}.txt ::: bam/clonal_DNAseq/*.bam
parallel -j 12 samtools coverage {} ">" alignment_stats/clonal_DNAseq/coverage_{/.}.txt ::: bam/clonal_DNAseq/*.bam

## CNVs

mkdir out_DNA out_DNA/cnv out_DNA/cnv/clonal_DNAseq

# generate binned read counts 
samtools faidx GRCm38_68_CAST_EiJ_N-masked.fa.gz
awk '!/^A|^K/ {OFS="\t"; print $1,$2}' GRCm38_68_CAST_EiJ_N-masked.fa.gz.fai > GRCm38_68.chromsizes.txt

bedtools makewindows -w 100000 -g GRCm38_68.chromsizes.txt | sort -k1,1 -k2,2n > cnv/clonal_DNAseq/genomebins_100kb.bed
((echo -e 'CHROM START END' bam/clonal_DNAseq/!(Undetermined).bam | tr ' ' '\t' ) && bedtools multicov -p -q 13 -bams bam/clonal_DNAseq/!(Undetermined).bam -bed cnv/clonal_DNAseq/genomebins_100kb.bed) | gzip > cnv/clonal_DNAseq/read.counts.100kb.bin.tsv.gz

# generate BAF counts (NOTE: slow)
bcftools mpileup --ignore-RG --threads 64 -R CAST.SNPs.superset.ENSEMBL.vcf.gz -Ou -a AD,DP --max-depth 8000 -f /huxley.hdd/reference_genomes/mouse/GRCm38.mgp/ref/GRCm38_68.fa bam/clonal_DNAseq/!(Undetermined).bam | bcftools call --threads 64 -mv -Oz -o out_dna/cnv/clonal_DNAseq/variants_cast.vcf.gz
bcftools index out_dna/cnv/clonal_DNAseq/variants_cast.vcf.gz
bcftools isec -Oz -o out_dna/cnv/clonal_DNAseq/variants_cast_matched.vcf.gz out_dna/cnv/clonal_DNAseq/variants_cast.vcf.gz CAST.SNPs.superset.ENSEMBL.vcf.gz -n =2 -w 1

(echo -e 'CHROM POS REF ALT' bam/clonal_DNAseq/!(Undetermined).bam  && bcftools query -f '%CHROM %POS %REF %ALT [%AD ]\n' out_dna/cnv/clonal_DNAseq/variants_cast_matched.vcf.gz) | tr ' '  '\t' | sed 's/\t$//' | gzip > out_dna/cnv/clonal_DNAseq/variants_cast_matched_counts.tsv.gz

##