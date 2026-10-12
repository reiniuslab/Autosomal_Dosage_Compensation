#!/bin/bash

# Run quantification for set1, set2 was quantified by Xpress genomics
ulimit -n 5000
zUMIs/zUMIs.sh -y zUMIs_clonal_xpress_1.yaml
gzip -k out_mRNA/clonal_xpress_1/clonal_xpress_1.final_annot.gtf

# Allelic expression for both sets
nohup Rscript get_variant_overlap_CAST_ALv3.R --yaml zUMIs_clonal_xpress_1.yaml --vcf CAST.SNPs.superset.ENSEMBL.vcf.gz
nohup Rscript get_variant_overlap_CAST_ALv3.R --yaml zUMIs_clonal_xpress_2.yaml --vcf CAST.SNPs.superset.ENSEMBL.vcf.gz
