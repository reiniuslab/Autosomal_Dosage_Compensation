#!/bin/bash

# Run quantification
zUMIs/zUMIs.sh -y zUMIs_minibulk_SS2_p0.yaml
gzip -k out_mRNA/minibulk_SS2_p0/minibulk_p0.final_annot.gtf

# Allelic expression
Rscript get_variant_overlap_CAST_ALv3.R --yaml zUMIs_minibulk_SS2_p0.yaml --vcf CAST.SNPs.superset.ENSEMBL.vcf.gz
