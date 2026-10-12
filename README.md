# Autosomal Dosage Compensation in Aneuploidy and Cancer Involves Transcriptional Bursting Regulation
*Natali Papanicolaou<sup>1†</sup>, Sebastian Wettersten<sup>1†</sup>, Guilherme Maia<sup>1</sup>, Antonio Lentini<sup>1,2,3#</sup> and Björn Reinius<sup>1#</sup>*

<sup>1</sup>Department of Medical Biochemistry and Biophysics, Karolinska Institutet, Stockholm, Sweden.

<sup>2</sup>Department of Biomedical and Clinical Sciences, Linköping University, Linköping, Sweden.

<sup>2</sup>Science for Life Laboratory (SciLifeLab), Linköping University, Linköping, Sweden.

##

This repository contains scripts and additional data needed to reproduce findings described in [PAPER](), [PRE-PRINT](https://doi.org/10.1101/2025.09.18.677044).

Raw and processed sequencing data is available through ArrayExpress: [[DNA-seq]](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-15634/), [[Bulk RNA-seq]](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-15633/) and [[Allelic scRNA-seq]](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-15744/), and proteomics data is available through PRIDE: [[Proteomics]](https://www.ebi.ac.uk/pride/archive/projects/PXD068580). Additional processed datasets are available in Zenodo [[Datasets]](https://doi.org/10.5281/zenodo.23055741).

Validated CAST/EiJ SNPs were obtained from Hagemann-Jensen et al., 2020 doi:10.1038/s41587-020-0497-0, [data](https://github.com/sandberg-lab/Smart-seq3/blob/5d5938475039f5c98d0d94faf89db917f66fe8ac/allele_level_expression/CAST.SNPs.validated.vcf.gz).

UMI-spike-in sequences were obtained from Ziegenhain & Hendricks et al., 2022 doi:10.1038/s41592-022-01446-x, [data](https://github.com/sandberg-lab/molecularSpikes/blob/d3d85cbb896941103a29827fbad047d26fb8d609/fasta_reference/molecularSpikes_complexset_5p.fa).

*txburst.py* is from Larsson et al., 2019 doi:10.1038/s41586-018-0836-1, [code](https://github.com/sandberg-lab/txburst/blob/6b6525f951f04833a0c59b2926a7e66836209bf0/txburstML.py).

## Abstract
Aneuploidy, resulting from segmental or whole-chromosome gains and losses, is generally detrimental to organismal fitness yet is strongly associated with tumorigenesis. While sex chromosome aneuploidies are relatively well tolerated due to dedicated allelic dosage compensation mechanisms, the existence of comparable regulatory processes for autosomes remain debated. Here, we investigate transcriptional responses to autosomal aneuploidy and identify widespread dosage compensation across affected chromosomes. Using high-sensitivity, allele-resolved single-cell RNA sequencing of monoclonal cell expansions spanning diverse ploidies and degrees of aneuploidy, we uncover consistent transcriptional compensation of both chromosomal gains and losses. This response operates in a region-specific manner across both complete and segmental aneuploidies, revealing an unexpected flexibility in dosage regulation. We show that chromosome losses are compensated through increased transcriptional burst frequency from the remaining allele. Complementary proteomics analyses demonstrate additional buffering at the protein level, resulting in a substantial, though still incomplete, restoration of stoichiometric balance across autosomal aneuploid chromosomes for both copy number gains and losses. Analyses of cancer transcriptomes confirm that compensatory mechanisms are active in human tumours. Together, our findings establish dosage compensation as a genome-wide and dynamic response to gene dose imbalance. This mechanism extends beyond sex chromosomes, supports transcriptional homeostasis, and represents a fundamental and evolutionarily conserved mode of gene regulation active across species and in aneuploid cancer cells.
##

### Session info
*R version 4.4.2 (2024-10-31)*

| Package | Version | Source |
| :------ | :------ | :----- |
| `ACE` | 1.17.0 | Github (tgac-vumc/ACE) |
| `HMMcopy` | 1.48.0 | Bioconductor |
| `scran` | 1.34.0 | Bioconductor |
| `scater` | 1.34.1 | Bioconductor |
| `scuttle` | 1.16.0 | Bioconductor |
| `ggplot2` | 4.0.2 | Bioconductor |
| `patchwork` | 1.3.0 | CRAN |
| `rtracklayer` | 1.66.0 | Bioconductor |
| `magrittr` | 2.0.3 | CRAN |
| `data.table` | 1.17.0 | CRAN |
| `UMIcountR` | 0.5.0 | Github (cziegenhain/UMIcountR) |
| `TOSTER` | 0.8.6 | CRAN |
