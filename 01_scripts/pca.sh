#!/bin/bash

# Launch with srun from 01_SNP directory :
# srun -c 1 --mem 50G -p medium --time 7-00:00 -J pca_SNPs -o log/pca_SNPs_%j.log /bin/sh 01_scripts/pca.sh &

# VARIABLES
OUT_DIR="05_calls"
MERGED_DIR="06_merged"
MERGED_VCF="$MERGED_DIR/merged.vcf.gz"
FILT_DIR="07_filtered"
FILT_VCF="$FILT_DIR/merged_NS30_AD1_2all.vcf.gz"

PCA_DIR="08_pca"

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

# LOAD REQUIRED MODULES
module load bcftools/1.15
module load vcftools
module load htslib/1.8

# 0. Create PCA dir if it does not already exist
if [[ ! -d $PCA_DIR ]]
then
  mkdir $PCA_DIR
fi

# 1. Generate genotype matrix using --012
vcftools --gzvcf $FILT_VCF --012 --out $PCA_DIR/"$(basename -s .vcf.gz $FILT_VCF)".geno_mat
	
## ".012" contains the genotypes of each individual on a separate line. Genotypes are represented as 0, 1 and 2, where the number represent that number of non-reference alleles. Missing genotypes are represented by -1. 
## ".012.indv" details the individuals included in the main file. 
## ".012.pos" details the site locations included in the main file.

# 2. Run PCA in R using script pca.R
MAT_file="$PCA_DIR/"$(basename -s .vcf.gz $FILT_VCF)".geno_mat.012"
Rscript ./01_scripts/pca.r "$MAT_file" "$ID_SEX_POP"
