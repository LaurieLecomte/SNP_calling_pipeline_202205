#!/bin/bash

# VARIABLES
OUT_DIR="05_calls"
MERGED_DIR="06_merged"
MERGED_VCF="$MERGED_DIR/merged.vcf.gz"
FILT_DIR="07_filtered"
FILT_VCF="$FILT_DIR/merged_NS30_AD1_2all.vcf.gz"

PCA_DIR="pca"
RDA_DIR="rda"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

CPU=6

# LOAD REQUIRED MODULES
module load R/4.0


#Rscript 01_scripts/RDA_pop.r
Rscript 01_scripts/RDA_pop.r $PCA_DIR/"$(basename -s .vcf.gz $FILT_VCF)".geno_mat $ID_SEX_POP $CPU

# srun -p medium -c 6 --mem=300G --time=7-00:00:00 -J RDA_pop -o log/RDA_pop_%j.log /bin/sh 01_scripts/RDA_pop.sh &