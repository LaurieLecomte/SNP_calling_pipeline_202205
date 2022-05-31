#!/bin/bash

# Filter SNPs 
# Run this script from SNP_calling_pipeline_202205 directory using :
# srun -c 1 --mem=10G -p medium --time=7-00:00:00 -J 03_filter -o log/03_filter_%j.log 01_scripts/03_filter.sh &

# VARIABLES
OUT_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

# LOAD REQUIRED MODULES
module load bcftools/1.12
module load vcftools
module load htslib/1.8


# 1. Filter with bcftools
bcftools filter -e 'MQ < 30' $MERGED_DIR/merged.vcf.gz -Ov > $FILT_DIR/filtered.tmp.vcf
## print number of filtered sites
zgrep -v ^\#\# $FILT_DIR/filtered.tmp.vcf | wc -l

# 2. Filter with vcftools
vcftools --gzvcf $FILT_DIR/filtered.tmp.vcf \
    --minQ 30 \
    --minGQ 20 \
    --minDP 5 \
    --mac 2 \
    --max-alleles 2 \
    --max-missing 0.7 \
    --maf 0.05 \
    --recode \
    --stdout > $FILT_DIR/filtered.vcf

# 3. Compress and tabix
bgzip $FILT_DIR/filtered.vcf
tabix -p vcf $FILT_DIR/filtered.vcf.gz

# For SNPs AND indels :
# Filter with bcftools
# bcftools filter -e 'MQ < 30' $MERGED_DIR/merged_SNPs_indels.vcf.gz -Ov > $FILT_DIR/filtered_SNPs_indels.tmp.vcf
## print number of filtered sites
#zgrep -v ^\#\# $FILT_DIR/filtered_SNPs_indels.tmp.vcf | wc -l

# 2. Filter with vcftools
#vcftools --gzvcf $FILT_DIR/filtered_SNPs_indels.tmp.vcf \
#    --minQ 30 \
#    --minGQ 20 \
#    --minDP 5 \
#    --mac 2 \
#    --max-alleles 2 \
#    --max-missing 0.7 \
#    --maf 0.05 \
#    --recode \
#    --stdout > $FILT_DIR/filtered_SNPs_indels.vcf

# 3. Compress and tabix
#bgzip $FILT_DIR/filtered_SNPs_indels.vcf
#tabix -p vcf $FILT_DIR/filtered_SNPs_indels.vcf.gz