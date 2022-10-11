#!/bin/bash

# Filter SNPs 
# Run this script from SNP_calling_pipeline_202205 directory using :
# srun -c 2 --mem=10G -p medium --time=7-00:00:00 -J 03_filter -o log/03_filter_%j.log 01_scripts/03_filter.sh &

# VARIABLES
OUT_DIR="05_calls"

MERGED_DIR="06_merged"
MERGED_VCF="$MERGED_DIR/merged.vcf.gz"

FILT_DIR="07_filtered"

CPU=2

# TO DO : auto-format file names with pre-defined variables : 
#MIN_DP=1
#N_SAMPLES=$(less 02_infos/bam_list.txt | wc -l)
#NS=$(echo "($(less 02_infos/bam_list.txt | wc -l) * 0.5)" | bc -l)
#MIN_MAF=0.05
#MAX_MAF=0.95

# LOAD REQUIRED MODULES
module load bcftools/1.15
#module load vcftools
module load htslib/1.8

# 1. Add tags
bcftools +fill-tags $MERGED_VCF -Oz -- -t all > $MERGED_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_tagged.vcf.gz

# Filter with same criteria as SVs
# 2. Filter for depth > 1 (FORMAT/AD > 1) in genotyped samples, where at least 50% of samples have been genotyped
bcftools filter -i "N_PASS(GT!='mis' & FMT/AD>0) > 30" $MERGED_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_tagged.vcf.gz -O z --threads $CPU > $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1.vcf.gz

# 3. Filter for genotyped in >50% of samples
#bcftools filter -i 'INFO/NS >= 30' $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_DP1.vcf.gz -O z --threads $CPU > $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_DP1_NS30.vcf.gz

# 4. Filter for max number of alleles = 2
bcftools view --max-alleles 2 $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1.vcf.gz -O z --threads $CPU > $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1_2all.vcf.gz

# 5. Filter for maf > 0.05 and < 0.95
#bcftools filter -i 'INFO/MAF >= 0.05' && 'INFO/MAF <= 0.95' $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1_2all.vcf.gz -O z --threads $CPU > $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1_2all_maf0.05.vcf.gz


# 6. Compress and tabix
#bgzip $FILT_DIR/filtered.vcf
#bgzip $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1_2all.vcf.gz
tabix -p vcf $FILT_DIR/"$(basename -s .vcf.gz $MERGED_VCF)"_NS30_AD1_2all.vcf.gz



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

# 3. Seperate SNPs and indels 
# bcftools filter -i "INFO/INDEL=1" $FILT_DIR/filtered_SNPs_indels.vcf > $FILT_DIR/filtered_indels.vcf
# bcftools filter -e "INFO/INDEL=1" $FILT_DIR/filtered_SNPs_indels.vcf > $FILT_DIR/filtered_SNPs.vcf

# 4. Compress and tabix
#bgzip $FILT_DIR/filtered_SNPs_indels.vcf
#tabix -p vcf $FILT_DIR/filtered_SNPs_indels.vcf.gz

#bgzip $FILT_DIR/filtered_SNPs.vcf
#tabix -p vcf $FILT_DIR/filtered_SNPs.vcf.gz

#bgzip $FILT_DIR/filtered_indels.vcf
#tabix -p vcf $FILT_DIR/filtered_indels.vcf.gz