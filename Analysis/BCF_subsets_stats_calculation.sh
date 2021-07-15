#!/bin/bash

## Compare Variant Statistics Across Subsets
##
## This script is designed for a very specific task - it takes a VCF file as
## input, and then sequentially pulls subsets of lines out of it by using a user-
## supplied file of line names. Finally, it measures several statistics for each
## variant, including frequency of missing data, minor allele frequency, and
## total depth.
##
## grep is used to find all lines matching a specific pattern, so the line
## names must contain patterns that easily allow for grepping many at once.
##
## Note that this will calculate MAF = 0 if a particular SNP happens to be
## monomorphic in a particular subset of lines.
################################################################################

###This file was originally recieved from Brian Ward, former member of the Brown-Guedira lab
###Heavily modified by N. Lara to work with input data sets and regional filtering
################################################################################
###requires bcftools to function


###REPLACE TO LOCATION ON HOME DEVICE
cd "Desktop/exome_project"

#Input files
vcf_in="Data/variance_files/rename_v1_het01_miss08_maf003_variants_biallelic_csq.bcf"
lines_file="Data/name_region.csv"

#filter by regional
#"ACROSS" will capture all lines
region=("North" "East" "Central" "PNWWinter" "PNWSpring" "ACROSS")

mkdir "output/Processed_subsets"
dir="output/Processed_subsets"

echo
echo "Start time:"
date

## Loop through region array
for i in "${region[@]}"; do
    ## Create output header
    echo "CHROM\tPOS\tTYPE\tREF\tALT\tF_MISSING\tMAF\tDP_SUM" |
        gzip -c -f > "${dir}/${i}_subsets_stats.txt.gz"
    ## Get lines in subset and count how many there are
    if [[ "$i" == "ACROSS" ]]; then
        cut -d "," -f 1 "$lines_file" > "${dir}/${i}_subset_samples.csv"
    else
        grep "$i" "$lines_file" | cut -d "," -f 1 > "${dir}/${i}_subset_samples.csv"
    fi
    nlines=$(cat "${dir}/${i}_subset_samples.csv" | wc -l)
	## Complicated pipe to subset out the lines, calculate MAF for each SNP,
    ## print out info on each SNP, and add in subset identifying info
    bcftools view "$vcf_in" -S "${dir}/${i}_subset_samples.csv" -Ou --force-samples |
        bcftools +fill-tags -Ou -- -t MAF,F_MISSING,'DP=sum(DP)' |
        bcftools query -f '%CHROM\t%POS\t%TYPE\t%REF\t%ALT\t%INFO/F_MISSING\t%INFO/MAF\t%INFO/DP\n' |
        gzip -c -f >> "${dir}/${i}_subsets_stats.txt.gz"
done
echo
echo "End time:"
date
