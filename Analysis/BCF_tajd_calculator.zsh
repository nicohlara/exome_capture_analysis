#built by N. Lara
#updated 2021-5-14

#set working directory of project
cd "Desktop/exome_project/"
mkdir "output/taj_d"
dir="output/taj_d"
vcf_in="Data/variance_files/rename_imp.vcf.gz"
lines_file="Data/name_region.csv"

#region=("North" "East" "Central" "PNWWinter" "PNWSpring")
region=("North" "PNWSpring" "PNWWinter" "East" "Central")

#loop through regions, create region v total population
for i in "${region[@]}"; do
  subpop="${dir}/../Processed_subsets/${i}_subset_samples.csv"
  vcftools --gzvcf ${vcf_in} \
    --keep $subpop \
    --TajimaD 10000 \
    --out "${dir}/${i}"
done
