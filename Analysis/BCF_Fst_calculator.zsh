#built by N. Lara
#updated 2021-5-14

#set working directory of project
cd "Desktop/exome_project"
#create folder to put output in
mkdir "output/fst"
dir="output/Processed_subsets"
dirout="output/fst"
#set input master .bcf file to pull from
vcf_in="Data/variance_files/rename_imp.vcf.gz"

#region=("North" "East" "Central" "PNWWinter" "PNWSpring")
region=("North" "PNWSpring" "PNWWinter" "East" "Central")

fst_maker () {
  pop1="${dir}/${1}_subset_samples.csv"
  pop2="${dir}/${2}_subset_samples.csv"
  echo $pop1
  echo $pop2
  pops="${1}_${2}_fst"
  echo $pops
  vcftools --gzvcf ${vcf_in} \
    --weir-fst-pop $pop1 \
    --weir-fst-pop $pop2 \
    --stdout > output
    #--out $pops
  return $output
}

#loop through regions, create region v total population
for i in "${region[@]}"; do
  subpop="${dir}/${i}_subset_samples.csv"
  bigpop="${dir}/ACROSS_subset_samples.csv"
  comm -13 $subpop $bigpop > "Processed_subsets/comp_sub.csv"
  across="${dir}/comp_sub.csv"
  #calculate weir fst value between the populations. Output
  vcftools --gzvcf ${vcf_in} \
    --weir-fst-pop $subpop \
    --weir-fst-pop $across \
    --out "${dirout}/${i}"
done

#loop through regions, compare disparate regions
reg_len=$((${#region}-1))
len=$((reg_len+1))
reg_com=()
for i in $(eval echo "{1..$reg_len}"); do
  j=$((i+1))
  for k in $(eval echo "{$j..$len}"); do
    reg_com+=("${region[$i]}/${region[$k]}")
    pop1="${dir}/${region[$i]}_subset_samples.csv"
    pop2="${dir}/${region[$k]}_subset_samples.csv"
    vcftools --gzvcf ${vcf_in} \
      --weir-fst-pop $pop1 \
      --weir-fst-pop $pop2 \
      --out "${dirout}/${region[$i]}_${region[$k]}"
  done
done
