#! /bin/bash

MASTER_OUT_FOLDER="/home/stefanos/new_storage/inference_estimation/"
MODEL_FOLDER="sequencies/no_migration_constant_size_scale_2"
VCF_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/vcf/"
EIGENSTRAT_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/eig/"
STATISTICS_FOLDER=$MASTER_OUT_FOLDER"statistics"

mkdir -p $MASTER_OUT_FOLDER
mkdir -p $VCF_FOLDER
mkdir -p $EIGENSTRAT_FOLDER
mkdir -p $STATISTICS_FOLDER

# * msprime scripts

python3 ../msprime_scripts/msprime_no_migration.py \
    -out_folder $VCF_FOLDER \
    -name no_migration_constant_size \
    -how_many 100 \
    -scale 2

LIST_OF_FILES=$(find $VCF_FOLDER -type f -name "*.vcf" -exec readlink -f {} \;)

# * Convert vcf output to eigenstrat
"${LIST_OF_FILES[@]}"
for file in ${LIST_OF_FILES[@]}; do
    prefix=$(basename "$file")
    prefix=${prefix%.vcf}
    python ../utility_scripts/gdc/vcf2eigenstrat.py -v $file -o $EIGENSTRAT_FOLDER$prefix
done

## Need to change CentiMorgan distance of the converted .snp files

LIST_OF_snps=$(find $EIGENSTRAT_FOLDER -type f -name "*.snp" -exec readlink -f {} \;)
for file in ${LIST_OF_snps[@]}; do
    ## Multiply position with recombination rate
    awk '{ $3 = $4*1.25e-7 } 1' $file > $file".tmp"
    # awk '{ $3 = (50/100)*log(1/(1-2*(1e-8)*$4)) } 1' $file > $file".tmp"
    mv $file".tmp" $file
done


# * Run ADMIXTOOLS and other Inference software

## Runs R wrapper for qpadm inference.
# Rscript ../qpadm_inference/qpadm_inference.R $EIGENSTRAT_FOLDER $STATISTICS_FOLDER
