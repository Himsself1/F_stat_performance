#! /bin/bash

MASTER_OUT_FOLDER="/home/stefanos/simulations/inference_estimation/"
MODEL_FOLDER="no_migration_constant_size"
VCF_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/vcf/"
EIGENSTRAT_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/eig/"

mkdir -p $MASTER_OUT_FOLDER
mkdir -p $VCF_FOLDER
mkdir -p $EIGENSTRAT_FOLDER


# * msprime scripts

python3 ../msprime_scripts/msprime_no_migration.py \
    -out_folder $VCF_FOLDER \
    -name no_migration_constant_size \
    -how_many 1
    
## Need to record output folders

LIST_OF_FILES=$(find $VCF_FOLDER -type f -name "*.vcf" -exec readlink -f {} \;)

# * Convert vcf output to eigenstrat
# "${LIST_OF_FILES[@]}"
for file in ${LIST_OF_FILES[@]}; do
    prefix=$(basename "$file")
    prefix=${prefix%.vcf}
    python ../utility_scripts/gdc/vcf2eigenstrat.py -v $file -o $EIGENSTRAT_FOLDER$prefix
done
# --renameScaff
# * Run ADMIXTOOLS and other Inference software
