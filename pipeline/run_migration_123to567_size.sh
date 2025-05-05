#! /bin/bash

RECOMB_RATE=1.25e-7
SEQ_LENGTH=1e+7

MASTER_OUT_FOLDER="/home/stefanos/new_storage/inference_estimation/"
MODEL_FOLDER="sequencies/migration_123to567_mig_02_constant_size_rec_${RECOMB_RATE}_seq_${SEQ_LENGTH}"
VCF_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/vcf/"
EIGENSTRAT_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/eig/"
STATISTICS_FOLDER=$MASTER_OUT_FOLDER"statistics"

SCRIPT_FOR_R_ANALYSIS="/home/stefanos/F_stat_performance/qpadm_inference/best_populations_plot_funtions.R"

mkdir -p $MASTER_OUT_FOLDER
mkdir -p $VCF_FOLDER
mkdir -p $EIGENSTRAT_FOLDER
mkdir -p $STATISTICS_FOLDER

# * msprime scripts

python3 ../msprime_scripts/msprime_migration_123to567.py \
    -out_folder $VCF_FOLDER \
    -name "migration_123to567_mig_02_constant_size_rec_${RECOMB_RATE}_seq_${SEQ_LENGTH}" \
    -how_many 100 \
    -rec $RECOMB_RATE \
    -mig_ratio 0.2
    
## Need to record output folders

LIST_OF_FILES=$(find $VCF_FOLDER -type f -name "*.vcf" -exec readlink -f {} \;)

# * Convert vcf output to eigenstrat

for file in ${LIST_OF_FILES[@]}; do
    prefix=$(basename "$file")
    prefix=${prefix%.vcf}
    python ../utility_scripts/gdc/vcf2eigenstrat.py -v $file -o $EIGENSTRAT_FOLDER$prefix
done

## Need to change CentiMorgan distance of the converted .snp files

LIST_OF_snps=$(find $EIGENSTRAT_FOLDER -type f -name "*.snp" -exec readlink -f {} \;)
for file in ${LIST_OF_snps[@]}; do
    ## Multiply position with recombination rate
    awk '{ $3 = $4 * $RECOMB_RATE } 1' $file > $file".tmp"
    mv $file".tmp" $file
done


# * Run ADMIXTOOLS and other Inference software

## Runs R wrapper for qpadm inference.

Rscript ../qpadm_inference/qpadm_inference.R \
	$EIGENSTRAT_FOLDER\
	$STATISTICS_FOLDER\
	$SCRIPT_FOR_R_ANALYSIS

