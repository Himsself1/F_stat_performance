#! /bin/bash

# RECOMB_RATE=(1.25e-8 1.25e-9)
# MUT_RATE=(1.25e-8 1.25e-9)

RECOMB_RATE=1.25e-8
MUT_RATE=1.25e-8

SEQ_LENGTH=249e+6

MASTER_OUT_FOLDER="/home/stefanos/new_storage/inference_estimation/"
PATH_TO_CUSTOM_PYTHON_FUNCTION="/home/stefanos/F_stat_performance/msprime_scripts/ts_to_eigenstrat.py"
SCRIPT_FOR_R_ANALYSIS="/home/stefanos/F_stat_performance/qpadm_inference/best_populations_plot_funtions.R"

MODEL_FOLDER="harney_model_mut_${MUT_RATE}_rec_${RECOMB_RATE}_seq_${SEQ_LENGTH}"
EIGENSTRAT_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/eig/"
OUTPUT_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER

mkdir -p $OUTPUT_FOLDER
mkdir -p $EIGENSTRAT_FOLDER

# * msprime scripts

# python3 ../msprime_scripts/harney_et_al_demography.py \
# 	-out_folder $MASTER_OUT_FOLDER \
#     -name "harney_model_mut_${MUT_RATE}_rec_${RECOMB_RATE}_seq_${SEQ_LENGTH}" \
#     -how_many 100 \
#     -rec $RECOMB_RATE \
# 	-mut $MUT_RATE \
#     -seq_length $SEQ_LENGTH \
# 	-nsamples 5

# ** Parallel loop
## Need to handple output folders as well!
# parallel -j 4 '
# 		 MODEL_FOLDER="sequencies/harney_model_mut_${1}_rec_${2}_seq_${SEQ_LENGTH}"
# 		 EIGENSTRAT_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/eig/"
# 		 STATISTICS_FOLDER=$MASTER_OUT_FOLDER$MODEL_FOLDER"/statistics/"
# 		 mkdir -p $MODEL_FOLDER
# 		 mkdir -p $STATISTICS_FOLDER
# 		 mkdir -p $EIGENSTRAT_FOLDER
# 		 python3 ../msprime_scripts/harney_et_al_demography.py \
# 		     -out_folder "$VCF_FOLDER" \
# 			 -name "harney_model_mut_{1}_rec_${2}_seq_${SEQ_LENGTH}" \
# 			 -how_many 100 \
# 			 -rec {2} \
# 			 -mut {1} \
# 			 -seq_length "$SEQ_LENGTH" \
# 			 -nsamples 5

# 		Rscript ../qpadm_inference/qpadm_inference.R \
# 				$EIGENSTRAT_FOLDER\
# 				$STATISTICS_FOLDER\
# 				$SCRIPT_FOR_R_ANALYSIS

# 		Rscript ../qpadm_inference/qpadm_inference_pop3.R \
# 				$EIGENSTRAT_FOLDER \
# 				$STATISTICS_FOLDER \
# 				$SCRIPT_FOR_R_ANALYSIS
# ' ::: "${MUT_RATES[@]}" ::: "${REC_RATES[@]}"

# * Run ADMIXTOOLS and other Inference software

## Runs R wrapper for qpadm inference.

# Rscript ../qpadm_inference/qpadm_inference_harney.R \
# 		$EIGENSTRAT_FOLDER \
# 		$OUTPUT_FOLDER \
# 		$SCRIPT_FOR_R_ANALYSIS

Rscript ../qpadm_inference/plot_qpadm_results.R \
		"${OUTPUT_FOLDER}/stats" \
		"${OUTPUT_FOLDER}/plots" \
		$SCRIPT_FOR_R_ANALYSIS \
		$MODEL_FOLDER


# Rscript ../qpadm_inference/qpadm_inference_pop3.R \
# 		$EIGENSTRAT_FOLDER \
# 		$STATISTICS_FOLDER \
# 		$SCRIPT_FOR_R_ANALYSIS
