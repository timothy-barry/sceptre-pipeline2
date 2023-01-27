#!/bin/bash

# Limit NF driver to 4 GB memory
export NXF_OPTS="-Xms500M -Xmx4G"

##########################
# REQUIRED INPUT ARGUMENTS
##########################
source ~/.research_config
gasp_data_dir=$LOCAL_SCEPTRE2_DATA_DIR"data/gasperini/at_scale/"
# 1) multimodal metadata file
multimodal_metadata_fp=$gasp_data_dir"multimodal_metadata.rds"
# 2) response ODM
response_odm_fp=$gasp_data_dir"gene/matrix.odm"
# 3) grna ODM
grna_odm_fp=$gasp_data_dir"grna_expression/matrix.odm"
# 4) undercover
undercover="true"
# 5) grna modality name
grna_modality_name="grna_expression"
# 6) response modality name
response_modality_name="gene"

##############
# OUTPUT FILE:
##############
result_fp="~/sceptre_result_highmoi_example.rds"

###############
# OPTIONAL ARGS
###############
# formula, threshold, B, side, n_pairs_to_sample, response_pod_size, grna_group_pod_size, pair_pod_size are optional args
n_pairs_to_sample=20
undercover_group_size=2
response_pod_size=5
pair_pod_size=10

#################
# Invoke pipeline
#################
nextflow run ../main.nf \
 --multimodal_metadata_fp  $multimodal_metadata_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --result_fp $result_fp \
 --grna_modality_name $grna_modality_name \
 --response_modality_name $response_modality_name \
 --response_pod_size $response_pod_size \
 --pair_pod_size $pair_pod_size \
 --result_fp $result_fp \
 --undercover $undercover \
 --undercover_group_size $undercover_group_size \
 --n_pairs_to_sample $n_pairs_to_sample
