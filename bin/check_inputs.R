#!/usr/bin/env Rscript

#######
# SETUP
#######

# get the command line args
args <- commandArgs(trailingOnly = TRUE)
multimodal_metadata_fp <- args[1] # multimodal metadata fp
response_odm_fp <- args[2] # response ODM backing file
grna_odm_fp <- args[3] # grna ODM backing file
pairs_to_analyze_fp <- args[4] # pairs to analyze df
response_modality_name <- args[5] # response modality name
grna_modality_name <- args[6] # grna modality name
grna_group_column_name <- args[7] # grna group column name
response_pod_size <- max(as.integer(args[8]), 2) # response pod size
pair_pod_size <- max(as.integer(args[9]), 2) # grna pod size
n_pairs_to_sample <- as.integer(args[10]) # n pairs to sample

# load multimodal ODM and pairs
mm_odm <- ondisc::read_multimodal_odm(c(response_odm_fp, grna_odm_fp), multimodal_metadata_fp)
pairs_to_analyze <- readRDS(pairs_to_analyze_fp)

#################################
# CARRY OUT THE PRELIMINARY STEPS
#################################

# step 1: check the inputs (called for side effects)
check_pipeline_inputs(mm_odm, pairs_to_analyze, response_modality_name, grna_modality_name, grna_group_column_name)

# step 2: 