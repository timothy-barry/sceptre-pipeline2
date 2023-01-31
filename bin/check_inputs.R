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
undercover <- as.logical(args[10]) # undercover analysis
n_pairs_to_sample <- as.integer(args[11]) # n pairs to sample
undercover_group_size <- as.integer(args[12]) # number of NTCs in undercover groups

# load multimodal ODM and pairs
mm_odm <- ondisc::read_multimodal_odm(c(response_odm_fp, grna_odm_fp), multimodal_metadata_fp)
pairs_to_analyze <- readRDS(pairs_to_analyze_fp)

#################################
# CARRY OUT THE PRELIMINARY STEPS
#################################
# step 1: check the inputs (called for side effects)
sceptre3:::check_pipeline_inputs(mm_odm, pairs_to_analyze, response_modality_name,
                                grna_modality_name, grna_group_column_name, undercover) |> invisible()

# step 2: process multimodal odm
mm_odm <- sceptre3:::process_multimodal_odm(mm_odm, response_modality_name,
                                           grna_modality_name, grna_group_column_name)
gc() |> invisible()

# step 3: if in low MOI, update the mm_odm, obtaining the index of each gRNA group
if (mm_odm@global_misc$moi == "low") {
  mm_odm <- sceptre3:::assign_grnas_to_cells_lowmoi(mm_odm = mm_odm, undercover = undercover)
}

# step 4: construct the negative control pairs (if undercover) # (IMPROVE)
if (undercover) pairs_to_analyze <- sceptre3:::construct_negative_control_pairs(mm_odm, n_pairs_to_sample,
                                                                                undercover_group_size)

# step 5: perform pairwise quality control, filtering for pairs with a sufficiently large sample size
# (IMPLEMENT)

# step 6: order pairs to analyze, and assign pair and response IDs
sceptre3:::assign_pod_ids(pairs_to_analyze, response_pod_size, pair_pod_size)
saveRDS(mm_odm, "mm_odm_new.rds")
