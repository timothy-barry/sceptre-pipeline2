#!/usr/bin/env Rscript

source("~/research_code/sceptre-pipeline2/bin/run_response_precomputation_funct.R")

# 0. get the command line args
args <- commandArgs(trailingOnly = TRUE)
multimodal_metadata_fp <- args[1]
response_odm_fp <- args[2]
grna_odm_fp <- args[3]
response_to_pod_id_map_fp <- args[4]
curr_pod <- as.integer(args[5])
use_mass <- as.logical(args[6])

# 1. load
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(response_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = multimodal_metadata_fp)
response_to_pod_id_map <- readRDS(response_to_pod_id_map_fp)

# 2. setup
response_modality <- mm_odm |> ondisc::get_modality("response")
curr_responses <- get_response_ids_from_pod(response_to_pod_id_map, curr_pod)
covariate_matrix <- mm_odm |> ondisc::get_cell_covariates()
low_moi <- mm_odm@global_misc$moi == "low"

# 3. if in low moi, obtain the idxs of the NT gRNAs
if (low_moi) {
  nt_idxs <- mm_odm@modalities$grna@misc$grna_group_idxs[["non-targeting"]]
  covariate_matrix <- covariate_matrix[nt_idxs,]
}
rm(mm_odm, response_to_pod_id_map); invisible(gc())

# 4. run the precomputations
out <- sceptre3:::run_response_precomputation_high_level(response_modality = response_modality,
                                                         covariate_matrix = covariate_matrix,
                                                         curr_responses = curr_responses,
                                                         nt_idxs = if (low_moi) nt_idxs else NULL,
                                                         use_mass = use_mass)

saveRDS(out$precomp_coef_df, file = "precomp_coef_df.rds")
saveRDS(out$precomp_info_df, file = "precomp_info_df.rds")
