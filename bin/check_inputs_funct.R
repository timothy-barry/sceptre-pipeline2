check_pipeline_inputs <- function(mm_odm, pairs_to_analyze, response_modality_name, grna_modality_name, grna_group_column_name) {
 
  grna_modality <- mm_odm |> ondisc::get_modality(grna_modality_name)
  response_modality <- mm_odm |> ondisc::get_modality(response_modality_name)
  grna_feature_covariates <- grna_modality |> ondisc::get_feature_covariates()
  odm_response_ids <- response_modality |> ondisc::get_feature_ids()
  
  # 1. check modality names
  modality_names <- names(mm_odm@modalities)
  if (!(response_modality_name %in% modality_names)) {
    stop(paste0("The response modality of the multimodal ondisc matrix must be named `", response_modality_name,"`. You can check the current modality names via `names(mm_odm@modalities)`, where `mm_odm` is the multimodal ondisc matrix."))
  }
  if (!(grna_modality_name %in% modality_names)) {
    stop(paste0("The grna modality of the multimodal ondisc matrix must be named `", grna_modality_name,"`. You can check the current modality names via `names(mm_odm@modalities)`, where `mm_odm` is the multimodal ondisc matrix."))
  }
  
  # 2. check that "grna_group" is a column of the feature covariate matrix of the grna modality
  if (!(grna_group_column_name %in% colnames(grna_feature_covariates))) {
    stop(paste0("The feature covariate data frame of the grna modality of the multimodal ondisc matrix must have a column named `", grna_group_column_name, "` specifying the group to which each gRNA belongs. You can view the names of the columns of this data frame via `mm_odm |> ondisc::get_modality('", grna_modality_name, "') |> ondisc::get_feature_covariates() |> colnames()`, where `mm_odm` is the multimodal ondisc matrix."))
  }
  
  # 3. check for the presence of "non-targeting" in the grna_group column
  nt_present <- "non-targeting" %in% grna_feature_covariates[[grna_group_column_name]]
  if (!nt_present) {
    stop(paste0("The string 'non-targeting' must be present in the `", grna_group_column_name, "` column of the feature covariate data frame of the grna modality of the multimodal ondisc matrix. You can access this column via `mm_odm |> ondisc::get_modality('", grna_modality_name, "') |> ondisc::get_feature_covariates() |> dplyr::pull('", grna_group_column_name ,"')`, where `mm_odm` is the multimodal ondisc matrix."))
  }
  
  # 4. if the pairs data frame has been specified...
  if (!is.null(pairs_to_analyze)) {
    # 4.a verify that response id is present as a column
    if (!("response_id" %in% colnames(pairs_to_analyze))) {
      stop("The `pairs_to_analyze` data frame must contain a column called `response_id`.")
    }
    # 4b. verify that grna_group is present as a column
    if (!("grna_group" %in% colnames(pairs_to_analyze))) {
      stop("The `pairs_to_analyze` data frame must contain a column called `grna_group`.")
    }
    # 4c. check that the response ids in the pairs_to_analyze data frame are a subset of the ids of the response modality
    pairs_response_ids <- as.character(pairs_to_analyze$response_id)
    if (!all(pairs_response_ids %in% odm_response_ids)) {
      stop(paste0("The column `response_id` of the `pairs_to_analyze` data frame must be a subset of the feature IDs of the response modality of the multimodal ondisc matrix. You can access the latter via `mm_odm |> ondisc::get_modality('", response_modality_name, "') |> ondisc::get_feature_ids()`, where mm_odm is the multimodal ondisc matrix. Update the `response_id` column of the `pairs_to_analyze` data frame to ensure that this requirement is satisfied."))
    }
    # 4d. check that the grna groups in the pairs_to_analyze data frame are a subset of the grna groups in the grna modality
    odm_grna_groups <- grna_feature_covariates[[grna_group_column_name]]
    pairs_grna_groups <- as.character(pairs_to_analyze$grna_group)
    if (!all(pairs_grna_groups %in% odm_grna_groups)) {
      stop(paste0("The column `grna_group` of the `pairs_to_analyze` data frame must be a subset of the gRNA groups in the gRNA modality of the multimodal ondisc matrix. You can access the latter via `mm_odm |> ondisc::get_modality('", grna_modality_name, "') |> ondisc::get_feature_covariates() |> dplyr::pull('", grna_group_column_name ,"')`, where `mm_odm` is the multimodal ondisc matrix. Update the `grna_group` column of the `pairs_to_analyze` data frame so that this requirement is satisfied."))
    }
  }
  
  # 5. verify that both modalities have unique response ids
  grna_ids <- grna_modality |> ondisc::get_feature_ids()
  if (length(unique(grna_ids)) != length(grna_ids) || length(grna_ids) != nrow(grna_modality)) stop("The gRNA modality of the multimodal ondisc matrix must have one unique ID per individual gRNA.")
  if (length(unique(odm_response_ids)) != length(odm_response_ids) || length(odm_response_ids) != nrow(response_modality)) stop("The response modality of the multimodal ondisc matrix must have one unique ID per individual response.")
   
  # 6. check that `formula` is present in the global_misc field of the `mm_odm`
  global_misc_names <- names(mm_odm@global_misc)
  if (!("formula" %in% global_misc_names)) {
    stop("A formula must be specified. Add a formula via `mm_odm@global_misc$formula <- my_formula`, where `my_formula` is a formula and `mm_odm` is the multimodal ondisc matrix.")
  }
  
  # 7. likewise check that `moi` is present in `mm_odm` and that this field is set to 'low' or 'high'
  if (!("moi" %in% global_misc_names)) {
    stop("The MOI of the data must be specified. Specify the moi via `mm_odm@global_misc$moi <- my_moi`, where `my_moi` is a string indicating the MOI of the data (either 'low' or 'high') and `mm_odm` is the multimodal ondisc matrix.")
  }
  if (!(mm_odm@global_misc$moi %in% c("low", "high"))) {
    stop("The MOI of the data must be either 'low' or 'high'. Access the MOI via `mm_odm@global_misc$moi`, where `mm_odm` is the multimodal ondisc matrix.")
  }
  return(NULL)
}
