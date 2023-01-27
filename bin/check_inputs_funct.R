check_pipeline_inputs <- function(mm_odm, pairs_to_analyze, response_modality_name, grna_modality_name, grna_group_column_name, undercover) {
 
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
  
  # 4. if we are doing a discovery analysis...
  if (!undercover) {
    # 4.a verify that a pairs data frame has been supplied
    if (is.null(pairs_to_analyze)) stop("A `pairs_to_analyze` data frame must be supplied when conducting a discovery analysis.")
    # 4.b verify that response id is present as a column
    if (!("response_id" %in% colnames(pairs_to_analyze))) {
      stop("The `pairs_to_analyze` data frame must contain a column called `response_id`.")
    }
    # 4.c verify that grna_group is present as a column
    if (!("grna_group" %in% colnames(pairs_to_analyze))) {
      stop("The `pairs_to_analyze` data frame must contain a column called `grna_group`.")
    }
    # 4.d check that the response ids in the pairs_to_analyze data frame are a subset of the ids of the response modality
    pairs_response_ids <- as.character(pairs_to_analyze$response_id)
    if (!all(pairs_response_ids %in% odm_response_ids)) {
      stop(paste0("The column `response_id` of the `pairs_to_analyze` data frame must be a subset of the feature IDs of the response modality of the multimodal ondisc matrix. You can access the latter via `mm_odm |> ondisc::get_modality('", response_modality_name, "') |> ondisc::get_feature_ids()`, where mm_odm is the multimodal ondisc matrix. Update the `response_id` column of the `pairs_to_analyze` data frame to ensure that this requirement is satisfied."))
    }
    # 4.e check that the grna groups in the pairs_to_analyze data frame are a subset of the grna groups in the grna modality
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
   
  # Also ensure that the amperstand symbol (&) is absent from the grna ids
  problematic_grna_ids <- grep(pattern = "&", x = grna_ids)
  if (length(problematic_grna_ids) >= 1) {
    stop(paste0("The ampersand character (&) cannot be present in the gRNA IDs. The following gRNA IDs contain an ampersand: ", paste0(grna_ids[problematic_grna_ids], collapse = ", ")))
  }
  
  # 6. check that `formula` is present in the global_misc field of the `mm_odm`; also, check that there are no offsets
  global_misc_names <- names(mm_odm@global_misc)
  if (!("formula" %in% global_misc_names)) {
    stop("A formula must be specified. Add a formula via `mm_odm@global_misc$formula <- my_formula`, where `my_formula` is a formula and `mm_odm` is the multimodal ondisc matrix.")
  }
  if (grepl("offset", as.character(mm_odm@global_misc$formula)[2])) stop("Offsets are not currently supported in formula objects.")
  
  # 7. likewise check that `moi` is present in `mm_odm` and that this field is set to 'low' or 'high'
  if (!("moi" %in% global_misc_names)) {
    stop("The MOI of the data must be specified. Specify the moi via `mm_odm@global_misc$moi <- my_moi`, where `my_moi` is a string indicating the MOI of the data (either 'low' or 'high') and `mm_odm` is the multimodal ondisc matrix.")
  }
  if (!(mm_odm@global_misc$moi %in% c("low", "high"))) {
    stop("The MOI of the data must be either 'low' or 'high'. Access the MOI via `mm_odm@global_misc$moi`, where `mm_odm` is the multimodal ondisc matrix.")
  }
  
  # 8. check the global cell covariate matrix for correctness
  global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()
  if (ncol(global_cell_covariates) == 0) {
    stop("The global cell covariate matrix must contain at least one column.")
  }
  
  return(NULL)
}


process_multimodal_odm <- function(mm_odm, response_modality_name, grna_modality_name, grna_group_column_name) {
  # 1. rename modalities and grna group column
  names(mm_odm@modalities)[names(mm_odm@modalities) == response_modality_name] <- "response"
  names(mm_odm@modalities)[names(mm_odm@modalities) == grna_modality_name] <- "grna"
  mm_odm@modalities$grna@feature_covariates <- dplyr::rename(mm_odm@modalities$grna@feature_covariates,
                                                             "grna_group" = !!grna_group_column_name)
  
  # 2. apply formula object to global cell covariates
  global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()
  global_cell_covariates_new <- stats::model.matrix(object = mm_odm@global_misc$formula,
                                                    data = global_cell_covariates) |> as.data.frame()
  
  # verify that the global cell covariates are OK after transformation
  for (col_name in colnames(global_cell_covariates_new)) {
    vect <- global_cell_covariates_new[[col_name]]
    if (any(vect == -Inf) || any(vect == Inf) || any(is.na(vect))) {
      stop(paste0("The column `", col_name, "` of the transformed global cell covariate matrix contains entries that are NA, -Inf, or Inf. You can access this column via `global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates(); stats::model.matrix(object = mm_odm@global_misc$formula, data = global_cell_covariates) |> as.data.frame() |> dplyr::pull('", col_name,"')`, where mm_odm is the multimodal ondisc matrix."))
    }
  }
  mm_odm@global_cell_covariates <- global_cell_covariates_new
  
  # 3. thin the multimodal odm
  mm_odm <- ondisc:::thin_multimodal_odm(mm_odm)
  return(mm_odm)
}


get_undercover_groups <- function(ntc_names, group_size, n_undercover_groups) {
  set.seed(4)
  n_ntcs <- length(ntc_names)
  total_possible_paritions <- choose(n_ntcs, group_size)
  n_undercover_groups <- min(total_possible_paritions, n_undercover_groups)
  
  my_undercover_groups <- character()
  ntc_names_copy <- ntc_names
  repeat {
    while (length(ntc_names_copy) >= group_size) {
      curr_grp <- sample(x = ntc_names_copy, size = group_size, 
                         replace = FALSE)
      curr_grp_string <- paste0(sort(curr_grp), collapse = "&")
      if (!(curr_grp_string %in% my_undercover_groups)) {
        my_undercover_groups <- c(my_undercover_groups, curr_grp_string)
      }
      ntc_names_copy <- ntc_names_copy[!(ntc_names_copy %in% curr_grp)]
    }
    if (length(my_undercover_groups) >= n_undercover_groups) break
    ntc_names_copy <- ntc_names
  }
  return(my_undercover_groups[seq(1, n_undercover_groups)])
}


##########################
# THIS FUNCTION NEEDS WORK
##########################
construct_negative_control_pairs <- function(mm_odm, n_pairs_to_sample, undercover_group_size) {
  n_genes_per_undercover_group <- 5
  
  nt_grnas <- mm_odm |> ondisc::get_modality("grna") |>
    ondisc::get_feature_covariates() |>
    dplyr::filter(grna_group == "non-targeting") |>
    row.names()
  response_ids <- mm_odm |> ondisc::get_modality("response") |> ondisc::get_feature_ids()
  undercover_grps <- get_undercover_groups(ntc_names = nt_grnas,
                                           group_size = undercover_group_size,
                                           n_undercover_groups = floor(length(nt_grnas)/undercover_group_size)) |> factor()
  # pair randomly with genes
  response_ids_samp <- sample(x = response_ids,
                              size = n_genes_per_undercover_group * length(undercover_grps),
                              replace = FALSE) |> factor()
  # return out data frame
  out <- data.frame(grna_group = rep(undercover_grps, times = n_genes_per_undercover_group),
                    response_id= response_ids_samp) |> dplyr::sample_n(n_pairs_to_sample)
}


get_id_vect <- function(df, pod_size) {
  breaks <- round(nrow(df)/pod_size)
  if (breaks >= 2) {
    as.integer(cut(seq(1, nrow(df)), breaks)) 
  } else {
    rep(1L, nrow(df))
  }
}


get_id_vect <- function(df, pod_size) {
  breaks <- round(nrow(df)/pod_size)
  if (breaks >= 2) {
    as.integer(cut(seq(1, nrow(df)), breaks)) 
  } else {
    rep(1L, nrow(df))
  }
}


write_vector <- function(file_name, vector) {
  file_con <- file(file_name)
  writeLines(as.character(vector), file_con)
  close(file_con)
}


assign_pod_ids <- function(pairs_to_analyze, gene_pod_size, pair_pod_size) {
  # create the response-to-pod ID data frame and update the pairs to analyze data frame with pod IDs
  response_pod_df <- data.frame(response_id = unique(as.character(pairs_to_analyze$response_id)))
  response_pod_df$pod_id <- get_id_vect(response_pod_df, gene_pod_size)
  pairs_to_analyze$pod_id <- get_id_vect(pairs_to_analyze, pair_pod_size)
  
  # write the data frames, alongside text files containing the integer indices, to disk
  saveRDS(response_pod_df, "response_to_pod_id_map.rds")
  saveRDS(pairs_to_analyze, "pairs_to_analyze.rds")
  write_vector("response_pods.txt", unique(response_pod_df$pod_id))
  write_vector("pair_pods.txt", unique(pairs_to_analyze$pod_id))
}
