// Use DSL 2
nextflow.enable.dsl = 2

// Set default values of optional parameters
params.B1 = 5000
params.B2 = 250000
params.p_thresh = 0.01
params.pairs_to_analyze_fp = "$baseDir/assets/dummy_pairs_to_analyze.rds"
params.test_stat = "full"
params.resampling_mech = "default"
params.side = "both"
params.grna_threshold = 3
params.n_pairs_to_sample = 0
params.result_fp = "$PWD/sceptre_result.rds"
params.response_pod_size = 200
params.pair_pod_size = 5000
params.undercover_n_pairs = 100000
params.undercover_group_size = 1
params.response_modality_name = "gene"
params.grna_modality_name = "grna"
params.grna_group_column_name = "grna_group"

// Mild command line argument processing; get the directory and name of the result file
File out_f = new File(params.result_fp)
result_file_name = out_f.getName()
result_dir = out_f.getAbsoluteFile().getParent()


// PROCESS 1: Check inputs; output the list of response IDs and grna groups
process check_inputs {
  debug true
  time "5m"
  memory "5 GB"

  input:
  path "multimodal_metadata_fp"
  path "response_odm_fp"
  path "grna_odm_fp"
  path "pairs_to_analyze_fp"

  /*
  output:
  path "response_to_pod_id_map.rds", emit: response_to_pod_id_map_ch
  path "pair_to_pod_id_map.rds", emit: pair_to_pod_id_map_ch
  path "response_pods.txt", emit: response_pods_ch
  path "pair_pods.txt", emit: pair_pods_ch
  path "mm_odm_new.rds", emit: multimodal_metadata_ch
  path "pairs_to_analyze_new.rds", emit: pairs_to_analyze_ch
  */

  """
  echo $multimodal_metadata_fp \
  $response_odm_fp \
  $grna_odm_fp \
  $pairs_to_analyze_fp \
  $params.response_modality_name \
  $params.grna_modality_name \
  $params.grna_group_column_name \
  $params.response_pod_size \
  $params.pair_pod_size \
  $params.n_pairs_to_sample
  """
}


workflow {
  check_inputs(Channel.fromPath(params.multimodal_metadata_fp, checkIfExists : true),
               Channel.fromPath(params.response_odm_fp, checkIfExists : true),
               Channel.fromPath(params.grna_odm_fp, checkIfExists : true),
               Channel.fromPath(params.pairs_to_analyze_fp, checkIfExists : true))
}
