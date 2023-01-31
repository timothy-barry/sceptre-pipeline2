get_response_ids_from_pod <- function(response_to_pod_id_map, curr_pod) {
  curr_responses <- response_to_pod_id_map |>
    dplyr::filter(pod_id == curr_pod) |>
    dplyr::pull(response_id) |>
    as.character()
}


#' Run response precomputation
#'
#' Runs the precomputation for a response. This consists of regressing the response expressions onto the covariate matrix.
#'
#' @param expressions the numeric vector of response expressions
#' @param covariate_matrix the covariate matrix on which to regress the expressions (NOTE: the matrix should contain an interecept term)
#' @param use_mass should the MASS packge be used to carry out the NB regression? If FALSE, we instead carry out an initial Poisson regression, estimate theta by using the residuals, and then carry out a final NB regression, treating theta (as estimated in the previous step) as known. Setting `use_mass` to FALSE is faster, but the resulting estimates are expected to be slightly less accurate.
#'  
#' @return a list containing the following elements:
#' (i) "precomp_str": a string summarizing the method used to fit the GLM and calculate the NB size parameter. The part of the string before the colon is either "nb" or "pois," indicating whether NB regression or Poisson regression was used. ("pois_{warn}" indicates that Poisson regression was used and that a warning was generated). The part of the string after the colon indicates the method used to estimate the size parameter ("mass" for the "mass" package, "resid_mle" for MLE on the residuals, and "resid_mm" for method of moments on the residuals.)  
#' (ii) "fitted_coefs": a vector of fitted coefficients; the final entry of this vector is the fitted theta. 
run_response_precomputation_low_level <- function(expressions, covariate_matrix, use_mass) {
  # backup: return fitted coefficients from Poisson regression
  backup_3 <- function(pois_fit, pois_warn) {
    list(fitted_coef_str = paste0("pois", if (pois_warn) "_(warn)" else NULL),
         fitted_coefs = stats::coef(pois_fit))
  }
  
  # Backup: method of moments
  backup_2 <- function(pois_fit) {
    list(theta_fit_str = "resid_mm",
         MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual))
  }
  
  # Backup: MLE on poisson reg
  backup <- function() {
    # fit Poisson model, tracking warning
    pois_warn <- FALSE
    wHandler <- function(w) {pois_warn <<- TRUE; invokeRestart("muffleWarning")}
    withCallingHandlers(expr = {
      pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    }, warning = wHandler)
    
    # get theta; save theta itself (in response_theta) and theta_fit_string (i.e., procedure used to fit theta)
    response_theta_list <- tryCatch({
      list(theta_fit_str = "resid_mle",
           response_theta = MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1])
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
    theta_fit_str <- response_theta_list$theta_fit_str
    response_theta <- response_theta_list$response_theta
    response_theta <- min(max(response_theta, 0.1), 1000)
    
    # obtain the fitted coefficients
    fitted_coefs_list <- tryCatch({
      fit_nb <- stats::glm(formula = expressions ~ . + 0,
                           family = MASS::negative.binomial(response_theta), data = covariate_matrix)
      list(fitted_coef_str = "nb",
           fitted_coefs = stats::coef(fit_nb))
    }, error = function(e) backup_3(pois_fit, pois_warn), warning = function(w) backup_3(pois_fit, pois_warn))
    fitted_coefs <- fitted_coefs_list$fitted_coefs
    precomp <- c(fitted_coefs, response_theta = response_theta)
    fitted_coef_str <- fitted_coefs_list$fitted_coef_str
    
    # construct the precomp string
    precomp_str <- paste0(fitted_coef_str, ":", theta_fit_str)
    list(precomp_str = precomp_str,
         precomp = c(fitted_coefs, response_theta = response_theta))
  }
  
  if (use_mass) {
    # try to fit a negative binomial GLM with unknown dispersion
    result <- tryCatch({
      fit_nb <- MASS::glm.nb(formula = expressions ~ . + 0, data = covariate_matrix)
      response_theta <- min(max(fit_nb$theta, 0.1), 1000)
      fitted_coefs <- stats::coef(fit_nb)
      list(precomp_str = "nb:mass",
           precomp = c(fitted_coefs, response_theta = response_theta))
    }, error = function(e) backup(), warning = function(w) backup()) 
  } else {
    result <- backup()
  }
  
  return(result)
}


run_response_precomputation_high_level <- function(response_modality, covariate_matrix, curr_responses, nt_idxs, use_mass = TRUE) {
  precomp <- lapply(curr_responses, function(curr_response) {
    # print(paste0("Regressing response ", curr_response, " onto the covariates."))
    # load the expression data
    expressions <- as.numeric(response_modality[[curr_response,]])
    if (low_moi) expressions <- expressions[nt_idxs] # subset expressions
    curr_precomp <- run_response_precomputation_low_level(expressions, covariate_matrix, use_mass)
  })
  
  precomp_coef_df <- sapply(X = precomp, FUN = function(l) l[["precomp"]]) |>
    t() |> data.table::as.data.table() |> dplyr::mutate(response_id = curr_responses)
  precomp_info_df <- data.frame(precomp_str = sapply(X = precomp, FUN = function(l) l[["precomp_str"]]),
                                response_id = curr_responses)
  
  return(list(precomp_coef_df = precomp_coef_df, precomp_info_df = precomp_info_df))
}
