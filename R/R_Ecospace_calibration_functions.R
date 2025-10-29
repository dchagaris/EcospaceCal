# ===================================================================
# AutoCal-Eco: An R Framework for Automated Calibration of Ecospace
#
# This script contains the core function `run_calibration` which
# encapsulates the logic for running different optimization methods.
# ===================================================================

# === Load All Required Libraries ===
library('GA')
library('cmaes')
library('rBayesianOptimization')
library('dplyr')
library('digest')
library('readxl')
library('doParallel')
library('R.utils')

# ===================================================================
# Main Calibration Function
# ===================================================================
#' Run an automated calibration for an Ecospace model.
#'
#' @param method The optimization method to use. One of "GA", "CMAES", or "BO".
#' @param sensitivity_file Path to the Master Vulnerability Table Excel file.
#' @param sheet_number The sheet number to use in the Excel file.
#' @param output_base The base directory name for saving output files.
#' @param calibration The calibration method to be used in the objective function. 1 = Temporal (default), 2 = Spatiotemporal
#' @param n_cores The number of cores to use. There is a default of one less than the maximum available and a minimum of 1.
#' @param config A list of configuration options specific to the chosen method. No response corresponds to the default.
#'
#' @return A list containing the best parameters and the full result object.
#' 
  
run_calibration <- function(method, predprey_pairs, output_base,calibration = 1,n_cores = max(detectCores()-1,1), config = list()) {
# 
  # calibration = 1
  # n_cores = max(detectCores()-1,1)
  # method = "GA"
  # predprey_pairs = predprey.sens[,1:2]
  # output_base = dir.out
  # config <- list(popSize = 200, run = 200, pmutation = 0.2, maxiter = 2000)
  # 
  # === 1. Default Configurations ===
  
  # These can be overridden by the user-provided 'config' list.
  defaults <- list(
    GA = list(popSize = 25, run = 10, pmutation = 0.2, maxiter = 2000),
    CMAES = list(stop.if.no.improvement = 250, sigma = NULL, maxit = Inf),
    BO = list(init_points = 50, n_iter_chunk = 10, stop.if.no.improvement = 250)
  )

  
  # Merge user config with defaults
  run_config <- defaults[[method]]
  run_config[names(config)] <- config
  
  # === 2. General Setup ===
  VULN_MIN <- 1.01
  VULN_MAX <- 1e6
  seed <- 42
  set.seed(seed)

  # Validate setup from sourced file
  if (!exists("cmd_base") || !exists("fn.runEwE")) {
    stop("Required variables 'cmd_base' and/or 'fn.runEwE' not found. Ensure setup.R is sourced.")
  }
  
  n_vars <- nrow(predprey_pairs)
  if (n_vars == 0) stop("No predator-prey pairs found")
  
  # Setup output directory
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(output_base, paste0(method, "_Run_", timestamp))
  dir.create(run_dir, recursive = TRUE)
  cat(sprintf("Output will be saved in: %s\n", run_dir))
  
  # Initialize cache
  cache <- new.env(hash = TRUE, parent = emptyenv())
  
  # === 3. Standardized Objective Function (Minimization) ===
  # This function is used by all methods. It returns the raw score (lower is better).
  objective_function <- function(log_vuln_vec) {
    #log_vuln_vec = log(predprey.sens$baseval+1)
    vuln_vec <- exp(log_vuln_vec)
    config_hash <- digest(vuln_vec, algo = "md5")
    
    
    if (exists(config_hash, envir = cache)) {
      return(cache[[config_hash]])
    }
    
    tags <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey_pairs$pred, " ", predprey_pairs$prey,
                   "), ", sprintf("%.5f", vuln_vec), ", Indexed.Single")
    this_dir <- tempfile(tmpdir = run_dir)
    dir.create(this_dir)
    #out_path <- paste0(this_dir, "-output")
    out_path <- paste0(this_dir)
    dir.create(out_path)
    cmd_j <- cmd_base
    cmd_j[startsWith(cmd_j, "<ECOSPACE_OUTPUT_DIR>")] <-
      sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", out_path)
    cmd_j <- c(cmd_j, tags)
    cmd_file <- file.path(this_dir, "cmd.txt")
    writeLines(cmd_j, cmd_file)
    
    # Temporal Calibration
    if (calibration == 1){
    
      result <- try(withTimeout(
        fn.runEwE(dir.cmdfile = cmd_file, do.obj = 1),
        timeout = 3600, onTimeout = "error", silent = TRUE
      ))
      
    # Spatio-temporal Calibration
    } else {
      result <- try(withTimeout(
        fn.runEwE(dir.cmdfile = cmd_file, do.obj = 2),
        timeout = 3600, onTimeout = "error", silent = TRUE
      ))
    }
    
    if (inherits(result, "try-error")) {
      score <- Inf
    } else {
      score <- result[1]
    }
    cat(sprintf("Run complete. Score = %.4f\n", score))
    cache[[config_hash]] <- score
    return(score)
  }
  
  
  # === 4. Method-Specific Execution ===
  result_object <- NULL
  best_vuln_log <- NULL
  

  # Setup parallel backend for GA and CMAES
  closeAllConnections()
  #if (method %in% c("GA", "CMAES")) {
  #   cl <- makeCluster(n_cores)
  #   registerDoParallel(cl)
  #   clusterExport(cl,
  #                 c("objective_function", "cache", "predprey_pairs", "run_dir",
  #                   "cmd_base", "fn.runEwE", "digest", "withTimeout", "calibration"),
  #                 envir = environment())
  # }
  
  # --- Genetic Algorithm ---
  if (method == "GA") {
    
    
    #initialize with base values
    ga.popfxn = function(object){
      init.vuln = log(predprey.sens$baseval[1:n_sens_kij])
      popSize = run_config$popSize
      nBits = n_sens_kij
      matrix(rep(init.vuln,popSize),ncol=nBits, nrow=popSize, byrow = T)
    }
    
    cat("Starting Genetic Algorithm optimization...\n")
    # GA maximizes, so we need a wrapper for our minimizing function
    fitness_wrapper <- function(log_vuln_vec) - objective_function(log_vuln_vec)
    
    result_object <- ga(
      type = "real-valued",
      fitness = fitness_wrapper,
      lower = rep(log(VULN_MIN), n_vars),
      upper = rep(log(VULN_MAX), n_vars),
      popSize = run_config$popSize,
      run = run_config$run,
      maxiter = run_config$maxiter,
      pmutation = run_config$pmutation,
      population = ga.popfxn,
      seed = seed,
      elitism = 10,
      parallel = TRUE,
      monitor = function(obj) cat(sprintf("Generation %d: Best fitness = %.4f\n", obj@iter, obj@fitnessValue))
    )
    best_vuln_log <- result_object@solution[1, ]
    summary(result_object)
    result_object@fitness
    # --- CMA-ES ---
  } else if (method == "CMAES") {
    cat("Starting CMA-ES optimization...\n")
    initial_params <- rep((log(VULN_MIN) + log(VULN_MAX)) / 2, n_vars)
    if(is.null(run_config$sigma)) run_config$sigma <- (log(VULN_MAX) - log(VULN_MIN)) / 4
    
    control_list <- list(
      sigma = run_config$sigma,
      maxit = run_config$maxit,
      stop.if.no.improvement = run_config$stop.if.no.improvement,
      parallel = TRUE,
      diag.verb = 10
    )
    result_object <- cma_es(
      par = initial_params,
      fn = objective_function,
      lower = rep(log(VULN_MIN), n_vars),
      upper = rep(log(VULN_MAX), n_vars),
      control = control_list
    )
    best_vuln_log <- result_object$par
    
    # --- Bayesian Optimization ---
  } else if (method == "BO") {
    cat("Starting Bayesian Optimization with custom early stopping...\n")
    # BO maximizes, so we need a wrapper
    bayes_opt_wrapper <- function(...) {
      score <- -objective_function(c(...))
      return(list(Score = score, Pred = 0))
    }
    bounds <- lapply(1:n_vars, function(i) c(log(VULN_MIN), log(100))) # Constrain for stability
    names(bounds) <- paste0("v", 1:n_vars)
    
    result_object <- BayesianOptimization(
      FUN = bayes_opt_wrapper, bounds = bounds, init_points = run_config$init_points,
      n_iter = 0, acq = "ucb", verbose = TRUE
    )
    
    evals_since_last_improvement <- 0
    while (evals_since_last_improvement < run_config$stop.if.no.improvement) {
      cat(sprintf("\n--- Running BO chunk. Evals since improvement: %d ---\n", evals_since_last_improvement))
      last_best_score <- result_object$Best_Value
      result_object <- BayesianOptimization(
        FUN = bayes_opt_wrapper, bounds = bounds, init_grid_dt = result_object$History,
        init_points = 0, n_iter = run_config$n_iter_chunk, acq = "ucb", verbose = TRUE
      )
      if (result_object$Best_Value > last_best_score) {
        evals_since_last_improvement <- 0
      } else {
        evals_since_last_improvement <- evals_since_last_improvement + run_config$n_iter_chunk
      }
    }
    best_vuln_log <- result_object$Best_Par
  } else {
    stop("Invalid method specified. Choose 'GA', 'CMAES', or 'BO'.")
  }
  
  # Stop parallel cluster if it was started
  if (exists("cl")) stopCluster(cl)
  
  # === 5. Save Final Results ===
  cat("\nOptimization complete. Saving results...\n")
  best_vuln <- exp(best_vuln_log)
  final_cmd <- cmd_base
  final_tags <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey_pairs$pred, " ", predprey_pairs$prey,
                       "), ", sprintf("%.5f", best_vuln), ", Indexed.Single")
  final_cmd <- c(final_cmd, final_tags)
  final_output_dir <- file.path(run_dir, "final_output")
  final_cmd[startsWith(final_cmd, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", final_output_dir)
  writeLines(final_cmd, file.path(run_dir, "final_cmd.txt"))
  
  results_df <- bind_cols(predprey_pairs, vuln = best_vuln)
  write.csv(results_df, file.path(run_dir, "optimized_vulnerabilities.csv"), row.names = FALSE)
  
  saveRDS(result_object, file.path(run_dir, "result_object.rds"))
  cat("Results saved successfully.\n")
  
  return(list(best_parameters = results_df, result_object = result_object))
}
