# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# AutoCal-Eco: An R Framework for Automated Calibration of Ecospace
#
# This script contains the core function `run_calibration` which
# encapsulates the logic for running different optimization methods.
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#SETUP--------------------------------------------------------------------------
library('GA')
library('cmaes')
library('rBayesianOptimization')
library('dplyr')
library('digest')
library('readxl')
library('doParallel')
library('R.utils')

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

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#OBJECTIVE FUNCTION-------------------------------------------------------------
objective_function <- function(log_par_vec, run_path) {
  par_vec <- exp(log_par_vec)-1
  config_hash <- digest(par_vec, algo = "md5")
  
  if (exists(config_hash, envir = cache)) {
    return(cache[[config_hash]])
  }
  
  # this_dir <- tempfile(tmpdir = run_dir)
  # dir.create(this_dir)
  # #out_path <- paste0(this_dir, "-output")
  # out_path <- paste0(this_dir)
  # dir.create(out_path)
  this_dir <- run_path
  
  tags.vul = tags.env = character()
  if(do.vuls){
    vuln_vec = par_vec[vul.par.idx]
    tags.vul <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey_pairs$pred, " ", ifelse(is.na(predprey_pairs$prey),"",predprey_pairs$prey),
                       "), ", sprintf("%.5f", vuln_vec), ", Indexed.Single")
  }
  
  if(do.env){
    env_vec = par_vec[env.par.idx]
    tags.env <- paste0("<ECOSIM_ENVIRONMENTAL_RESPONSE_INDEXED>(", respfxn_num,")",sprintf("%.5f", env_pars), ", Indexed.Single[]")
  }
  
  tags = c(tags.vul,tags.env)
  
  cmd_j <- cmd_base
  cmd_j[startsWith(cmd_j, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", this_dir)
  cmd_j <- c(cmd_j, tags)
  cmd_file <- file.path(this_dir, "cmd.txt")
  writeLines(cmd_j, cmd_file)
  
  
  result <- tryCatch({
    withTimeout({
      fn.runEwE(dir.cmdfile=cmd_file, do.obj=calibration)
    }, timeout = 3600, onTimeout = "error", silent = TRUE)
  }, error = function(e) {
    cat("Error in fn.runEwE:", conditionMessage(e), "\n")
    return(Inf)
  })
  
  score <- if (inherits(result, "try-error") || !is.finite(result[1])) Inf else result[1]
  
  if (!is.finite(score)) {
    cat("Model failed for parameters:", paste(round(par_vec, 3), collapse = ", "), "\n")
  }
  
  
  cat(sprintf("Run complete. Score = %.4f\n", score))
  cache[[config_hash]] <- score
  unlink(this_dir, recursive=TRUE)
  return(score)
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# CALIBRATION FUNCTION----------------------------------------------------------
fn.calibrate_ecospace<- function(method, do.vuls, do.env, predprey_pairs, env_respfxns, calibration = 1,
                                  n_cores = max(detectCores()-1,1), config = list(), output_base=dir.out) {
   
  # calibration = 1
  # n_cores = max(detectCores()-1,1)
  # method = "GA"
  # do.vuls = TRUE
  # do.env = FALSE
  # predprey_pairs = predcol_vuls #predprey.sens[,c(1:2,5)]
  # env_respfxns = NULL
  # output_base = dir.out
  # config <- myconfig #list(popSize = 24, run = 10, pmutation = 0.2, maxiter = 50)
  
  seed <- 42
  set.seed(seed)
  
  ## validate setup-------------------------------------------------------------
  if (!exists("cmd_base") || !exists("fn.runEwE")) {
    stop("Required variables 'cmd_base' and/or 'fn.runEwE' not found. Ensure setup.R is sourced.")
  }
  
  if(do.vuls){
    pdcol = unique(predprey_pairs$pred[is.na(predprey_pairs$prey)])
    pdpy = c(0,unique(predprey_pairs$pred[!is.na(predprey_pairs$prey)]))
    if(length(which(pdcol %in% pdpy))>0){
      stop("Trying to estimate ki and kij for at least on predator. Check predprey_pairs input.")
    }
  }
  
  ## set configurations----------------------------------------------------------
  # These can be overridden by the user-provided 'config' list.
  defaults <- list(
    GA = list(popSize = 25, run = 10, pmutation = 0.2, maxiter = 2000, elitism=2),
    CMAES = list(stop.if.no.improvement = 250, sigma = NULL, maxit = Inf),
    BO = list(init_points = 50, n_iter_chunk = 10, stop.if.no.improvement = 250)
  )

  # Merge user config with defaults
  run_config <- defaults[[method]]
  run_config[names(config)] <- config
  
  ## make parameter vector---------------------------------------------------
  n_vuls <- n_env <- 0
  if(do.vuls) n_vuls <- nrow(predprey_pairs)
  if(do.env) n_env <- nrow(env_respfxn)  #need to melt the env resp fxn dataframe
  #medations functions
  
  n_pars = n_vuls + n_env
  if (n_pars == 0) stop("n_pars==0: No predator-prey pairs or environmental responses parameters to estimate.  Check inputs.")
  
  log_vuln_vec = log_env_vec = numeric()
  if(do.vuls) log_vuln_vec = log(predprey_pairs$baseval+1)
  if(do.env) log_env_vec = log(env_respfxns$baseval+1)
  log_par_vec = c(log_vuln_vec,log_env_vec)
  
  # index parameter types-------------------------------------------------------
  #need to think about how to index the parameter vector as it continues to grow
  vul.par.idx = env.par.idx = numeric()
  if(do.vuls & !do.env){ 
    vul.par.idx = 1:n_vuls
  } else if(do.env & !do.vuls) {
    env.par.idx = 1:n_env
  } else {
    vul.par.idx = 1:n_vuls
    env.par.idx = (n_vuls+1):(n_vuls+n_env)
  }  
  length(log_par_vec)
  log_par_vec
  
  # set parameter bounds--------------------------------------------------------
  VULN_MIN <- 1.01
  VULN_MAX <- 1e6
  
  lower.vuls <- upper.vuls <- lower.env <- upper.env <- numeric()
  lower.vuls <- rep(log(VULN_MIN+1),n_vuls)
  upper.vuls <- rep(log(VULN_MAX+1),n_vuls)
  lower.env <- rep(0,n_env)
  upper.env <- rep(1,n_env)
  
  # create output directory------------------------------------------------------
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(output_base, paste0(method, "_Run_", timestamp))
  dir.create(run_dir, recursive = TRUE)
  cat(sprintf("Output will be saved in: %s\n", run_dir))
  
  # Initialize cache------------------------------------------------------------
  cache <- new.env(hash = TRUE, parent = emptyenv())
  
  # GA SETUP----------------------------------------------------------------------
  #=== 3. Standardized Objective Function (Minimization) ===
  # This function is used by all methods. It returns the raw score (lower is better).
  ###initialize GA----
  ga.popfxn = function(object){
    init.pars = log_par_vec #log(predprey_pairs$baseval)
    popSize = run_config$popSize
    nBits = n_pars
    matrix(rep(init.pars,popSize),ncol=nBits, nrow=popSize, byrow = T)
  }
  
  ### fitness wrapper-------------------------------------------------------------
  fitness_wrapper <- function(log_par_vec) {
    tryCatch({
    #create unique run directory for each thread
    run_id <- paste0("run_", digest(log_par_vec, algo=c("md5")))
    run_path <- file.path(run_dir, run_id)
    dir.create(run_path, showWarnings = TRUE)
    
    
    # Write input files, run model, read output
    score <- objective_function(log_par_vec, run_path=run_path )
    
    # Clean up if needed
    unlink(run_path, recursive = TRUE)
    
    return(-score)
    }, error = function(e) {
      cat("Error in fitness_wrapper:", conditionMessage(e), "\n")
      return(-Inf)
    })
  }


  ### setup clusters----   
  closeAllConnections()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl,
                 c("objective_function", "cache", "predprey_pairs", "run_dir","do.vuls","do.env","env_respfxns","vul.par.idx","env.par.idx","log_par_vec",
                   "cmd_base", "fn.runEwE", "digest", "withTimeout", "calibration", "lower.vuls","upper.vuls","lower.env","upper.env",
                   "ga.popfxn","run_config","seed", "fitness_wrapper","file.console","obs.ts","group.names","df.names","startyear","endyear_sens"),
                 envir = environment())
  
  clusterEvalQ(cl, {
    library(digest)
    library(R.utils)
    library('GA')
    library('cmaes')
    library('rBayesianOptimization')
    library('dplyr')
    library('digest')
    library('readxl')
    library('doParallel')
    library('R.utils')
    # Add any other packages used in fn.runEwE or fn.objfxn1
  })
  
  
  
  # # === 4. Method-Specific Execution ===
  # result_object <- NULL
  # best_par_log <- NULL
  
  ### run ga--------------------------------------------------------------------
  cat("Starting Genetic Algorithm optimization...\n")
  result_object <- ga(
      type = "real-valued",
      fitness = fitness_wrapper,
      lower = c(lower.vuls, lower.env), #rep(log(VULN_MIN), n_vars),
      upper = c(upper.vuls, upper.env), #rep(log(VULN_MAX), n_vars),
      popSize = run_config$popSize,
      run = run_config$run,
      maxiter = run_config$maxiter,
      pmutation = run_config$pmutation,
      population = ga.popfxn,
      seed = seed,
      elitism = run_config$elitism,
      parallel = cl,
      monitor = function(obj) cat(sprintf("Generation %d: Best fitness = %.4f\n", obj@iter, obj@fitnessValue))
    )
    
    best_par_log <- result_object@solution[1, ]
    summary(result_object)
    result_object@fitness

  # Stop parallel cluster if it was started
  if (exists("cl")) stopCluster(cl)
  closeAllConnections()  

  # Save Final Results----------------------------------------------------------
  cat("\nOptimization complete. Saving results...\n")
  best_pars <- exp(best_par_log)-1
  final_cmd <- cmd_base
  final_tags <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey_pairs$pred, " ", predprey_pairs$prey,
                       "), ", sprintf("%.5f", best_pars), ", Indexed.Single")
  final_cmd <- c(final_cmd, final_tags)
  final_output_dir <- file.path(run_dir, "final_output")
  final_cmd[startsWith(final_cmd, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", final_output_dir)
  writeLines(final_cmd, file.path(run_dir, "final_cmd.txt"))
  
  results_df <- bind_cols(predprey_pairs, vuln = best_pars)
  write.csv(results_df, file.path(run_dir, "optimized_vulnerabilities.csv"), row.names = FALSE)
  
  saveRDS(result_object, file.path(run_dir, "result_object.rds"))
  cat("Results saved successfully.\n")
  
  return(list(best_parameters = results_df, result_object = result_object))
}
