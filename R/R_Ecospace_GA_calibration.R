
#prepare parameter vector----
fn.makeparvec <- function(
  do.vuls = TRUE, 
  do.env = FALSE, 
  predprey_pairs = predcol_vuls,
  env_respfxns = NULL,
  vul.min = 1.01,
  vul.max = 1e6){
  
  #validate setup----
  if(do.vuls){
    pdcol = unique(predprey_pairs$pred[is.na(predprey_pairs$prey)])
    pdpy = c(0,unique(predprey_pairs$pred[!is.na(predprey_pairs$prey)]))
    if(length(which(pdcol %in% pdpy))>0){
      stop("Trying to estimate ki and kij for at least on predator. Check predprey_pairs input.")
    }
  }
  
  #make parameter vector
  n_vuls <- n_env <- 0
  if(do.vuls) n_vuls <- nrow(predprey_pairs)
  if(do.env) n_env <- nrow(env_respfxn)  #need to melt the env resp fxn dataframe
  #medations functions
  
  n_pars <<- n_vuls + n_env
  if (n_pars == 0) stop("n_pars==0: No predator-prey pairs or environmental responses parameters to estimate.  Check inputs.")
  
  log_vuln_vec = log_env_vec = numeric()
  if(do.vuls) log_vuln_vec = log(predprey_pairs$baseval+1)
  if(do.env) log_env_vec = log(env_respfxns$baseval+1)
  log_par_vec <<- c(log_vuln_vec,log_env_vec)
  
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
  vul.par.idx <<- vul.par.idx
  env.par.idx <<- env.par.idx
  
  # index pred-prey vulnerabilities
  predprey_pairs <<- predprey_pairs
  # set parameter bounds--------------------------------------------------------
  lower.vuls <- upper.vuls <- lower.env <- upper.env <- numeric()
  lower.vuls <- rep(log(vul.min+1),n_vuls)
  upper.vuls <- rep(log(vul.max+1),n_vuls)
  lower.env <- rep(0,n_env)
  upper.env <- rep(1,n_env)
  L.bounds <<- c(lower.vuls,lower.env)
  U.bounds <<- c(upper.vuls, upper.env)
}

#make GA populations----
fn.GApop = function(object){
  run_config <- myconfig
  matrix(runif(run_config$popSize*n_pars, L.bounds, U.bounds), nrow=run_config$popSize, byrow=T)
  #matrix(rep(log_par_vec, run_config$popSize), nrow=run_config$popSize, byrow=T)
}

#write command files----
fn.parvec2cmd <- function(log_par_vec){
  par_vec <- exp(log_par_vec)-1
  
  #output directory----
  run_id <- paste0("run_", digest(log_par_vec, algo=c("md5")))
  run_path <- file.path(run_dir, run_id)
  dir.create(run_path, showWarnings = TRUE)

  #parameter tags----
  tags.vul = tags.env = character()
  if(length(vul.par.idx)>0){
    vuln_vec = par_vec[vul.par.idx]
    tags.vul <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey_pairs$pred, " ", ifelse(is.na(predprey_pairs$prey),"",predprey_pairs$prey),
                       "), ", sprintf("%.5f", vuln_vec), ", Indexed.Single")
  }
  
  if(length(env.par.idx)>0){
    env_vec = par_vec[env.par.idx]
    tags.env <- paste0("<ECOSIM_ENVIRONMENTAL_RESPONSE_INDEXED>(", respfxn_num,")",sprintf("%.5f", env_pars), ", Indexed.Single[]")
  }
  
  tags = c(tags.vul,tags.env)
  
  
  #command files----
  cmd_j <- cmd_base
  cmd_j[startsWith(cmd_j, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", run_path)
  cmd_j <- c(cmd_j, tags)
  cmd_file <- file.path(run_path, "cmd.txt")
  writeLines(cmd_j, cmd_file)
  
}

#run the population of models----
fn.runEwE.gapop <-  function(
    files.cmd, 
    obj.fxn=1, 
    cl.export = list("files.cmd", "obs.ts")
){
  #source(file.setup)

  clusterExport(cl,append(cl.export,list("file.console", "fn.runEwE", "fn.objfxn1", "fn.objfxn2","startyear","endyear_sens","group.names","df.names")))

  #runlist=runlist[1:20,]
  pbar <- winProgressBar("Running Ecospace GApop",label=paste0("Simulation 0 of ",length(files.cmd)),max=100)
  prog <- function(n) setWinProgressBar(pbar,(n/length(files.cmd)*100),label=paste("Simulation Run", n,"of", 
                                                                                   length(files.cmd),"Completed"))
  opts <- list(progress=prog)
  
  #print(paste('Running',length(files.cmd),'Ecospace simulations'))
  t1 <- Sys.time()
  runs <- foreach(i = 1:length(files.cmd), .errorhandling = 'pass', .options.snow = opts) %dopar% {
    fn.runEwE(dir.cmdfile=files.cmd[i], do.obj = obj.fxn)
  }
  close(pbar)
  #print(paste('Run time',round(as.numeric(Sys.time()-t1),2)))
  
  ##missing runs----
  filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
  erruns <- which(filecheck<=1)  
  
  while(length(erruns)>=1){
    #print(paste0('Redo missing runs: n=',length(erruns)))
    
    pbar <- winProgressBar("Running Ecospace GApop: Missing Runs",label=paste0("Simulation 0 of ",length(erruns)),max=100)
    prog <- function(n) setWinProgressBar(pbar,(n/length(erruns)*100),label=paste("Simulation Run", n,"of", length(erruns),"Completed"))
    opts <- list(progress=prog)
    
    runs.erruns <- foreach(i=1:length(erruns),.errorhandling='pass',.options.snow=opts) %dopar% {
      fn.runEwE(dir.cmdfile=files.cmd[erruns[i]], do.obj=obj.fxn)
    }
    close(pbar)
    
    for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])
    
    filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
    erruns <- which(filecheck<=1)  #if there are many missing runs, then need to do this in parallel
  }
  
  fitness <- sapply(runs, function(x) x[1])   #cbind(runlist, do.call(rbind, runs))
  #print('All runs completed')
  #stopCluster(cl);
  #closeAllConnections()
  unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
  return(-fitness)
}

# === Selection ===
select_parents <- function(gapop, fitness) {
  ranks <- rank(fitness)
  probs <- ranks / sum(ranks)
  #cbind(ranks,fitness,probs)
  selected <- gapop[sample(1:nrow(gapop), pop_size, replace = TRUE, prob = probs), ]
  return(selected)
}

# === Crossover ===
crossover <- function(parents) {
  offspring <- parents
  for (i in seq(1, pop_size - 1, by = 2)) {
    if (runif(1) < 0.8) {
      point <- sample(1:(n_pars - 1), 1)
      temp <- offspring[i, (point + 1):n_pars]
      offspring[i, (point + 1):n_pars] <- offspring[i + 1, (point + 1):n_pars]
      offspring[i + 1, (point + 1):n_pars] <- temp
    }
  }
  return(offspring)
}

# === Mutation ===
mutate <- function(population) {
  #population=gapop
  low = apply(population,2,min)
  upp = apply(population,2,max)
  for (i in 1:nrow(population)) {
    for (j in 1:n_pars) {
      if (runif(1) < mutation_rate) {
        population[i, j] <- runif(1, low[j], upp[j])
      }
    }
  }
  return(population)
}


#GA function----
fn.GA <- function(myconfig){
  pop_size <<- myconfig$popSize
  n_generations <<- myconfig$n_gen
  mutation_rate <<- myconfig$pmutation
  elitism <<- myconfig$elitism
  
  #initial population
  gapop.init <- fn.GApop()
  apply(gapop.init,1,function(x) fn.parvec2cmd(x)) 
  files.cmd <- list.files(path=run_dir,pattern="cmd.txt", full.names=T, recursive=T)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, cl.export = list("obs.ts"))
  gapop <- gapop.init
  
  for (gen in 1:n_generations) {
    #gen=1
    cat(sprintf("Generation %d: Best fitness = %.4f\n", gen, max(fitness)))
    
    # Elitism
    elite_idx <- order(fitness, decreasing = TRUE)[1:elitism]
    elite <- gapop[elite_idx, ]
    
    # Selection, Crossover, Mutation
    parents <- select_parents(gapop, fitness)
    offspring <- crossover(parents)
    offspring <- mutate(offspring)
    
    # Evaluate new population
    apply(offspring,1,function(x) fn.parvec2cmd(x)) 
    files.cmd <- list.files(path=run_dir,pattern="cmd.txt", full.names=T, recursive=T)
    new_fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, cl.export = list("obs.ts"))
    
    # Combine elite + offspring
    offspring.rank = rank(new_fitness)
    drop.idx = which(offspring.rank<=elitism)
    gapop <- rbind(elite, offspring[-drop.idx,])
    fitness <- c(fitness[elite_idx], new_fitness[-drop.idx])
  }
}
##parents----
##crossover----
##mutate----
#write command files