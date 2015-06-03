#### Load C Library ####
path.string <- strsplit(getwd(),'/')[[1]]

if(path.string[length(path.string)-1]!='exchange'){
  warning('Warning: you must load simulation_R-ext.so manually unless you are in a directory one below exchange.\nExample: exchange/R, exchange/preston, exchange/rundir\n')
}else{ dyn.load('../C/eta_simulation/SE.simulation_R-ext.so')}

rm(path.string)

calc.moments.eta <- function(gene.list,trna.dat,A_1 = 4,A_2 = 4,B = 0.0025){
  #Purpose: Analytically calculate first and second moments of the distribution of eta values based
  #         on math by Russ. 
  #Inputs: gene.list - gene list with codon index located at gene.list$gene.dat$c_index
  #        trna.dat - data structure containing codon specific parameters 
  #        A.1 - Translation initiation cost
  #        A.2 - Translation elongation cost
  #        B - background nonsense error rate
  #Outputs: Named vector containing:
  #           [1] - eta.mean
  #           [2] - eta.var
  #           [3] - eta.min
  #           [4] - eta.max
  #           [5] - eta.obs
  
  ############################
  # SEQUENCE DATA/PARAMETERS #
  ############################
  
  CODON_INDEX = gene.list$gene.dat$c_index
  if(CODON_INDEX[length(CODON_INDEX)]==99||CODON_INDEX[length(CODON_INDEX)]==61||CODON_INDEX[length(CODON_INDEX)]==62||CODON_INDEX[length(CODON_INDEX)]==63){
    CODON_INDEX <- CODON_INDEX[-length(CODON_INDEX)] #remove the last codon(stop codon) -- stop codon is not used in c
  }
  AA_COUNT = length(CODON_INDEX)
  PHI = 0#gene.list$phi.value
  
  AA.VEC = trna.dat$aa[sort.list(trna.dat$c_index)]
  CODON.VEC = trna.dat$codon[sort.list(trna.dat$c_index)]
  
  
  #####################
  # GENOME PARAMETERS #
  #####################
  
  ELONGATION_RATES = trna.dat$elong_rate[sort.list(trna.dat$c_index)]
  MUTATION_RATES = trna.dat$mut_rate[sort.list(trna.dat$c_index)]
  AT_BIAS=0.5 #at_bias parameter (currently unused)
  GAMMA = 1 #transition/transversion parameter (currently unused)
  IGNORE = 0 #Integer number to ignore last # amino acids
  
  #####################################
  # Initialize output dummy variables #
  #####################################
  
  LIK = 0
  ETA_MEAN = 0
  ETA_VAR = 0
  ETA_MIN = 0
  ETA_MAX = 0
  ETA_OBS = 0
  
  ###################
  # Call C function #
  ###################
  
  OUTPUT <- .C("calc_moments", as.integer(CODON_INDEX),     AA_COUNT = length(CODON_INDEX),   as.double(PHI),
               as.double(ELONGATION_RATES), as.double(MUTATION_RATES),        as.double(1),
               as.double(A_1),              as.double(A_2),                   as.double(AT_BIAS),
               as.double(B),                as.integer(IGNORE),               as.double(GAMMA),
               as.double(1),                as.double(ETA_MEAN),              as.double(ETA_VAR),
               as.double(ETA_MIN),          as.double(ETA_MAX),               as.double(ETA_OBS),
               as.character(AA.VEC),        as.character(CODON.VEC),          as.integer(length(trna.dat$codon)))
  
  ##################
  # Extract Output #
  ##################
  
  eta_mean <- OUTPUT[[14]]
  eta_var <- OUTPUT[[15]]
  eta_min <- OUTPUT[[16]]
  eta_max <- OUTPUT[[17]]
  eta_obs <- OUTPUT[[18]]
  
  c(eta.mean=eta_mean,eta.var=eta_var,eta.min=eta_min,eta.max=eta_max,eta.obs=eta_obs)
}