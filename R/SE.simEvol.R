#### Load C Library ####
path.string <- strsplit(getwd(),'/')[[1]]

if(path.string[length(path.string)-1]!='exchange'){
  warning('Warning: you must load simulation_R-ext.so manually unless you are in a directory one below exchange.\nExample: exchange/R, exchange/preston, exchange/rundir\n')
}else{ dyn.load('../C/eta_simulation/SE.simulation_R-ext.so')}

rm(path.string)

simEvol <- function(phi,c_index,trna.dat, BIS, GMT=0, MES=0,SIMULATION.METHOD = 'M',
                    A.1 = 4, A.2 = 4, B = 0.0025, Ne = 1, Q = 1, IGNORE=FALSE,
                    AT.BIAS=0.5,GAMMA=1,BENCH=0,RANDOM=0){
  #Purpose: Simulate codon sequence based on a set of parameters. 
  #Inputs: phi - gene expression rate
  #        c_index - codon index 
  #        trna.dat - trna.dat structure. Can be created with SE.semppr.data.load
  #        BIS - Number of simulation steps per amino acid
  #        GMT - Amount of simulation time (only use with Evolution simulation method)
  #        MES - Number of simulation steps (you can use this or BIS or both)
  #        SIMULATION.METHOD - Choose between MCMC or Evolutionary simulation to simulate codon sequence.
  #                            MCMC is ~10x faster.
  #        RANDOM - Option to start the codon sequence at a random sequence
  #        A.1 - Translation initiation cost
  #        A.2 - Translation elongation cost
  #        B - background nonsense error rate
  #        Ne - effective population size
  #        Q - scaling factor
  #        IGNORE - number of amino acids to ignore at tail end of c_index. some theorize that 
  #                 proteins can function properly missing the last few AA's in sequence. For 
  #                 this reason, NSE may not cause strong CUB at tail end of codon sequence.
  #        AT.BIAS - DEPRECATED, codon specific mutation rates replaces AT.BIAS model
  #        GAMMA - DEPRECATED, codon specific mutation rates replaces GAMMA transition/transversion model
  #        BENCH - DEPRECATED
  #Outputs: c_index - simulated codon index
  
  
  
  #Debug
  #SIMULATION.METHOD = 'M';A.1 = 4; A.2 = 4; B = 0.0025; Ne = 1; Q = 1; AT.BIAS=0.5;BENCH=0;MES=0;IGNORE=0;RANDOM=0;GAMMA=1;N
  
  
  
  #######################
  # SEQUENCE PARAMETERS #
  #######################
  
  CODON.INDEX = c_index  
  if(CODON.INDEX[length(CODON.INDEX)]==99||CODON.INDEX[length(CODON.INDEX)]==61||CODON.INDEX[length(CODON.INDEX)]==62||CODON.INDEX[length(CODON.INDEX)]==63){
    CODON.INDEX <- CODON.INDEX[-length(CODON.INDEX)] #remove the last codon(stop codon) -- stop codon is not used in c
  }
  
  AA.COUNT = length(CODON.INDEX)
  PHI = phi
  
  #####################
  # GENOME PARAMETERS #
  #####################
  
  ELONGATION.RATES = trna.dat$elong_rate[sort.list(trna.dat$c_index)]
  MUTATION.RATES = trna.dat$mut_rate[sort.list(trna.dat$c_index)]
  AA.VEC = trna.dat$aa[sort.list(trna.dat$c_index)]
  CODON.VEC = trna.dat$codon[sort.list(trna.dat$c_index)]
  
  #Call C function
  OUT.DATA <- .C("CES_new",c_index = as.integer(CODON.INDEX),         as.double(ELONGATION.RATES),      as.double(MUTATION.RATES),
                           AA.COUNT = length(CODON.INDEX),            as.double(PHI),                   as.double(Ne), 
                           as.double(A.1),                            as.double(A.2),                   as.double(AT.BIAS),
                           as.integer(BENCH),                         as.double(B),                     as.double(GMT),
                           as.integer(IGNORE),                        as.double(GAMMA),                 as.double(Q),
                           as.integer(MES),                           as.integer(BIS),                  as.character(AA.VEC), 
                           as.character(CODON.VEC),          as.integer(length(trna.dat$codon)),        as.character(SIMULATION.METHOD),
                           PACKAGE='SE.simulation_R-ext')
    
  
  #return 
  OUT.DATA$c_index
  
}
