main.semppr <- function(obs.genome,trna.dat,pop.parms=list(Q = 1, Ne = 1, A1 = 4, A2 =4, B = 0.0025),out.prefix,min.num.steps=500,max.num.steps=3000,BIS=4,GMT=0,MES=0,proposal.type = 'RN',simulation.type = 'M',n.cores = 2,simulate.dataset=TRUE,phi.only=FALSE){
  # Purpose: Fit an inter-ribosomal model of nonsense errors in protein translation described in Gilchrist et al. 2009 to genome data using an evolutionary
  #          model described in Sella and Hirsh 2005 in order to estimate gene-specific expression rates and genome-wide parameters like
  #          codon mutation rates and nonsense error rates. The enormous state space created by this model renders the posterior distribution of these parameters 
  #          doubly-intractible. Inference requires the MCMC Exchange method described in Murray et al. 2012 to circumvent calculation of parameter-dependent 
  #          normalization constants in the parameter acceptance ratios. This method requires the simultaneous proposal of an auxiliary variable along with
  #          parameter proposal. The auxiliary variable manifests itself in this model as a simulated gene or simulated genome, depending on whether the 
  #          proposed parameter is gene-specific, such as expression rate, or genome-wide, such as codon mutation rate.
  # Inputs: obs.genome - list of lists containing gene-specific information
  #                                        gene.list - list containting 
  #                                             $phi.value - expression rate
  #                                             $name - gene ID
  #                                             $gene.dat, which contains
  #                                                 $codon - nucleotide char sequence
  #                                                 $aa - amino acid char sequence
  #                                                 $c_index - representation of codon sequence using integers
  #                                                 $mut_rate - vector of mutation rates corresponding to c_index
  #                                                 $elong_rat - vector of elong rates...
  #                                                 $count - amino acid/codon count
  #         pop.parms - list containing genome-wide parameters
  #                                         $Q - Arbitrary scaling factor, default = 1
  #                                         $Ne - Effective population size, default = 1
  #                                         $A1 - Cost of translation initiation (in ATP's), default = 4
  #                                         $A2 - Cost of peptide elongation (in ATP's), default = 4
  #                                         $B - background Nonsense Error Rate (in s^-1), default = 0.0025 
  #         trna.dat - dataframe with rows corresponding to codon data, including: 
  #                                         [,1] - AA - Character corresponding to translated amino acid 
  #                                         [,2] - codon - string containing NT sequence of codon 
  #                                         [,3] - c_index - integer denoting codon index 
  #                                         [,4] - mut_rate - mutation rate to codon 
  #                                         [,5] - elong_rate - elongation rate of codon
  #                                         [,6] - count - of synonyms for every codon
  #         out.prefix - desired folder to write output files
  #         min.num.steps - minimum number of MCMC steps for each gene
  #         max.num.steps - maximum number of MCMC steps for each gene
  #         BIS - For Aux variable generation: Number of simulation steps per amino acid
  #         GMT - For Aux variable generation: Amount of simulation time (only use with Evolution simulation method)
  #         MES - For Aux variable generation: Number of simulation steps (you can use this or BIS or both)
  #         proposal.type - Phi proposal type, options:
  #                                      'RN' - reflecting normal proposal centered around current phi value. adapts to 
  #                                             acceptance ratio >0.25 and <0.35
  #                                      'LN' - Lognormal proposal centered around current phi value. adapts to accept
  #                                             ratio >0.25 and <0.35
  #         simulation.type - For Aux variable generation: simulation method, options:
  #                                      'M' - MCMC
  #                                      'E' - Evolution simulation
  #         n.cores - number of cores used for mclapply
  #         simulate.dataset - Flag to simulate datset 
  #         phi.only - Flag to estimate phi only
  # Outputs: Files located in out.prefix
  # Usage: See 'R/test_script1.R' for a usage example
  
  #For Debugging
  #min.num.steps=500;max.num.steps=3000;BIS=10;GMT=0;MES=0;proposal.type = 'RN';simulation.type = 'M';n.cores = 2
  
  
  #load required libraries and files
  require(multicore)
  source('../R/SE.semppr.data.load.R')
  source('../R/SE.simEvol.R')
  source('../R/SE.gen.seq.R')
  source('../R/SE.calcEta.R')
  source('../R/SE.calcMu.R')
  source('../R/SE.gelman.conv.R')
  dyn.load('../C/eta_simulation/SE.simulation_R-ext.so')
  
  #create output directory if it doesn't exist
  dir.create(out.prefix)
  
  #make sure out.prefix string ends in a '/'
  paste(paste(strsplit(out.prefix,'/')[[1]],collapse='/'),'/',sep='')
  
  #begin timing 
  start.time = proc.time() 
  
  #Set some parameters and other stuff
  pop.parms = list(Q = 1, Ne = 1, A1 = 4, A2 =4, B = 0.0025) #population parameters
  genome.end=list() #When a gene converges on the correct phi distribution, take 
                    #the gene out of genome list and put it in this list
  
  for(i in 1:length(obs.genome)){
    obs.genome[[i]]$phi.obs = obs.genome[[i]]$phi.value  
  }
  
  trna.dat.obs=trna.dat
  
  ##########################################################
  # Step 1) Burn in aux sequence to initial lambda and phi #
  ##########################################################
  
  genome = mclapply(obs.genome, function(gene.list,trna.dat,simulation.type,BIS,GMT,MES){
     gene.list$aux_c_index = simEvol(phi=gene.list$phi.obs,c_index=gene.list$gene.dat$c_index,trna.dat=trna.dat,BIS=10*BIS,GMT=GMT,MES=MES,SIMULATION.METHOD=simulation.type);
     gene.list$gene.dat$codon=NULL        # These next few lines reduce
     gene.list$gene.dat$aa=NULL           # the size of the genome object
     gene.list$gene.dat$elong_rate=NULL;  
     gene.list$gene.dat$mut_rate=NULL
     gene.list$gene.dat$count=NULL
     gene.list$gene.dat$c_index=as.integer(gene.list$gene.dat$c_index)
     gene.list$aux_c_index = as.integer(gene.list$aux_c_index)
     gene.list
  },trna.dat,simulation.type,BIS,GMT,MES,mc.cores=n.cores,mc.preschedule=TRUE)
  num.genes = length(genome)
  
  ##########################################################################################
  # Step 2) USE SIMULATED SEQUENCE--SET THE OBSERVED SEQUENCE = INITIAL BURNED IN SEQUENCE #
  ##########################################################################################
  
  if(simulate.dataset){
    genome = lapply(genome, function(gene){gene$gene.dat$c_index[-length(gene$gene.dat$c_index)] =gene$aux_c_index;gene})#
  }
  
  ###################################
  # Step 3) Guess phi based on SCUO #
  ###################################
  
  phi.dat = rlnorm(length(genome),0,1.5)
  
  #set max phi = 400. If phi gets too high, 
  #numerical precision becomes a problem because fitness=exp(-qNe*phi*eta)
  
  for(i in 1:length(genome)){
    if(phi.dat[i]>400) phi.dat[i]=400
  }
  
  scuo = calcScuo(genome,trna.dat,2)
  scuo.order = order(scuo,decreasing=TRUE)
  
  phi.dat=phi.dat[order(phi.dat,decreasing=TRUE)]
  phi.dat[scuo.order] = phi.dat
  
  #Set the phi value in the genome data structure
  for(i in 1:length(genome)){
    genome[[i]]$phi.value=phi.dat[i]
  }
  
  #TEST: start elong rate of codon 5 at a value different from observed

  if(!phi.only){ 
    genome.accept.history=rep(0,50)
    genome.prop.var=0.2
    trna.dat[5,4]=1
    trna.dat[5,5]=2
    elr.trace=0
    mut.trace=0
  }

  #Set parameters for mcmc
  genome = lapply(genome, function(gene){gene$MCMC$phi.proposal.var = gene$phi.value;gene$MCMC$phi.accept.rate = 0;gene$MCMC$phi.accept.hist=rep(x=0,50);gene})
  
  
  after.burnin.time = proc.time()
  
  ###########################
  #Step 4) MCMC starts here #
  ###########################
  
  i=0
  while(i<max.num.steps&&length(genome)>0){
   i=i+1
    ######################################
    # Step 4.1) Update phi for each gene #
    ######################################
    
    genome = mclapply(genome,function(gene.list,trna.dat,pop.parms,simulation.type,BIS,MES,GMT){
      
        prior.ratio = 1 #assume flat prior
        mh.ratio = 1 #metropolis-hastings ratio -- assume symmetric jumping distribution
      
        #4.1a) Propose new phi 
          # Use random walk with lognormal distribution or reflecting normal. 
          # Note: reflecting normal seems to work better because when using a lognormal
          # distbn, the phi traces can get stuck in low-phi "wells" 
        
        if(proposal.type == 'LN'){ #Lognormal Distribution has symmetric q'/q
          P.log.phi = rnorm(1,mean=log(gene.list$phi.value),sd=gene.list$MCMC$phi.proposal.var)
          P.phi = exp(new.log.phi)
        }else if(proposal.type == 'RN'){ #Reflecting Normal also has symmetric q'/q
          P.phi = rnorm(1,mean=gene.list$phi.value,sd=gene.list$MCMC$phi.proposal.var)
          if(P.phi<0) P.phi=-P.phi #reflecting normal distribution
        }else{
          error(paste('Error in SempprExchange.main.R: Invalid phi proposal type',proposal.type))
        }
        
        #4.1b) Generate auxiliary variable (simulate codon sequence)
        #      Russ thinks we need to start at a random sequence to make
        #      sure our samples are "exact and independent" like the 
        #      exchange algorithm requires, but it takes a very long time to burn in a random sequence
        #      for a high expression gene.
          
        #gene.list$aux_c_index <- rand.seq(gene.list,trna.dat)
        gene.list$aux_c_index <- simEvol(phi=P.phi,c_index=gene.list$aux_c_index,trna.dat=trna.dat,BIS=BIS,GMT=GMT,MES=MES,SIMULATION.METHOD=simulation.type)
        
        #4.1c) Calculate acceptance rate
          #delta eta is defined as eta_obs - eta_aux
        eta_obs = calcEta.NSE(codon_index=gene.list$gene.dat$c_index,trna.dat,pop.parms)
        eta_aux = calcEta.NSE(gene.list$aux_c_index,trna.dat,pop.parms)
        delta.eta = eta_obs - eta_aux
        
        #NOTE: Mu cancels out of this acceptance ratio 
        accept.ratio = exp(pop.parms$Q*pop.parms$Ne*delta.eta*(gene.list$phi.value - P.phi))*
          prior.ratio*
          mh.ratio
        
        #4.1d) Draw random number and accept/reject
        
        rn = runif(1)
        #debug output
        #cat(paste('phi: ',gene.list$phi.value,' prop.phi: ',P.phi,' eta.obs: ',eta_obs,'eta.aux: ',eta_aux,' accept.rat: ',accept.ratio,'\n'))
        if(rn < accept.ratio){
          gene.list$phi.value = P.phi
          gene.list$MCMC$phi.accept.hist = c(gene.list$MCMC$phi.accept.hist[2:50],1)
        }else{
          gene.list$MCMC$phi.accept.hist = c(gene.list$MCMC$phi.accept.hist[2:50],0)
        }
        
        #4.1e) Adjust the proposal variance every 50 steps
        if(i%%50==0){
          gene.list$MCMC$phi.accept.rate[i/50] = sum(gene.list$MCMC$phi.accept.hist)/50
          
          if(gene.list$MCMC$phi.accept.rate[i/50] > 0.45){ #if accept rate is too big, make jumps bigger
            gene.list$MCMC$phi.proposal.var = gene.list$MCMC$phi.proposal.var*1.1
          }else if(gene.list$MCMC$phi.accept.rate[i/50] < 0.35){ #if accept rate is too small, make jumps smaller
            gene.list$MCMC$phi.proposal.var = gene.list$MCMC$phi.proposal.var/1.1
          }
          
        }
        
        gene.list
    },trna.dat,pop.parms,simulation.type,BIS,MES,GMT,mc.cores=n.cores,mc.preschedule=TRUE)
    #Note: mc.preschedule = TRUE may work better when the num.gene/n.cores ratio is large
    #                       
    #      mc.preschedule = FALSE may work better when there are only ~5 genes per core
	
    #4.2) Save phi values for this step
    if(i==1){
      phi.dat = unlist(lapply(genome,function(list) list$phi.value))
    }else{
      phi.dat = rbind(phi.dat,unlist(lapply(genome,function(list) list$phi.value)))
    }
   
    #4.3) If we are not estimating genome-wide parameters, we can stop  
    #     taking phi draws from genes whose phi traces have converged. 
    #     estimation of genome-wide parameters
    #     like elong rate or NSE pr requries all genes
    #     Test for convergence using gelman convergence criteria. If phi has converged for a gene, move the gene 
    #     from genome to genome.end so mclapply will skip gene. Also remove
    #     the gene's phi trace from phi.dat
   
    if(phi.only){
      cat(paste('i',i,'length',length(genome),'\n')) 
      if(i%%20==0&&i>=min.num.steps){
        #save(phi.dat,file=paste(out.prefix,'phi.dat.temp'))  
        j=1
        if(length(genome)>1){	
          while(j <=length(genome)){
            cat(paste('j',j,'length',length(genome),'col',ncol(phi.dat),'\n'))
            R.hat=gelman.conv(phi.dat[,j])$R.hat
            if(j>1&&!is.na(R.hat)&&R.hat<1.1){
              genome[[j]]$phi.trace = phi.dat[,j]
              phi.dat=phi.dat[,-j];
              genome.end[[length(genome.end)+1]]=genome[[j]]
              genome[[j]]=NULL          
            }else{
              j=j+1
            }
          }
        }else{
          R.hat=gelman.conv(phi.dat)$R.hat
          if(!is.na(R.hat)&&R.hat<1.1){
            genome[[1]]$phi.trace = phi.dat
            genome.end[[length(genome.end)+1]]=genome[[1]]
            genome[[j]]=NULL
          }
        }
      }
    }else{
      #Start by estimating elong rate and mutation rate of one 2-codon amino acid
    
      #Amino Acid 2 is cysteine, which has 2 codons
    
      #Propose new trna.dat
      #TODO: this is where we should propose a new elong rate and new mu 
      #       at the same time based on a multivariate normal distribution.
      #NOTE: I'm proposing a new elongation rate and new mu simultaneously
      #      from independent normal distributions. This is the same as 
      #      proposing from a bivariate normal with 0 correlation
      
      P.trna.dat=trna.dat
      P.trna.dat[5,4]=exp(rnorm(1,log(as.numeric(trna.dat[5,4])),genome.prop.var))
      P.trna.dat[5,5]=exp(rnorm(1,log(as.numeric(trna.dat[5,5])),genome.prop.var))

      #Simulate aux sequence and calculate accept rate for the gene
      
    
      accept.rat.vec = mclapply(genome,function(gene.list,trna.dat,P.trna.dat){
        
        prior.ratio = 1 #assume flat prior
        mh.ratio = 1 #assume jumping distributions are symmetric--metropolis-hastings ratio=1      
        
        gene.list$aux_c_index <- simEvol(phi=gene.list$phi.value,c_index=gene.list$aux_c_index,trna.dat=P.trna.dat,BIS=BIS,GMT=GMT,MES=MES,SIMULATION.METHOD=simulation.type)
        eta_obs = calcEta.NSE(gene.list$gene.dat$c_index,trna.dat,pop.parms)
        eta_obs_p = calcEta.NSE(gene.list$gene.dat$c_index,P.trna.dat,pop.parms)
        eta_sim  = calcEta.NSE(gene.list$aux_c_index,trna.dat,pop.parms)
        eta_sim_p = calcEta.NSE(gene.list$aux_c_index,P.trna.dat,pop.parms)
        
        mu_obs = calcMu(gene.list$gene.dat$c_index,trna.dat)
        mu_obs_p = calcMu(gene.list$gene.dat$c_index,P.trna.dat)
        mu_sim = calcMu(gene.list$aux_c_index,trna.dat)
        mu_sim_p= calcMu(gene.list$aux_c_index,P.trna.dat)
        
        mu_obs*mu_sim_p/(mu_obs_p*mu_sim)*
          exp(-pop.parms$Q*pop.parms$Ne*gene.list$phi.value*(eta_obs_p+eta_sim-eta_obs-eta_sim_p))*
          prior.ratio*
          mh.ratio
          
          
      },trna.dat=trna.dat,P.trna.dat=P.trna.dat,mc.cores=n.cores,mc.preschedule=TRUE)
      
      #Since the entire genome shares these parameters, we must accept for all 
      #genes or reject for all genes. The product of acceptance ratios for
      #each gene gives the acceptance ratio for all genes since the acceptance
      #ratios are independent of other genes.
      accept.ratio=prod(unlist(as.numeric(accept.rat.vec)))
      cat(paste("i=",i,'elr:',trna.dat[5,4],'mut:',trna.dat[5,5],'p.elr:',P.trna.dat[5,4],'p.mut:',P.trna.dat[5,5],'accept:',accept.ratio,'accept.rate:',mean(genome.accept.history),'\n'))   
      if(runif(1)<accept.ratio){
        trna.dat=P.trna.dat
        genome.accept.history[1:49]=genome.accept.history[2:50]
        genome.accept.history[50]=1
      }else{
        genome.accept.history[1:49]=genome.accept.history[2:50]
        genome.accept.history[50]=0
      }
      elr.trace[i] = trna.dat[5,4]
      mut.trace[i] = trna.dat[5,5]
      
      #Adjust proposal variance every 50 steps
      if(i%%10==0){
        genome.accept.rate=mean(genome.accept.history)
        if(genome.accept.rate>0.45){
          genome.prop.var = genome.prop.var*1.1
        }else if(genome.accept.rate < 0.35){
          genome.prop.var = genome.prop.var/1.1
        }
        cat(paste('genome.prop.var: ',genome.prop.var,'\n'))
      }
    } #end if(!phi.only){}else{}
  }
  
  if(length(genome)>0){
    for(i in 1:length(genome)){
      genome[[i]]$phi.trace=phi.dat[,i]
      genome.end[[length(genome.end)+1]]=genome[[i]]
    }
    genome=genome.end
  }else{genome=genome.end}
  
  phi.dat = lapply(genome,function(list) list$phi.trace)
  phi.obs = unlist(lapply(genome,function(list)list$phi.obs))
  phi.ID = unlist(lapply(genome,function(list)list$name))
  
  phi.obs.dat = data.frame("GeneID"=phi.ID,"PhiObs"=phi.obs)

  finish.time = proc.time()
  
  Date = strsplit(as.character(Sys.time()),' ')[[1]][1]
  
  if(!phi.only) save(elr.trace,mut.trace,file=paste(out.prefix,Date,'_parm.trace',sep=''))
  save(start.time,after.burnin.time,finish.time,file=paste(out.prefix,Date,'_time.dat',sep=''))
  save(phi.dat,phi.obs.dat,file=paste(out.prefix,Date,'_phi.dat',sep=''))
  save(obs.genome,trna.dat.obs,file=paste(out.prefix,Date,'_genome.dat',sep=''))
}
