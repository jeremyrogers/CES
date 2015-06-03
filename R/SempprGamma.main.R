
main.semppr <- function(obs.genome,trna.dat,out.prefix,proposal.type = 'RN',n.cores = 2){
  require(multicore)
  
  source('../R/calc.llik.GammApprox.R')
  source('../R/calcScuo.R')
  source('../R/SE.simEvol.R')
  source('../R/post.mcmc.R')
  
  dyn.load('../C/eta_simulation/SE.simulation_R-ext.so')
  
  #create output directory if it doesn't exist
  dir.create(out.prefix)
  
  start.time = proc.time() 
  
  pop.parms = list(Q = 1, Ne = 1, A1 = 4, A2 =4, B = 0.0025) #population parameters
  genome.end=list()
  for(i in 1:length(obs.genome)){
    obs.genome[[i]]$phi.obs = obs.genome[[i]]$phi.value  
    
  }
  
  genome = mclapply(obs.genome, function(gene.list,trna.dat,simulation.type,BIS,GMT,MES){gene.list$aux_c_index = simEvol(phi=gene.list$phi.obs,c_index=gene.list$gene.dat$c_index,trna.dat=trna.dat,BIS=100,GMT=0,MES=0,SIMULATION.METHOD='M');
                                                                                         gene.list$gene.dat$codon=NULL        # These next few lines reduce
                                                                                         gene.list$gene.dat$aa=NULL           # the size of the genome object
                                                                                         gene.list$gene.dat$elong_rate=NULL;  
                                                                                         gene.list$gene.dat$mut_rate=NULL
                                                                                         gene.list$gene.dat$count=NULL
                                                                                         gene.list$gene.dat$c_index=as.integer(gene.list$gene.dat$c_index)
                                                                                         gene.list$aux_c_index = as.integer(gene.list$aux_c_index)
                                                                                         gene.list},trna.dat,simulation.type,BIS,GMT,MES,mc.cores=n.cores,mc.preschedule=TRUE)
  num.genes = length(genome)
  
  ########################################################################################################################
  #1) USE SIMULATED SEQUENCE--SET THE OBSERVED SEQUENCE = INITIAL BURNED IN SEQUENCE                                   #
  #                                                                                                                      #
  genome = lapply(genome, function(gene){gene$gene.dat$c_index[-length(gene$gene.dat$c_index)] =gene$aux_c_index;gene})  #
  #                                                                                                                      #
  ########################################################################################################################
  
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
  
  for(i in 1:length(genome)){
    genome[[i]]$phi.value=phi.dat[i]
  }
  
  mcmc.outp=lapply(7,function(i) post.mcmc(y.start=genome[[i]]$phi.obs,prop.var=1,LIK.FUN=calc.llik.GammApprox,gene.list=genome[[i]],trna.dat=trna.dat,tot.steps=1000))
  
  #optim.list=lapply(genome,function(gene.list) optimize(calc.llik.GammApprox,interval=c(0,500),gene.list=gene.list,trna.dat=trna.dat,maximum=TRUE))
  outp=list()
  

  
  phi.est.mle=unlist(lapply(optim.list,function(list)list$maximum))
  phi.real=unlist(lapply(genome,function(list)list$phi.value))
  
  Date = strsplit(as.character(Sys.time()),' ')[[1]][1]
  save(phi.est.mle,phi.real,file=paste(out.prefix,Date,'_phi.dat',sep=''))
}