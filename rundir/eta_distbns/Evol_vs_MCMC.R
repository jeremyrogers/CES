# Compare eta distributions
setwd('~/WorkingCopies/ces3/branches/exchange/rundir/')
rm(list=ls())
source('../R/SE.simEvol.R')
source('../R/SE.calcEta.R')
load('../R/data/rand.gene.list.50.dat')
library(ggplot2)

plot.list = list()
hist.dat.list = list()
counter =1

for(phi in c(0.1,0,1,10)){
  num_steps=10000
  
  
  #get eta distribution for evolution simulation
  system.time(t_c_index <- simEvol(phi,rand.gene.list.50.decreasing[[44]]$gene.dat$c_index,trna.dat,40,-10,SIMULATION.METHOD='E'))
  eta_evol = calcEta.NSE(t_c_index,trna.dat)
  for(i in 1:(num_steps-1)){
    t_c_index = simEvol(phi,t_c_index,trna.dat,1,-1,SIMULATION.METHOD='E')
    eta_evol = c(eta_evol,calcEta.NSE(t_c_index,trna.dat))
  }
  
  
  
  #get eta distribution for mcmc simulation
  t_c_index = simEvol(phi,rand.gene.list.50.decreasing[[44]]$gene.dat$c_index,trna.dat,40,0,SIMULATION.METHOD='M')
  eta_mcmc = calcEta.NSE(t_c_index,trna.dat)
  for(i in 1:(num_steps-1)){
    t_c_index = simEvol(phi,t_c_index,trna.dat,1,0,SIMULATION.METHOD='M')
    eta_mcmc = c(eta_mcmc,calcEta.NSE(t_c_index,trna.dat))
  }
  
  hist.dat = rbind(data.frame("eta"=eta_mcmc,"category"=rep('MCMC',length.out=num_steps)),data.frame("eta"=eta_evol,"category"=rep('Evol',length.out = num_steps)))
  #hist.dat2 = data.frame("mcmc"=eta_mcmc,'evol'=eta_evol)
  hist.dat.list[[counter]]=hist.dat
  
  
  plot.list[[counter]] = ggplot(data=hist.dat,aes(x=eta,fill=category))+geom_histogram(alpha=0.6,position='identity')+ggtitle(paste('Simulated Eta for Random Gene\naa_count = ',length(rand.gene.list.50.decreasing[[44]]$gene.dat$c_index),' phi = ',phi,' num_steps = ',num_steps,'\nMCMC mean = ',round(mean(eta_mcmc),2),' Evol mean = ',round(mean(subset(hist.dat,category=='Evol')[,1]),2),sep=''))
  counter = counter+1
  
}