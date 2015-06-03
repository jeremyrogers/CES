#### runtime vs. length for evolution and MCMC gene simulation ####
rm(list=ls())
load('../R/data/rand.gene.list.50.dat')
source("../R/SE.simEvol.R")
library(ggplot2)

gene.lengths = gene.lengths = unlist(lapply(rand.gene.list.50.decreasing,function(list) length(list$gene.dat$c_index)))

runtime.evol = rep(0,50)

for(i in 1:50){
  
  time = system.time(simEvol(phi.dat.decreasing[i],rand.gene.list.50.decreasing[[i]]$gene.dat$c_index,trna.dat,BIS=0,GMT=0,MES=10000,SIMULATION.METHOD='E'))
  runtime.evol[i] = time[3]
  
}


runtime.mcmc = rep(0,50)

for(i in 1:50){
  time = system.time(simEvol(phi.dat.decreasing[i],rand.gene.list.50.decreasing[[i]]$gene.dat$c_index,trna.dat,BIS=0,GMT=0,MES=10000,SIMULATION.METHOD='M'))
  runtime.mcmc[i] = time[3]  
}


plot.data = data.frame("GeneLength"=c(gene.lengths,gene.lengths),"Runtime_per_1000_steps"=c(runtime.evol,runtime.mcmc),"SimulationType"=c(rep("Evol",50),rep("MCMC",50)))
plot.obj = ggplot(plot.data,aes(x=GeneLength,y=Runtime_per_1000_steps)) + geom_point(aes(color=SimulationType))
lmfit.evol = lm(runtime.evol~gene.lengths)
coef.evol = coef(lmfit.evol)[1]
coef.evol[2] = coef(lmfit.evol)[2]
lmfit.mcmc = lm(runtime.mcmc~gene.lengths)
coef.mcmc = coef(lmfit.mcmc)[1]
coef.mcmc[2] = coef(lmfit.mcmc)[2]

lm_eqn = function(df){
  m = lm(RT ~ GL, df);
  eq <- substitute(italic(RT) == a + b %.% italic(GL)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

evol.eqn = lm_eqn(data.frame('GL'=gene.lengths,'RT'=runtime.evol))
mcmc.eqn = lm_eqn(data.frame('GL'=gene.lengths,'RT'=runtime.mcmc))

plot.obj=plot.obj + geom_abline(aes(slope=coef.evol[2],intercept=coef.evol[1])) + geom_abline(aes(slope=coef.mcmc[2],intercept=coef.mcmc[1]))

jpeg(file='runtime/Runtime_plot_random_50.jpg',width=725,height=650)
plot.obj + geom_text(size=4,x=350,y=1.5,label=paste("Evol:",evol.eqn),parse=TRUE) + geom_text(size=4,x=350,y=1.4,label=paste("MCMC:",mcmc.eqn),parse=TRUE)+ggtitle(paste('Runtime vs. Gene Length for 50 randomly chosen genes.\nRuntime of Evolution simulation scales ',round(coef.evol[2]/coef.mcmc[2],digits=2),' faster than MCMC simulation with respect to gene length.',sep=''))
dev.off()
