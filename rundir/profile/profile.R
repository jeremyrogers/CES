#/usr/bin/Rscript

source('../../R/SE.simEvol.R')
dyn.load('../../C/eta_simulation/SE.simulation_R-ext.so')
load('../../data/Rdata/chosen.100.philength.percentiles.dat')


Rprof(interval=0.0001)
for(i in 1:100){
out_c = simEvol(1,chosen.100.philength.percentiles[[1]]$gene.dat$c_index,trna.dat,BIS=10,GMT=0)
}
Rprof(NULL)

