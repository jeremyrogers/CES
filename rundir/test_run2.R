source('../R/SempprExchange.main.R')
source('../R/calcScuo.R')

#Rprof()
#load('../data/Rdata/genome.dat')
load('../data/Rdata/glf.test.1000.dat')
#load('../data/Rdata/chosen.100.philength.percentiles.dat')
load('../data/Rdata/CETS.mut.trna.dat')
#sort for decreasing length
#dec.order = order(unlist(lapply(genome,function(list)length(list$gene.dat$c_index))),decreasing=TRUE)
#genome.dec.length = lapply(1:length(genome),function(i)genome[[dec.order[i]]])

glf = glf.test.1000

#Start main function
sink('output2.txt')
main.semppr(glf,trna.dat,out.prefix='1000_test/parms/',min.num.steps=500,max.num.steps=3000,simulation.type='M',proposal.type='RN',BIS=10,n.cores=22,simulate.data=TRUE,phi.only=FALSE)

sink()

#Rprof(NULL)
