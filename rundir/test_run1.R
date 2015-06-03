source('../R/SempprExchange.main.R')
source('../R/calcScuo.R')

#Rprof()
load('../data/Rdata/genome.dat')
#load('../data/Rdata/glf.test.1000.dat')
load('../data/Rdata/CETS.mut.trna.dat')
#sort for decreasing length
dec.order = order(unlist(lapply(genome,function(list)length(list$gene.dat$c_index))),decreasing=TRUE)
genome.dec.length = lapply(1:length(genome),function(i)genome[[dec.order[i]]])

glf = genome.dec.length

#Start main function
sink('output.txt')
main.semppr(glf,trna.dat,out.prefix='genome/real/',min.num.steps=500,max.num.steps=3000,simulation.type='M',proposal.type='RN',BIS=10,n.cores=22,simulate.data=FALSE,phi.only=TRUE)

sink()

#Rprof(NULL)
