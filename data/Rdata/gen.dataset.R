source('../R/SE.semppr.data.load.R')
outp = outp = semppr.data.load(dna.file='fasta/observed/S.cerevisiae.S288c.fasta',elong.file='elong_rates/S.cerevisiae.2007.splitSer.tRNA',mut.file='mut_rates/flat_mutation.tsv',phi.file='expression_rates/beyer.phat.values.csv')

genome = outp[[1]]
trna.dat=outp[[2]]

phi = unlist(lapply(genome,function(list)list$phi.value))
phi.med = median(phi)
genome = lapply(genome,function(list){list$phi.value=list$phi.value/phi.med;list})

dec.order = order(phi,decreasing=TRUE)
inc.order = order(phi)

genome.decreasing.phi = list()
genome.increasing.phi = list()
for(i in 1:length(genome)){
  genome.decreasing.phi[[i]] = genome[[dec.order[i]]]
  genome.increasing.phi[[i]] = genome[[inc.order[i]]]
  
}
percentiles = seq(0,1,length.out=101)
phi.percentiles = quantile(phi,probs=percentiles)

genome.percentiles = list()
counter = 1
prcentile = 1
for(i in 1:length(genome)){
  if(genome.increasing.phi[[i]]$phi.value > phi.percentiles[prcentile]){
    counter = 1
    prcentile = prcentile+1
  }
  
  
  genome.percentiles[[paste("phi.percentile_",percentiles[prcentile],sep='')]][[counter]]=genome.increasing.phi[[i]]
  
  counter=counter+1
}

genome.percentiles.inc.length = list()
genome.percentiles.dec.length = list()
for(i in 1:length(genome.percentiles)){
  len.order = order(unlist(lapply(genome.percentiles[[i]],function(list) length(list$gene.dat$c_index))))
  len.order.dec = order(unlist(lapply(genome.percentiles[[i]],function(list) length(list$gene.dat$c_index))),decreasing=TRUE)
  for(j in 1:length(genome.percentiles[[i]])){
    genome.percentiles.inc.length[[paste("phi.percentile_",percentiles[i],sep='')]][[j]] = genome.percentiles[[paste("phi.percentile_",percentiles[i],sep='')]][[len.order[j]]]
    genome.percentiles.dec.length[[paste("phi.percentile_",percentiles[i],sep='')]][[j]] = genome.percentiles[[paste("phi.percentile_",percentiles[i],sep='')]][[len.order.dec[j]]]
    
  }
}

save(genome,file='Rdata/genome.dat')
save(genome.decreasing.phi,file='Rdata/genome.decreasing.phi.dat')
save(genome.increasing.phi,file='Rdata/genome.increasing.phi.dat')
save(genome.percentiles.dec.length,file='Rdata/genome.percentiles.dec.length.dat')
save(genome.percentiles.inc.length,file='Rdata/genome.percentiles.inc.length.dat')

glf.test.1000 = list()
counter = 1


for(i in 2:101){
  counter=0
  while(counter < 10){
    len = 0
    while(len > 1000||len < 200){
      j = ceiling(55*runif(1))
      glf.test.1000[[length(glf.test.1000)+1]] = genome.percentiles.inc.length[[i]][[j]]
      len = length(glf.test.1000[[i-1]]$gene.dat$c_index)
    }
    counter=counter+1
  }
}
glf.test.1000.dec=list()
dec.order = order(unlist(lapply(glf.test.1000,function(list)length(list$gene.dat$c_index))),decreasing=TRUE)

for(i in 1:100){
  glf.test.1000.dec[[i]] = glf.test.1000[[dec.order[i]]]
}

hist(unlist(lapply(glf.test.1000.dec,function(list) length(list$gene.dat$c_index))))

glf.test.1000 = glf.test.1000.dec

save(glf.test.1000,trna.dat,file='Rdata/glf.test.1000.dat')
