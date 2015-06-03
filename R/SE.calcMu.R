calcMu <- function(c_index,trna.dat){
  n.codons=61
  #1) Calculate codon counts 
  codon.counts = unlist(lapply(1:n.codons,FUN=function(i) length(which(c_index==i))))
  #2) calc mu
  mut.rates = as.numeric(trna.dat$mut_rate[sort.list(trna.dat$c_index)])
  prod(mut.rates[1:n.codons]^codon.counts)
  
}