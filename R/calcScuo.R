calcScuo <- function(gene.list.list,trna.dat,n.cores=2){
  #Purpose: Calculate SCUO of a list of genes based on Wan et al 2006
  #Inputs: gene.list.list - list of genes--get this from semppr.data.load
  #        trna.dat - trna dataframe--get this from semppr.data.load
  #        n.cores - number of cores for mclapply
  #Outputs: SCUO - vector of SCUO values corresponding to each element of gene.list.list
  
  require(multicore)
  
  scuo = unlist(mclapply(gene.list.list,function(gene.list){
    x = matrix(nrow=19,ncol=6)
    i=1 #iterates through each codon
    j=1 #iterates through each aa
    k=1 #iterates through each syn codon for each aa
    while(i<=nrow(trna.dat)){
      
      if(trna.dat$count[i]==1 || trna.dat$aa[i]=='X'){
        i=i+1 #skip this codon
      }else{
        k=1
        n.syn.codon = trna.dat$count[i] 
        while(k<=n.syn.codon){
          codon.count = length(which(gene.list$gene.dat$c_index==trna.dat$c_index[i])) #number of codon i in codon index
          x[j,k]=codon.count
          i=i+1;k=k+1 #go onto next codon and next column in x
        }
        j=j+1 #go onto next row in x
      }
    }
    
    
    
    row.sums = rowSums(x,na.rm=TRUE)
    p = x/row.sums
    
    h = -rowSums(p*log(p),na.rm=TRUE)
    h_max = -log(1/rowSums(!is.na(x))) 
    
    o = (h_max-h)/h_max
    
    f = row.sums/sum(row.sums)
    
    sum(f*o)
    
  }))
  
  return(scuo)
  
}
