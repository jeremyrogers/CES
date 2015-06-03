opt.seq <- function(gene.list,trna.dat){
  aa.seq <- gene.list$gene.dat$aa
  
  if(is.na(aa.seq[length(aa.seq)])){
    aa.seq <- aa.seq[-length(aa.seq)]
  }
  
  #### This basically looks for the codon index with the best elongation prs 
  c.seq.opt <- unlist(lapply(aa.seq,FUN = function(aa) (trna.dat$c_index[which(trna.dat$aa==aa)])[which(trna.dat$elong_pr[which(trna.dat$aa==aa)] == max(trna.dat$elong_pr[which(trna.dat$aa==aa)]))[1]]))
  
  #return
  c.seq.opt
}

pess.seq <- function(gene.list,trna.dat){
  aa.seq <- gene.list$gene.dat$aa
  
  if(is.na(aa.seq[length(aa.seq)])){
    aa.seq <- aa.seq[-length(aa.seq)]
  }
  
  #### This basically looks for the codon index with the worst elongation prs 
  c.seq.pess <- unlist(lapply(aa.seq,FUN = function(aa) (trna.dat$c_index[which(trna.dat$aa==aa)])[which(trna.dat$elong_pr[which(trna.dat$aa==aa)] == min(trna.dat$elong_pr[which(trna.dat$aa==aa)]))[1]]))
  
  #return
  c.seq.pess
}

rand.seq <- function(gene.list,trna.dat){
  aa.seq <- gene.list$gene.dat$aa
  
  if(is.na(aa.seq[length(aa.seq)])){
    aa.seq <- aa.seq[-length(aa.seq)]
  }
  
  rand.seq <- unlist(lapply(aa.seq,FUN = function(aa){ candidates <- which(trna.dat$aa==aa)
                                                       (trna.dat$c_index[candidates])[ceiling(runif(n=1)*length(candidates))]
                                                       }))
  
  #return
  rand.seq
}