semppr.data.load <- function(dna.file,elong.file,mut.file,phi.file=NA,norm.expr=FALSE,Q=1,Ne=1,A1=4,A2=4,B=0.0025) {
  # Purpose: Easily load data in correct format for use with    
  #    calc.llik functions.                                     
  # Inputs: dna.file = name specifying dna file in fasta format. This file is required.
  #         nse.file = name specifying elongation rate file containing AA, codon nt's, 
  #           and elongation rates. NOTE: this is still elongation rate. 
  #         mut.file = name specifying mutation rate file containing AA, codon nt's, and
  #           elgonation rates. 
  #         norm.expr = flag to decide whether to scale         
  #           expression to have median 1.                        
  # Outputs: list(gene.list.full,trna.dat)                      
  # Usage: > Semppr_data_load(dna.name='/path/to/a.fasta',      
  #        + trna.dat='/path/to/b.tsv')      
  require(seqinr)
  require(stringr)
  require(plyr)
  
  #1) Read in fasta files using mikes approach.
  sc.fasta = read.fasta(file = dna.file, as.string = TRUE)
  
  num.genes = length(sc.fasta)
  
  #2) Read in parameters and other shit.  Read tRNA Files.

    ##################################### 
    # Read from elongation and mutation #
    # files                             #
    #####################################
    
    #Assume we have a header
    elong.dat = read.delim(elong.file,header=TRUE,sep = "\t",stringsAsFactors = FALSE)
    names(elong.dat) = c("AA","Codon","ElongRate")
    
    #Check if we were right--A should be first aa. if no header, one codon was chopped off
    if(length(which(elong.dat[,1]=='A')) < 4){ #Alanine has 4 codons
      elong.dat = read.delim(elong.file,header=FALSE,sep="\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","ElongRate")
    }
    
    # Check for errors
    if(length(elong.dat[,1])!=64){
      if(length(elong.dat[,1])<61){
        stop("Error: Too few codons included in mutation file.\n")
      }else if(length(elong.dat[,1])==64&&any(elong.dat[62:64,1]!=c('X','X','X'))){
        stop("Error: Stop codons must occur at the end of the elongation file.\n")
      }else if(length(which(elong.dat[,1]=='X'))==0&&length(elong.dat[,1]==61)){
        elong.dat[62,] = c('X','TAA',0)
        elong.dat[63,] = c('X','TGA',0)
        elong.dat[64,] = c('X','TAG',0)
      }else{
        stop("Unknown Error in elongation file. Elongation file must contain 64 codons with stop codons or 61 codons without stop codons. Stop codons must appear at the end of the file.")
      }
    }
    
    #Again, assume header
    mut.dat = read.delim(mut.file,header=TRUE,sep="\t",stringsAsFactors = FALSE)
    names(mut.dat) = c("AA","Codon","MutRate")
    
    #Check
    if(length(which(mut.dat[,1]=='A')) < 4){
      mut.dat = read.delim(mut.file,header=FALSE,sep="\t",stringsAsFactors=FALSE)
      names(mut.dat) = c("AA","Codon","MutRate")
    }
    
    # Check for errors
    if(length(mut.dat[,1]!=64)){
      if(length(mut.dat[,1])<61){
        stop("Error: Too few codons included in mutation file.\n")
      }else if(length(mut.dat[,1])==64&&any(mut.dat[62:64,1]!=c('X','X','X'))){
        stop("Error: Stop codons must occur at the end of the elongation file.\n")
      }else if(length(which(mut.dat[,1]=='X'))==0&&length(mut.dat[,1]==61)){
        mut.dat[62,] = c('X','TAA',0)
        mut.dat[63,] = c('X','TGA',0)
        mut.dat[64,] = c('X','TAG',0)
      }else{
        stop("Unknown Error in mutation file. Mutation file must contain 64 codons with stop codons or 61 codons without stop codons. Stop codons must appear at the end of the file.")
      }
    }
    
    #Assemble dataframe from the two dataframes
    trna.dat = elong.dat[,1:2]
    trna.dat[,3] = 0:63 #codon indices
    
    #Make sure codons are in the same order
    merged.dat = merge(elong.dat,mut.dat,by="Codon",sort=FALSE)

    trna.dat[,4] = merged.dat$ElongRate
    trna.dat[,5] = merged.dat$MutRate
    
    names(trna.dat) = c("aa","codon","c_index","elong_rate","mut_rate")
    
  
  
  # Still working with elongation rates.  trna.dat$p =trna.dat$elong/(b+trna.dat$elong)
  
  #tabulate the number of codons associated with each amino acid.
  aa.count = as.data.frame(table(trna.dat$aa))
  names(aa.count)=c("aa","count")
  trna.dat = merge(trna.dat,aa.count,by.x="aa",by.y="aa")

  #Read in the z = phi data.
  
  if(!is.na(phi.file)){
    phi.dat = read.csv(phi.file,header=TRUE,stringsAsFactors=FALSE)
    phi = phi.dat[,2]
    if(norm.expr){
      phi.dat.norm = phi.dat[,2]/median(phi.dat[,2]) # transform data so that it has median = 1
      phi=phi.dat.norm
    }
  }else{
    phi.dat=NULL #initialize
    for(i in 1:num.genes){
      name = attr(sc.fasta[[i]],'name')
      phi.dat[i] = as.numeric(str_split(name,'\t')[[1]][3])
      if(is.na(phi.dat[i])){
        stop("ERROR: Usage- must specifiy expression rate file or include expression rates in fasta according to the syntax:\n\t>GENE_ID\\tphi=\\tphi_value")
      }
      phi = phi.dat
      if(norm.expr){
        phi.dat.norm = phi.dat/median(phi.dat)
        phi=phi.dat.norm
      }
    }
  }
    
  #?  This transformation
  #should be phi.dat.norm =Ne*q*phi.  Since this is not known??? we try standardizing this way.  Seems like
  #this makes a huge difference and is a real weakness of the model.  Need to look back at old 
  #code and paper to see what this was.
  
  
  #3)Create a data frame to store codons in.  Within each gene it records each codon as a row 
  #with the coding amino acid, the codon, and the parameter values of the elongation and 
  #the translation rates.  At the optimization step will need to replace parameter values
  #with updated versions.
  
  
  #computing the size of the genome.
  lcodon = 0
  for(i in 1:num.genes){
    temp  = noquote(sc.fasta[[i]][1])
    lcodon = lcodon+nchar(temp)/3  
  }
  
  codon.paste = function(A){
    
    answer = paste(A[1],A[2],A[3],sep="")  
    
  }
  
  gene.matrix = function(sc.fasta,trna.dat,lcodon,phi){  
    #recasting as a function so that we can compile
    
    num.genes=length(sc.fasta)
    gene.list=vector(mode="list",length=num.genes)
    
    for(i in 1:num.genes){
      
      #computing number of amino acids in a codon
      temp  = s2c(toupper(noquote(sc.fasta[[i]][1])))  #make base string uppercase and and then divide into separate words.
      ltemp = length(temp)
      numb.aa  = ltemp/3
      codon.id= as.data.frame(apply(matrix(temp,numb.aa,3,byrow=TRUE),1,codon.paste))
      names(codon.id) = "codon"
      gene.dat = join(codon.id,trna.dat,by ="codon",type="left",match='all')
      
      gene.dat$c_index[numb.aa] = 99  #last row is all missing values
      name = attr(sc.fasta[[i]],'name')
      gene.list[[i]] = list(name=name,phi.value=phi[i],gene.dat=gene.dat)
      if(i%%50 ==0) print(i) 
      
    }
    answer = gene.list
  } 
  
  gene.list.full = gene.matrix(sc.fasta,trna.dat,lcodon,phi)
  
  phi.dat = unlist(lapply(gene.list.full,function(list) list$phi))
  
  
  answer = list(gene.list.full,trna.dat,phi.dat)
}