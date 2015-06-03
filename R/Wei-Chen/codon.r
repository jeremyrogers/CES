# Create codon count dataframe.
# From sequences, build dataframes with codon usage for each sequence.
# Input argument = Standard FASTA file 
# Returns = codon count dataframe:
  # Standard Format 
  ##   aa    ORF  1  2  3  4  5  6
  ## 11  A YAL001C 15  9  5 23 NA NA
  ## 12  A YAL002W 18 13  3 16 NA NA
  ## 13  A YAL003W  2  6  0 24 NA NA
  ## 14  A YAL005C  1 13  0 46 NA NA
  ## 15  A YAL007C  3  5  2  6 NA NA
createCodonCount <- function(fasta.file=NULL, verbose=FALSE)
{
  if (!is.character(fasta.file)){
    stop("input 'fasta.file' must be a character string")
  } else { 
    if (verbose) 
      cat("Creating codon table from fasta file: ", fasta.file, "\n")
    
    # Read in fasta file 
    rawseqdata <- read.fasta(file = fasta.file)
  }
  
  #Calculate codon usage indices for the sequence data
  codonseq.usage <- as.data.frame(sapply(rawseqdata, seqinr::uco))
  
  ##Create subindices for codons using 'translate' function
  codon <- row.names(codonseq.usage)
  
  #codon <- as.character(unique(codonseq.usage[,1]))
  
  #aa -> amino acid
  aacodonseq.df <- data.frame(translate(sapply(as.vector(codon), s2c)), codon)
  colnames(aacodonseq.df)[1] <- "amino.acid"
  
  #Split serine into two sets and add 'Z' to set of possible factors
  levels(aacodonseq.df$amino.acid) <- c(levels(aacodonseq.df$amino.acid), "Z")
  ##Find where split serine codons are
  z <- codon %in% c("agc", "agt")
  #aacodonseq.df[z, 'amino.acid'] = 'Ser2'
  aacodonseq.df[z, "amino.acid"] = "Z"
  
  ##create list for matching codons to aa w/split ser
  codon2aa <- aacodonseq.df$amino.acid
  names(codon2aa) <- aacodonseq.df$codon
  
  ##sort by amino acid, then codon
  aacodonseq.df <- aacodonseq.df[order(aacodonseq.df$amino.acid, aacodonseq.df$codon), ]
  
  ##Create a list mapping amino acid to subindicies
  aa2subindx <- (by(aacodonseq.df$codon, aacodonseq.df$amino.acid, order))
  
  ##create list of subindexes ordered by codon
  subindx <- stack(aa2subindx[unique(aacodonseq.df$amino.acid)])
  
  ##associate with codons
  aacodonsub <- data.frame(aacodonseq.df, subindx$values)
  names(aacodonsub)[3] <- "subindx"
  
  ##Reorder based on codon so it matches original data
  aacodonsub <- aacodonsub[order(aacodonsub$codon), ]
  
  ##Create list for matching codons to subindx values
  codon2subindx <- aacodonsub$subindx
  names(codon2subindx) <- toupper(aacodonsub$codon)
  
  
  ##Transpose seqtab, add row with ORF name, and rename rows
  tseqtab <- data.frame(names(codonseq.usage), t(codonseq.usage))
  names(tseqtab)[1] <- "ORF"
  row.names(tseqtab) <- NULL
  
  ##Melt data frame and add a rows for aa and subindx
  mseqtab <- melt(tseqtab, id = "ORF")
  names(mseqtab)[c(2, 3)] <- c("codon", "count")
  
  ##map codons to aa w/split ser and codon indices
  aaval <- codon2aa[mseqtab$codon]
  subval <- codon2subindx[mseqtab$codon]
  mseqtab <- cbind(aa = aaval, subindx = subval, mseqtab)
  
  ##Remove codon column
  ##mseqtab <- mseqtab[-4];
  
  
  # preserve ordering
  levels <- unique(mseqtab$ORF)
  mseqtab$ORF <- factor(mseqtab$ORF, levels(mseqtab$ORF)[levels])
  
  
  ##cast data so that counts are in columns based on subindx
  cseqtab <- cast(mseqtab, aa + ORF ~ subindx, value = "count")
  
  ########### START SELECTING AA WITH >1 CODON #######
  ##create list of aa w >1 codon and exclude stops
  aa <- aacodonseq.df$amino.acid
  selectaa <- names(summary(aa)[(summary(aa) > 1) & ("*" != names(summary(aa)))])
  
  ##select data
  codoncount.df <- cseqtab[is.element(cseqtab$aa, selectaa), ]
  row.names(codoncount.df) <- NULL
  
  return(codoncount.df)
}



# Match gene names in codonCounts or aaData df with gene names in phi file
match_codonCount_toPhi <- function(rawcodoncounts, phi) 
{
  codoncounts <- rawcodoncounts[with(rawcodoncounts, order(ORF)), ]
  
  #Order phi df based on the 'ID' column
  phi <- phi[with(phi, order(phi[[names(phi)[1]]])), ]
  
  #Match gene names 
  m1 <- match(codoncounts$ORF, phi[, 1])
  
  #Remove NA values from codoncounts df
  filtered.codoncounts <- codoncounts[!is.na(m1), ]
  
  m2 <- match(phi[, 1], filtered.codoncounts$ORF)
  
  #Remove NA values from phi df
  filtered.phi <- phi[!is.na(m2), ]
  
  #Construct the final ouput list object 
  out.list <- list(codoncounts = filtered.codoncounts, phi = filtered.phi)
  
 return(out.list)
}



# Read in template genetic code
# sort = (logical) sort AA list by codon number
readGeneticCode <- function(file=NULL, maxAA=19, sort=FALSE, verbose=FALSE)
{
  if (is.null(file)){
    return(template.genetic.code)
  }
  
  genetic.code <- read.table(file, header=TRUE)
  
  len <- length(genetic.code)
  
  #Add column header names
#  if (len == 2)
#    names(genetic.code) <- c("aa", "codon")
#  else if (len == 3)
    names(genetic.code) <- c("aa", "codon", "value")
#  else 
#    stop("Invalid genetic code")
  
#  #replace any NAs with zeros
#  genetic.code$value[is.na(genetic.code$value)] <- 0L
  
  if (verbose)
    cat("File read and genetic code loaded ...", "\n")
  
  n.codons <- table(genetic.code$aa)
  
  #Select aa with count > 1
  n.codons <- n.codons[n.codons > 1]
  
  #Limit #of AA's analyzed.  Useful for debugging.
  if (maxAA != 19) {
    
    tmpMaxAA <- length(n.codons)
    
    if (maxAA > tmpMaxAA) {
      warning(paste("TrEff:\n  readGeneticCode(): maxAA is set to", maxAA, "which is > the expected value from the genetic code:", tmpMaxAA, ". Setting maxAA to", tmpMaxAA))
      maxAA <- tmpMaxAA
    }
    
    n.codons <- n.codons[1:maxAA]
  }
  
  ##sort AA list by number of codons for each aa
  if (sort)
    n.codons <- base::sort(n.codons, decreasing = TRUE)
  
  aa.list <- names(n.codons)
  
  genetic.code <- new("genetic.code", code=genetic.code, aa.list=aa.list, n.codons=n.codons)
  
  return( genetic.code )
}




# split codonCounts df into list df's of codon counts by aa
split_codonCounts <- function(aaList, n.codons, codoncounts)
{
  #Create list structure
  codon.counts.list <- vector("list", length(aaList))
  names(codon.counts.list) <- aaList
  
  ##extract all the codonCounts on an amino acid by amino acid basis
  for (i in seq_along(aaList)){
    codon.counts.list[[i]] <- codoncounts[codoncounts[, "aa"] == aaList[i], ]
    codon.counts.list[[i]] <- codon.counts.list[[i]][, c(2, (2 + (1:n.codons[[i]])))]
  }
  
  return(codon.counts.list)
}



# Purpose: Select a random subset of genes and return their phi and codonCountsList values
# Input: codonCountsList dataframe, phi dataframe, and scalar fraction
# Output: list of  codonCountsList (order retained from input) and phi values.
sample_codonCountsList <- function(codon.counts.list, scaled.phi, fraction=0.5, resort=TRUE, verbose=FALSE)
{
  if (verbose) 
    print("selecting fraction of codon count list...")
  
  if (fraction > 1) 
    stop(paste("ERROR in select_fraction_of_codonCountsList_and_phi()\nfraction cannot be > 1. Value given:", fraction))
  
  if (fraction < 1) {
    nOrfs <- length(codon.counts.list[[1]][, 1])
    gindx <- sample(1:nOrfs, size=round(nOrfs * fraction), replace=FALSE)
    
    ##restore original ordering if so desired
    if (resort) 
      gindx <- sort(gindx)
    
    ##extract subset of codonCountsList to be analyzed (based on gindx)
    ## on an amino by amino acid basis
    for (i in 1:length(names(codon.counts.list))) {
      codon.counts.list[[i]] <- codon.counts.list[[i]][gindx, ]
    }
    
    ##order phi appropriately
    scaled.phi <- scaled.phi[gindx, ]
  }  ## if(fraction ==1) then just return data as passed
 
 return(list(scaled.phi = scaled.phi, codon.counts.list = codon.counts.list))
}



