gelman.conv <- function(parm.trace,LOG=TRUE){
  #Purpose: This function assesses convergence of a single chain of MCMC posterior draws based on 
  #         a modified version of criteria suggested by Gelman et al. 2004 (pg. 294). The suggested criteria
  #         starts several indpendent chains at random points and compares variance between and within
  #         chains. This method differs by dividing one chain into 10 chains that are treated as 
  #         independent. This method reduces runtime because starting multiple mcmc chains can be conputationally
  #         intensive. However, since these chains are not independent, this method can underestimate between-chain
  #         variance and determine convergence prematurely. This problem has not manifested in my experience and 
  #         can be remedied by altering the R.hat necessary to determine convergence, so this method remains a 
  #         useful alternative to starting multiple chains.
  #Inputs: parm.trace - vector containing raw phi values for one gene from MCMC
  #Outputs: List containing:
  #          $R.hat - factor used to assess convergence. Gelman suggests R.hat threshold of 1.1, meaning
  #                 you can stop simulation when R.hat goes below 1.1. Gelman also notes that R.hat 
  #                 threshold can be lowered for more stringent convergence criteria.
  #          $mean.all - mean across all traces
  #          $var.all - variance across all traces
  
  len = length(parm.trace)
  discard = len/2 + len%%20/2 #discard first half plus additional to allow the trace
                              #to divide evenly into 10 traces
  if(LOG){
    parm.mat = log(matrix(parm.trace[-(1:discard)],len/20))
  }else{
    parm.mat = matrix(parm.trace[-(1:discard)],len/20)
  }
    
  m = ncol(parm.mat) # number of traces
  n = nrow(parm.mat) # number of samples in each trace
  
  within.mean <- apply(parm.mat,2,mean)
  within.sq.mean <- apply(parm.mat^2,2,mean)
  
  mean.all <- mean(within.mean)
  mean.sq.all <- mean(within.sq.mean)
  
  #Within Variance
  W = sum(within.sq.mean - within.mean^2)/(m)
  
  #Between variance  
  B = n/(m-1)*sum((within.mean-mean.all)^2)
  
  #var^+
  var.plus = W + 1/n*B
  
  R.hat = sqrt(var.plus/W)
  
  list(R.hat = R.hat,mean.all=mean.all,var.all=mean.sq.all - mean.all^2)
}

