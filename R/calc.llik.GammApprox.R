path.string <- strsplit(getwd(),'/')[[1]]

if(path.string[length(path.string)-1]!='exchange'){
  warning('Warning: you must load simulation_R-ext.so manually unless you are in a directory one below exchange.\nExample: exchange/R, exchange/preston, exchange/rundir\n')
}else{ dyn.load('../C/eta_simulation/SE.simulation_R-ext.so')}

rm(path.string)

calc.llik.GammApprox <- function(phi,gene.list,trna.dat,A_1 = 4,A_2 = 4,B = 0.0025,Ne = 1,Q = 1,MES=0){
  #A_1 = 4;A_2 = 4;B = 0.0025;Ne = 1;Q = 1;MES=0
  source('../R/SE.calc.eta.moments.R')
  require(gsl)
  
  outp<-calc.moments.eta(gene.list,trna.dat,A_1,A_2,B,Ne,Q,MES)
  
  eta_mean = outp[1]
  eta_var = outp[2]
  eta_min = outp[3]
  eta_max = outp[4]
  eta_obs = outp[5]
  
  ########################
  # Gamma Distbn Moments #
  ########################
  
  alpha = as.numeric((eta_mean-eta_min)^2/eta_var)
  bta = as.numeric((eta_mean-eta_min)/eta_var)
  
  ##########################
  # Likelihood Calculation #
  ##########################
  
  #Taken from distn.of.eta.values.pdf
  # Here, z is equal to y, which is q*Ne*phi. 
  # f(z|\eta,\vec{\lambda})=\frac{e^{(-z-\beta)(\eta-\eta_{min})}(z+\beta)^\alpha(\eta-\eta_{min})^{\alpha+1}}{\Gamma{1+\alpha,\beta(\eta-\eta_{min})}}
  # Calculated values verified with msaum's semppr code--they are consistent with his calculations
  
  z = Q*Ne*phi # We also call this y
  #incomp.gamma = pgamma(1+alpha,bta*(eta_obs-eta_min),lower.tail=FALSE)*gamma(1+alpha) #See ?pgamma --> Note
  #the above line used to work but now gamma(1+alpha) returns Inf sometimes. Replaced this line with gamma_inc 
  #function from gsl package
  
  llik = (-z-bta)*(eta_obs-eta_min) + alpha*log(z+bta) + (alpha+1)*log(eta_obs - eta_min) - log(gamma_inc(bta*(eta_obs-eta_min),1+alpha))
  
  ###########################
  # Conditional Dist of Eta #
  ###########################
  #Evaluate f(eta|lambda) at midpoint of each bin
  #num.bins = 20
  #bin.lims = seq(eta_min,eta_max,length.out=num.bins+1)
  #intv = (eta_max-eta_min)/20
  #eta.midpts = seq(eta_min+1/2*intv,eta_max-1/2*intv,length.out=num.bins)
  #dim(eta.midpts) <- c(20,1)
  
  #f.eta.g.lambda <- function(eta,eta.min,alpha,bta,z) {
  #  num = exp((-z-bta)*(eta-eta.min))*(z+bta)^alpha*(eta-eta.min)^(alpha+1)
  #  denom = pgamma(1+alpha,bta*(eta-eta.min),lower.tail=FALSE)*gamma(1+alpha)
  #  
  #  num/denom
  #}
  
  #bins = apply(X=eta.midpts,MARGIN=1,FUN=f.eta.g.lambda,eta.min=eta_min,alpha=alpha,bta=bta,z=z)
  #bins = bins/sum(bins)

  
  #answer = list(lik=c(llik = llik,eta.obs = eta_obs,eta.mean = eta_mean,eta.var = eta_var,eta.min=eta_min,eta.max=eta_max),hist=list(bins=bins,lims=bin.lims))
  c(llik=llik)
}
