setwd('~/WorkingCopies/ces3/branches/exchange/rundir/')
library(ggplot2)

load('chosen_100/3000_steps/BIS10/reflnorm_prop/norand/2014-04-08_genome.dat')
load('chosen_100/3000_steps/BIS10/reflnorm_prop/norand/2014-04-08_time.dat')
load('chosen_100/3000_steps/BIS10/reflnorm_prop/norand/2014-04-08_phi.dat')
source('../R/SE.calc.eta.moments.R')

real.phi = phi.dat[1,]
#take off first half of the samples and thin by every fifth sample
phi.dat.thin = phi.dat[ceiling(length(phi.dat[,1])/2):length(phi.dat[,1]),]
phi.dat.thin = phi.dat.thin[seq(1,length(phi.dat.thin[,1]),by=5),]

pred.phi = apply(phi.dat.thin,2,mean)

lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

#counter = 0
#ind = 0
#for(i in 1:length(genome)){
#  output = calc.moments.eta(genome[[i]],trna.dat)
#  
#  if((output[1]-output[5])/output[2] > 0.2 && length(genome[[i]]$gene.dat$c_index) > 300){
#    counter = counter+1
##    ind[counter] = i
#  }
#  
#}

#real.phi = phi.dat[1,ind]
#pred.phi = apply(phi.dat,2,mean)[ind]

plot.data = data.frame("y"=(pred.phi),"x"=(real.phi))
eqn.txt = lm_eqn(plot.data)
plot.data = data.frame("Predicted_Phi"=(pred.phi),"Phi_obs"=(real.phi),"length"=unlist(lapply(genome,function(list)length(list$gene.dat$c_index))))
ggplot(plot.data,aes(x=Phi_obs,y=Predicted_Phi))+geom_point(aes(color=length))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_text(aes(x=20,y=60),label=eqn.txt,parse=TRUE)+
  ggtitle(label=paste('Predicted vs Observed phi for 100 simulated genes. \nRuntime:',(finish.time-start.time)[3],'sec'))+
  geom_abline(intercept=0,slope=1,color='green')

plot.data = data.frame("y"=log(pred.phi),"x"=log(real.phi))
eqn.txt = lm_eqn(plot.data)
plot.data = data.frame("Log_Predicted_Phi"=log(pred.phi),"Log_Phi_obs"=log(real.phi),"length"=unlist(lapply(genome,function(list)length(list$gene.dat$c_index))))
ggplot(plot.data,aes(x=Log_Phi_obs,y=Log_Predicted_Phi))+geom_point(aes(color=length))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_text(aes(x=0,y=4),label=eqn.txt,parse=TRUE)+
  ggtitle(label=paste('Predicted vs Observed Log(phi) for 100 simulated genes. \nRuntime:',(finish.time-start.time)[3],'sec'))+
  geom_abline(intercept=0,slope=1,color='green')

#Phi traces and hist
for(i in 1:100){
  output = calc.moments.eta(genome[[i]],trna.dat)
  
  png(file=paste('hist/',i,'_phi_hist.png',sep=''),width=530,height=530)
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
  hist((phi.dat[,i]),main=paste('Phi Hist'),xlab = 'Phi')
  abline(v=(mean(phi.dat[1,i])),col='red')
  hist(log(phi.dat[,i]),main=paste('Log Phi Hist'))
  plot(log(phi.dat[,i]),main=paste('Log Phi Trace'),xlab='',sub=paste('Length:',length(genome[[i]]$gene.dat$c_index),'\tPhi_obs:',round(real.phi[i],2),'\tPhi_pred',round(pred.phi[i],2),'\t(eta.mean-eta.obs)/eta.var:',round((output[1]-output[5])/output[2],3),'\neta.mean\teta.var\teta.min\teta.max\teta.obs\n',round(output[1],2),'\t',round(output[2],2),'\t',round(output[3],2),'\t',round(output[4],2),'\t',round(output[5],2)))
  dev.off()
  
}
