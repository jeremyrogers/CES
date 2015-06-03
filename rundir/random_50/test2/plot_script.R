plot(log(phi.pred),log(phi.dat[1,]),xlab='log(Predicted Phi)\n',ylab='log(Real Phi)',main='Predicted Phi for 50 randomly chosen simulated genes\nPhi prior = Normal, sd = 1/10 mean',sub='Phi was predicted with other parameters constant for 50 randomly chosen genes using doubly-\n intractible MCMC methods. Evolution simulation was used to generate auxiliary variables.')
abline(lmfit <- lm(log(phi.dat[1,])~log(phi.pred)))
legend("topleft", bty="n", legend=paste("R2 =", format(summary(lmfit)$adj.r.squared, digits=4)))

for(i in 1:50){
  png(file=paste(i,'_log_phi_trace.png',sep=''))
  plot(log(phi.dat[,i]),main=paste('Phi trace for gene',i),xlab ='MCMC Step',ylab = 'Log Phi')
  abline(h=log(phi.dat[1,i]))
  dev.off()
  
}
