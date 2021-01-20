# LEARNINGPLOT

learningPlot = function(burnin,i)
{
  burnin = 500
  oddsMeans = (1:(length(phi_names) / 2) * 2) - 1
  samples    = phi[(i*2)-1,,(burnin+1):length(phi[1,1,])]
  priorMu    = mu_mean_vec[i]
  priorSigma = mu_sd_vec[i]
  loci  = priorMu - 2*priorSigma
  hici  = priorMu + 2*priorSigma
  losamp= min(samples)
  hisamp= max(samples)
  lowx = min(loci,losamp)
  highx= max(hici,hisamp)
  x=seq(lowx,highx,.001)
  y=dnorm(x,priorMu,priorSigma)
  hist(samples,xlim=c(lowx,highx),freq=FALSE,col="red",main=(paste0("Histogram with Prior Curve, ", phi_names[(i*2)-1])),xlab="Posterior Samples (Post-Burn-In)")
  arrows(mean(samples),0,mean(samples),10,len=0,lwd=1,col="white")
  points(y~x,xlim=c(lowx,highx),lty=3,col="gray",type="l",lwd=1)
  arrows(priorMu,0,priorMu,10,len=0,col="gray",lwd=1,lty=3)
}
