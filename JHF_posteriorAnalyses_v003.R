setwd("/Users/OUConline/Desktop/DissertationEnds")
#bring in dataset of interest

str(phi)
str(theta)

IDTAG="JHF_end_v003_A"

Nparameters = dim(phi)[1]
Nchains = dim(phi)[2]
Niters  = dim(phi)[3]

for (i in 1:Nparameters)
{
    png(paste0(IDTAG,"_hyper_chains",phi_names[i],".png"))
    par(mfrow=c(1,1))
    plot(phi[i,1,]~c(1:Niters),type="l",col=sample(c("red","grey","black","lightgray","darkgray","brown"),1),xlab=phi_names[i],ylab="Estimate in Log Space")
    for (j in 2:Nchains)
      {
   	   points(phi[i,j,]~c(1:Niters),type="l",col=sample(c("red","grey","black","lightgray","darkgray","brown"),1))
      }
   dev.off()
}

longchains = matrix(rep(NA,Nparameters*Nchains*(Niters-burnin)),nrow=Nchains*(Niters-burnin))
burnin = 1000
firstgood = burnin+1
for (i in 1:Nparameters)
{
	for (j in 1:Nchains)
	{
        start=(j-1)*(Niters-burnin)+1
        print(paste0("start=",start))
        end=start+(Niters-burnin)-1
        print(paste0("end=",end))
		longchains[start:end,i]=phi[i,j,firstgood:Niters]
	}
}
colnames(longchains)=phi_names
oddsMeans=(1:(length(phi_names)/2)*2)-1
longchainsMeans = longchains[,oddsMeans]
colnames(longchainsMeans)=phi_names[oddsMeans]
pairs(longchainsMeans)

library(psych)
png(paste0("pairsplot_",IDTAG,".png"))
pairs.panels(longchainsMeans, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )
dev.off()             
             