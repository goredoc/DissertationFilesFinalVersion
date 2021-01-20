setwd("/Users/OUConline/Desktop/DissertationEnds")
source("MEANHYPERCHAINPLOT.R")
source("LEARNINGPLOT.R")
source("PAIRSPLOT.R")
source("HYPERMUTABLE.R")
source("POSTERIORPREDICTIONPLOT.R")

IDTAG="JHF_end_v006_A"
load( "JHF_end_v006_A.RData")

Nparameters = dim(phi)[1]
Nchains = dim(phi)[2]
Niters  = dim(phi)[3]

burnin=500

hyper_params_starts=c(log(v1_A),
                      log(v1_B),
                      log(v1_C),
                      log(v2),
                      log(b_acc),
                      log(b_speed),
                      log(rho / (1 - rho)))
estimation_space=c(rep("Log",6),"Logit")

library(asbio)
hyper.v1_A.mu    = phi[ 1,,]
hyper.v1_A.sd    = phi[ 2,,]
hyper.v1_B.mu    = phi[ 3,,]
hyper.v1_B.sd    = phi[ 4,,]
hyper.v1_C.mu    = phi[ 5,,]
hyper.v1_C.sd    = phi[ 6,,]
hyper.v2.mu      = phi[ 7,,]
hyper.v2.sd      = phi[ 8,,]
hyper.b_acc.mu   = phi[ 9,,]
hyper.b_acc.sd   = phi[10,,]
hyper.b_speed.mu = phi[11,,]
hyper.b_speed.sd = phi[12,,]
hyper.rho.mu     = phi[13,,]
hyper.rho.sd     = phi[14,,]

R.hat (hyper.v1_A.mu    , burn.in=burnin/Niters)
R.hat (hyper.v1_B.mu    , burn.in=burnin/Niters)
R.hat (hyper.v1_C.mu    , burn.in=burnin/Niters)
R.hat (hyper.v2.mu      , burn.in=burnin/Niters)
R.hat (hyper.b_acc.mu   , burn.in=burnin/Niters)
R.hat (hyper.b_speed.mu , burn.in=burnin/Niters)
R.hat (hyper.rho.mu     , burn.in=burnin/Niters)

R.hat (hyper.v1_A.sd    , burn.in=burnin/Niters)
R.hat (hyper.v1_B.sd    , burn.in=burnin/Niters)
R.hat (hyper.v1_C.sd    , burn.in=burnin/Niters)
R.hat (hyper.v2.sd      , burn.in=burnin/Niters)
R.hat (hyper.b_acc.sd   , burn.in=burnin/Niters)
R.hat (hyper.b_speed.sd , burn.in=burnin/Niters)
R.hat (hyper.rho.sd     , burn.in=burnin/Niters)


library(loo)
postBurninLL = hyper_weight[,(burnin+1):Niters]
loo::waic(postBurninLL)
loo::loo(postBurninLL)

meanHyperChainPlot(phi,phi_names,burnin,IDTAG,hyper_params_starts,estimation_space,TRUE,TRUE,10)
meanHyperChainPlot(phi,phi_names,burnin,IDTAG,hyper_params_starts,estimation_space,TRUE,FALSE,10)

### making a LearningPlot
### samples are a histogram
### prior is a normal curve overlay in the logs
png(
  paste0(IDTAG, "_learningPlots_means.png"),
  width = 6,
  height = 8,
  units = "in",
  res = 300
)
par(mfrow = c(4, 2))
for (j in 1:length(oddsMeans))
{
  learningPlot(burnin,j) 
}  
dev.off()

pairsPlot(phi,burnin,IDTAG,Nparameters,Niters,Nchains)




#### MAKE HYPERPARAMETER ESTIMATES TABLE ####
source("HYPERMUTABLE.R")
hyperMuTable(burnin,phi,phi_names,Niters,estimation_space)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,1)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,2)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,3)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,4)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,5)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,6)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,7)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,8)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,9)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,10)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,11)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,12)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,13)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,14)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,15)
hyperMuTable(burnin,theta,phi_names,Niters,estimation_space,16)










#### MAKE POSTERIOR PREDICTION PLOTS ####
posteriorPredictionPlot(1,"v1_A",500,2000)











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
burnin = 500
longchains = matrix(rep(NA,Nparameters*Nchains*(Niters-burnin)),nrow=Nchains*(Niters-burnin))
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
      


library(kableExtra)
# make a table of mean estimates of each parameter across all chains
# and medians, SD, 2.5 and 97.5 percentiles
# discarding burnin
estimation_space = c(rep("Log",6),"Logit")
inverseLogit = function(x) {
  return(exp(x) / (1 + exp(x)))
}
# subject estimates table

pnum    =6
arraynum=5

hyper_mu_table   = matrix(rep(NA, 7 * length(oddsMeans)), ncol = 7)
subject_mu_table = matrix(rep(NA, 7 * length(oddsMeans)), ncol = 7)
colnames(hyper_mu_table) = c("Parameter",
                             "Mean",
                             "Median",
                             "M-1.96SD",
                             "M+1.96SD",
                             "2.5ptile",
                             "97.5ptile")
colnames(subject_mu_table) = c("Parameter",
                               "Mean",
                               "Median",
                               "M-1.96SD",
                               "M+1.96SD",
                               "2.5ptile",
                               "97.5ptile")

roundrule = 3
for (i in 1:length(oddsMeans))
{
  newRow  = numeric() # for hyper parameters
  newRow2 = numeric() # for subject parameters
  newRow[1]  = phi_names[oddsMeans[i]]
  newRow2[1] = phi_names[oddsMeans[i]]
  rawValue  = mean(  phi[oddsMeans[i]         ,,(burnin + 1):Niters])
  rawValue2 = mean(theta[          i ,arraynum,,(burnin + 1):Niters])
  newRow[2] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  
  newRow2[2] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue2),
    inverseLogit(rawValue2)
  ),
  roundrule)
  
  rawValue   = median(phi  [oddsMeans[i]         , , (burnin + 1):Niters])
  rawValue2  = median(theta[           i,arraynum, , (burnin + 1):Niters])  
  newRow[3] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  
  newRow2[3] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue2),
    inverseLogit(rawValue2)
  ),
  roundrule)
  
  
  rawValue  = mean(phi[oddsMeans[i], , (burnin + 1):Niters]) - 1.96 * sd(phi[oddsMeans[i], , (burnin +
                                                                                                1):Niters])
  rawValue2 = mean(theta[i,arraynum, , (burnin + 1):Niters]) - 1.96 * sd(theta[i,arraynum, , (burnin +
                                                                                                                      1):Niters])
  newRow[4] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  newRow2[4] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue2),
    inverseLogit(rawValue2)
  ),
  roundrule)
  
  rawValue  = mean(phi[oddsMeans[i], , (burnin + 1):Niters]) + 1.96 * sd(phi[oddsMeans[i], , (burnin +
                                                                                                1):Niters])
  rawValue2 = mean(theta[i,arraynum, , (burnin + 1):Niters]) + 1.96 * sd(theta[i,arraynum, , (burnin +
                                                                                                                      1):Niters])
  
  newRow[5] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  
  newRow2[5] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue2),
    inverseLogit(rawValue2)
  ),
  roundrule)
  
  
  rawValue  = as.numeric(quantile((phi[oddsMeans[i], , (burnin + 1):Niters]), c(.025)))
  rawValue2 = as.numeric(quantile((theta[i,arraynum, , (burnin + 1):Niters]), c(.025)))
  newRow[6] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  
  newRow2[6] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue2),
    inverseLogit(rawValue2)
  ),
  roundrule)
  
  rawValue   = as.numeric(quantile((phi[oddsMeans[i], , (burnin + 1):Niters]), c(.975)))
  rawValue2  = as.numeric(quantile((theta[i,arraynum, , (burnin + 1):Niters]), c(.975)))  
  newRow[7] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  newRow2[7] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue2),
    inverseLogit(rawValue2)
  ),
  roundrule)
  
  hyper_mu_table[i, ] = newRow
  subject_mu_table[i,] = newRow2
}
kbl(hyper_mu_table, "latex")
#kbl(subject_mu_table, "latex")

