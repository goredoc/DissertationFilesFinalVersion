setwd("/Users/OUConline/Desktop/DissertationEnds")
source("MEANHYPERCHAINPLOT.R")
source("LEARNINGPLOT.R")
source("PAIRSPLOT.R")
source("HYPERMUTABLE.R")
IDTAG = "RH_end_v002_A"

load("/Users/OUConline/Desktop/DissertationEnds/RH_end_v002_A.RData")

Nparameters = dim(phi)[1]
Nchains = dim(phi)[2]
Niters  = dim(phi)[3]
burnin = 500

hyper_params_true = c(log(v1_A),
                      log(v1_B),
                      log(v1_C),
                      log(v2),
                      log(b_acc),
                      log(b_speed),
                      log(rho / (1 - rho)))
estimation_space = c(rep("Log",6),"Logit")

#### END SETUP OF VARIABLES ####

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




#### MAKE HYPERPARAMETER TRACE PLOTS ####
meanHyperChainPlot(phi,phi_names,1,IDTAG,hyper_params_true,estimation_space,TRUE,TRUE,2)
meanHyperChainPlot(phi,phi_names,1,IDTAG,hyper_params_true,estimation_space,TRUE,FALSE,2)
meanHyperChainPlot(theta,phi_names,1,IDTAG,hyper_params_true,estimation_space,FALSE,TRUE,2)
meanHyperChainPlot(theta,phi_names,1,IDTAG,hyper_params_true,estimation_space,FALSE,FALSE,2)


#### MAKE BAYESIAN LEARNING PLOTS ####
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








# for (i in 1:Nparameters)
# {
#     png(paste0(IDTAG,"_hyper_chains_",phi_names[i],".png"))
#     par(mfrow=c(1,1))
#     plot(phi[i,1,]~c(1:Niters),type="l",col=sample(c("red","grey","black","lightgray","darkgray","brown"),1),xlab=phi_names[i],ylab="Estimate in Log Space")
#     for (j in 2:Nchains)
#       {
#    	   points(phi[i,j,]~c(1:Niters),type="l",col=sample(c("red","grey","black","lightgray","darkgray","brown"),1))
#       }
#    dev.off()
# }

# this needs to be checked in each dataset until fixed in priors_init
hyper_params_true = c(log(v1_A),
                      log(v1_B),
                      log(v1_C),
                      log(v2),
                      log(b_acc),
                      log(b_speed),
                      log(rho / (1 - rho)))

estimation_space = c(rep("Log", 6), "Logit")

oddsMeans = (1:(length(phi_names) / 2) * 2) - 1
evensSDs = oddsMeans + 1
png(
  paste0(IDTAG, "_hyper_chains_means.png"),
  width = 6,
  height = 8,
  units = "in",
  res = 200
)
par(mfrow = c(4, 2))
for (i in 1:length(oddsMeans))
{
  #par(mfrow=c(1,1))
  plot(
    phi[oddsMeans[i], 1, ] ~ c(1:Niters),
    type = "l",
    col = sample(c("red", "grey", "darkgray", "lightgray"), 1),
    xlab = phi_names[oddsMeans[i]],
    ylab = paste0("Estimate in ", estimation_space[i], " Space")
  )
  for (j in 2:Nchains)
  {
    points(phi[oddsMeans[i], j, ] ~ c(1:Niters), type = "l", col = sample(c(
      "red", "grey", "lightgray", "darkgray"
    ),
    1))
  }
  arrows(
    0,
    hyper_params_true[i],
    Niters,
    hyper_params_true[i],
    len = 0,
    lwd = 3,
    col = "black"
  )
  arrows(
    (burnin + 1),
    mean(phi[oddsMeans[i], , burnin:Niters]),
    Niters,
    mean(phi[oddsMeans[i], , burnin:Niters]),
    len = 0,
    lwd = 3,
    col = "black",
    lty = 3
  )
  
}
dev.off()

longchains = matrix(rep(NA, Nparameters * Nchains * (Niters - burnin)), nrow =
                      Nchains * (Niters - burnin))
firstgood = burnin + 1
for (i in 1:Nparameters)
{
  for (j in 1:Nchains)
  {
    start = (j - 1) * (Niters - burnin) + 1
    print(paste0("start=", start))
    end = start + (Niters - burnin) - 1
    print(paste0("end=", end))
    longchains[start:end, i] = phi[i, j, firstgood:Niters]
  }
}
colnames(longchains) = phi_names

longchainsMeans = longchains[, oddsMeans]
colnames(longchainsMeans) = phi_names[oddsMeans]
pairs(longchainsMeans)

library(psych)
png(
  paste0("pairsplot_", IDTAG, ".png"),
  width = 6,
  height = 6,
  units = "in",
  res = 200
)
pairs.panels(
  longchainsMeans,
  method = "pearson",
  # correlation method
  hist.col = "#00AFBB",
  density = TRUE,
  # show density plots
  ellipses = TRUE # show correlation ellipses
)
dev.off()

library(kableExtra)
# make a table of mean estimates of each parameter across all chains
# and medians, SD, 2.5 and 97.5 percentiles
# discarding burnin

inverseLogit = function(x) {
  return(exp(x) / (1 + exp(x)))
}
hyper_mu_table = matrix(rep(NA, 7 * length(oddsMeans)), ncol = 7)
colnames(hyper_mu_table) = c("Parameter",
                             "Mean",
                             "Median",
                             "M-1.96SD",
                             "M+1.96SD",
                             "2.5ptile",
                             "97.5ptile")
roundrule = 3
for (i in 1:length(oddsMeans))
{
  newRow = numeric()
  newRow[1] = phi_names[oddsMeans[i]]
  rawValue  = mean(phi[oddsMeans[i], , (burnin + 1):Niters])
  newRow[2] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  rawValue  = median(phi[oddsMeans[i], , (burnin + 1):Niters])
  newRow[3] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  rawValue  = mean(phi[oddsMeans[i], , (burnin + 1):Niters]) - 1.96 * sd(phi[oddsMeans[i], , (burnin +
                                                                                                1):Niters])
  newRow[4] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  rawValue  = mean(phi[oddsMeans[i], , (burnin + 1):Niters]) + 1.96 * sd(phi[oddsMeans[i], , (burnin +
                                                                                                1):Niters])
  newRow[5] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  rawValue  = as.numeric(quantile((phi[oddsMeans[i], , (burnin + 1):Niters]), c(.025)))
  newRow[6] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  rawValue  = as.numeric(quantile((phi[oddsMeans[i], , (burnin + 1):Niters]), c(.975)))
  newRow[7] = round(ifelse(
    estimation_space[i] == "Log",
    exp(rawValue),
    inverseLogit(rawValue)
  ),
  roundrule)
  hyper_mu_table[i, ] = newRow
}
kbl(hyper_mu_table, "latex")


##### Posterior Prediction Plots #####
##### WARNING FIRST STEP CLEARS NAMESPACE #####
detach(package:psych)
rm(list = ls())
library(tidyverse)
#bring in dataset of interest
load("/Users/OUConline/Desktop/DissertationEnds/jonesHF_v004_fit001_line524.RData")


# the following has to be changed very carefully depending on the exact fit to be plotted
# looking at phi_names will help


# note that because of case 2 fitpnum does not equal pnum

fitpnum2pnum = c(1, 3:17) # index of this indicates which case in target dataset (pnum)

GETPARAMVALUE = function (pnum, paramname, burnin, lastiter)
{
  fitpnum = ifelse(pnum == 1, 1, pnum - 1)
  firstiter = burnin + 1
  # current coding strategy is risky
  # if the same param name is mentioned more than once in param_names
  # then this code will potentially fail silently
  # this should not happen
  estimatedParamYN = 0
  if (sum(param_names == paramname) == 1) {
    estimatedParamYN = 1
  }
  if (estimatedParamYN == 0)
  {
    if (paramname == "v1_A")    {
      return(v1_A)
    }
    else if (paramname == "v1_B")    {
      return(v1_B)
    }
    else if (paramname == "v1_C")    {
      return(v1_C)
    }
    else if (paramname == "v2")    {
      return(v2)
    }
    else if (paramname == "b_acc")   {
      return(b_acc)
    }
    else if (paramname == "b_speed") {
      return(b_speed)
    }
    else if (paramname == "rho")     {
      return(rho)
    }
    else if (paramname == "A")       {
      return(A)
    }
  }
  else
  {
    if (paramname == "rho")
    {
      MX = mean(theta["rho", fitpnum, , firstiter:lastiter])
      return(exp(MX) / (1 + exp(MX)))
    }
    else
    {
      MX = mean(theta[paramname, fitpnum, , firstiter:lastiter])
      return(exp(MX))
    }
  }
}


#OXI # parameter vector with c(omega,xi,iota)
# 1,0,0 nonsequential BALBA
# 1,0,1 GORE updating rule
# 0,1,1 VAN-ZANDT updating rule
RBALBAinR = function(v1, v2, A, alpha, beta, b, s1, s2)
{
  set.seed(get_seed())
  library(truncnorm)
  done = 0
  while (done == 0)
  {
    v1i   = rtruncnorm(
      n = 1,
      a = 0,
      b = Inf,
      mean = v1,
      sd = s1
    )
    v2i   = rtruncnorm(
      n = 1,
      a = 0,
      b = Inf,
      mean = v2,
      sd = s2
    )
    z1i   = rbeta(1, alpha, beta) * min(A, b) # this is an ad hoc rule to deal with the fact
    z2i   = rbeta(1, beta, alpha) * min(A, b) # that DEMCMC is estimating cases where b < A
    fpt1i = (b - z1i) / v1i
    fpt2i = (b - z2i) / v2i
    DTT    = min(fpt1i, fpt2i)
    Ch    = ifelse(DTT == fpt1i, 1, 0)
    if (DTT > .180 & DTT < 3) {
      done = 1
    }
  }
  return(list(rt = DTT, Ch = Ch))
}

jonesAugmentedBase = read.table (
  "SpeedAccData.txt",
  header = FALSE,
  col.names = c(
    "pnum",
    "blocknum",
    "practiceYN",
    "speedaccuracyYN",
    "stimulusID",
    "freqCondition",
    "Ch",
    "RTsec",
    "censorYN"
  )
)
jonesAugmentedBase = jonesAugmentedBase[c("pnum",
                                          "blocknum",
                                          "practiceYN",
                                          "speedaccuracyYN",
                                          "freqCondition",
                                          "RTsec",
                                          "Ch",
                                          "censorYN")]
jonesAugmentedBase = filter(jonesAugmentedBase, pnum != 2)
jonesAugmentedBase = filter(jonesAugmentedBase, practiceYN == 0)
summary(jonesAugmentedBase)
attach(jonesAugmentedBase)
table(pnum)

### set this each time based on the data file
OXI = c(1, 0, 0) # nonsequential BALBA

jonesAugmentedBase$OMEGA = OXI[1]
jonesAugmentedBase$XI = OXI[2]
jonesAugmentedBase$IOTA = OXI[3]

### We need to add the estimated value of RHO before getting ALPHA,BETA for each trial
for (i in 1:length(jonesAugmentedBase$pnum))
{
  jonesAugmentedBase$RHO[i] = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "rho", 4000, 5000)
}
summary(jonesAugmentedBase)

jonesAugmentedBase$ALPHA[1] = jonesAugmentedBase$OMEGA[1] + jonesAugmentedBase$XI[1]
jonesAugmentedBase$BETA[1] = jonesAugmentedBase$OMEGA[1] + jonesAugmentedBase$XI[1]

for (i in 2:length(jonesAugmentedBase$pnum))
{
  if (jonesAugmentedBase$pnum[i] == jonesAugmentedBase$pnum[i - 1])
  {
    if (jonesAugmentedBase$blocknum[i] == jonesAugmentedBase$blocknum[i - 1])
    {
      if (jonesAugmentedBase$freqCondition[i] < 4)
      {
        jonesAugmentedBase$ALPHA[i] = jonesAugmentedBase$RHO[i] * (jonesAugmentedBase$ALPHA[i -
                                                                                              1] - jonesAugmentedBase$OMEGA[i]) + jonesAugmentedBase$OMEGA[i] + jonesAugmentedBase$IOTA[i]
        jonesAugmentedBase$BETA[i]  = jonesAugmentedBase$RHO[i] * (jonesAugmentedBase$BETA[i -
                                                                                             1]  - jonesAugmentedBase$OMEGA[i]) + jonesAugmentedBase$OMEGA[i]
      }
      else
      {
        jonesAugmentedBase$ALPHA[i] = jonesAugmentedBase$RHO[i] * (jonesAugmentedBase$ALPHA[i -
                                                                                              1] - jonesAugmentedBase$OMEGA[i]) + jonesAugmentedBase$OMEGA[i]
        jonesAugmentedBase$BETA[i]  = jonesAugmentedBase$RHO[i] * (jonesAugmentedBase$BETA[i -
                                                                                             1]  - jonesAugmentedBase$OMEGA[i]) + jonesAugmentedBase$OMEGA[i] + jonesAugmentedBase$IOTA[i]
      }
    }
    else
    {
      jonesAugmentedBase$ALPHA[i] = jonesAugmentedBase$BETA[i] = jonesAugmentedBase$OMEGA[i] +
        jonesAugmentedBase$XI[i]
    }
  }
  else
  {
    jonesAugmentedBase$ALPHA[i] = jonesAugmentedBase$BETA[i] = jonesAugmentedBase$OMEGA[i] +
      jonesAugmentedBase$XI[i]
  }
}

# Now we need to add the parameters as fitted or defined initially
for (i in 1:length(jonesAugmentedBase$pnum))
{
  jonesAugmentedBase$V1A[i]     = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "v1_A",   4000, 5000)
  jonesAugmentedBase$V1B[i]     = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "v1_B",   4000, 5000)
  jonesAugmentedBase$V1C[i]     = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "v1_C",   4000, 5000)
  jonesAugmentedBase$V2[i]      = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "v2",     4000, 5000)
  jonesAugmentedBase$B_ACC[i]   = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "b_acc",  4000, 5000)
  jonesAugmentedBase$B_SPEED[i] = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "b_speed", 4000, 5000)
  jonesAugmentedBase$A[i]       = GETPARAMVALUE(jonesAugmentedBase$pnum[i], "A",      4000, 5000)
  jonesAugmentedBase$Sv[i]      = 1
}

detach(jonesAugmentedBase)
jonesAugmented = jonesAugmentedBase

# Now we can loop through and create a version of jonesAugmentedBase for each element of a sampling distribution of predicted values
NsimsFromParms = 100
fillDTCH = function()
{
  set.seed(get_seed())
  predDT = predCh = rep(NA, length(jonesAugmented$pnum))
  taboo1=taboo2=taboo3=taboo4=taboo5=taboo6=taboo7=taboo8=taboo9=tabooYN=rep(NA,length(jonesAugmented$pnum))
  jonesAugmented$tabooYN = tabooYN
  for (i in 1:length(jonesAugmented$pnum))
  {
    print(paste0("i=",i))
    if (jonesAugmented$freqCondition[i] == 1 |
        jonesAugmented$freqCondition[i] == 4) {
      v1 = jonesAugmented$V1A[i]
    } else if (jonesAugmented$freqCondition[i] == 2 |
               jonesAugmented$freqCondition[i] == 5) {
      v1 = jonesAugmented$V1B[i]
    } else if (jonesAugmented$freqCondition[i] == 3 |
               jonesAugmented$freqCondition[i] == 6) {
      v1 = jonesAugmented$V1C[i]
    }
    v2 = jonesAugmented$V2[i]
    if (jonesAugmented$freqCondition[i] < 4)
    {
      alpha = jonesAugmented$ALPHA[i]
      beta = jonesAugmented$BETA[i]
    } else if (jonesAugmented$freqCondition[i] > 3)
    {
      alpha = jonesAugmented$BETA[i]
      beta = jonesAugmented$ALPHA[i]
    }
    if (jonesAugmented$speedaccuracyYN[i] == 0)
    {
      b = jonesAugmented$B_ACC[i]
    } else if (jonesAugmented$speedaccuracyYN[i] == 1)
    {
      b = jonesAugmented$B_SPEED[i]
    }
    taboo1[i] = ifelse(jonesAugmented$V1A[i] < jonesAugmented$V1B[i], 1, 0)
    taboo2[i] = ifelse(jonesAugmented$V1A[i] < jonesAugmented$V1C[i], 1, 0)
    taboo3[i] = ifelse(jonesAugmented$V1B[i] < jonesAugmented$V1C[i], 1, 0)
    taboo4[i] = ifelse(jonesAugmented$V1A[i] < jonesAugmented$V2[i], 1, 0)
    taboo5[i] = ifelse(jonesAugmented$V1B[i] < jonesAugmented$V2[i], 1, 0)
    taboo6[i] = ifelse(jonesAugmented$V1C[i] < jonesAugmented$V2[i], 1, 0)
    taboo7[i] = ifelse(jonesAugmented$B_ACC[i] < jonesAugmented$B_SPEED[i], 1, 0)
    taboo8[i] = ifelse(jonesAugmented$B_SPEED[i] < jonesAugmented$A[i], 1, 0)
    taboo9[i] = ifelse(jonesAugmented$B_ACC[i] < jonesAugmented$A[i], 1, 0)
    tabooYN[i] = max(c(
      taboo1[i],
      taboo2[i],
      taboo3[i],
      taboo4[i],
      taboo5[i],
      taboo6[i],
      taboo7[i],
      taboo8[i],
      taboo9[i]
    ))
    jonesAugmented$tabooYN[i]=tabooYN[i]
    if (tabooYN[i] == 1)
    {
      predDT[i] = NA
      predCh[i] = NA
    }
    else
    {
      pair = RBALBAinR(
        v1 = v1,
        v2 = v2,
        A = jonesAugmented$A[i],
        alpha = alpha,
        beta = beta,
        b = b,
        s1 = 1,
        s2 = 1
      )
      predDT[i] = as.numeric(pair$rt)
      predCh[i] = as.numeric(pair$Ch)
    }
  }
  return(list(predDT=predDT,predCh=predCh,tabooYN=tabooYN))
}
outlist = fillDTCH()
jonesAugmented$tabooYN = outlist[[3]]
jonesAugmented$RTsim1 = outlist[[1]]
jonesAugmented$Chsim1 = outlist[[2]]
## repeat above by hand using cut and paste methods to generate a small sampling distribution of RTsim and Chsim

### posterior prediction plots ###
### pnum = 5 is a good case

pnum=5
quantileList = seq(.05,.95,.05)
attach(jonesAugmented)
par(mfrow=c(2,2))
qqplot(RTsec[censorYN==0 & speedaccuracyYN==1 & Ch==0],RTsim1[censorYN==0& speedaccuracyYN==1 & Ch==0],xlim=c(0,3),ylim=c(0,3),pch=1,col="darkgray",main="Speed Condition Errors",xlab="Observed RT", ylab="Predicted RT")
arrows(0,0,3,3,len=0,lty=3)
qqplot(RTsec[censorYN==0 & speedaccuracyYN==1 & Ch==1],RTsim1[censorYN==0& speedaccuracyYN==1 & Ch==1],xlim=c(0,3),ylim=c(0,3),pch=1,col="darkgray",main="Speed Condition Correct",xlab="Observed RT", ylab="Predicted RT")
arrows(0,0,3,3,len=0,lty=3)
qqplot(RTsec[censorYN==0 & speedaccuracyYN==0 & Ch==0],RTsim1[censorYN==0& speedaccuracyYN==0 & Ch==0],xlim=c(0,3),ylim=c(0,3),pch=1,col="darkgray",main="Accuracy Condition Errors",xlab="Observed RT", ylab="Predicted RT")
arrows(0,0,3,3,len=0,lty=3)
qqplot(RTsec[censorYN==0 & speedaccuracyYN==0 & Ch==1],RTsim1[censorYN==0& speedaccuracyYN==0 & Ch==1],xlim=c(0,3),ylim=c(0,3),pch=1,col="darkgray",main="Accuracy Condition Correct",xlab="Observed RT", ylab="Predicted RT")
arrows(0,0,3,3,len=0,lty=3)
