rm(list = ls())#
set.seed(1)#
setwd("C:/Users/Robert Gore/Desktop/dissertation/20200930")#

load("JHF_end_v006_A.RData")#

library(Rcpp)#
library(tidyverse)#
library(truncnorm)#

sourceCpp("BALBA.cpp")#
source("JHF_end_BALBA2_v007.R")#
source("JHF_end_DEMCMC_v007.R")#
source("prepareJonesData_edition_003.R")#

delogit = function(x) {return(exp(x)/(1+exp(x)))}
arrayNumFromPNUM = function(x) 
{
  if (x==1) {return(x)}
  else {return(x-1)}
}

pnum4PPP=17
arraynum4PPP=arrayNumFromPNUM(pnum4PPP)
burnin=500
firstiter=burnin+1
lastiter=1000

rho4PPP = delogit(mean(theta[7,arraynum4PPP,,firstiter:lastiter]))


jonesData = read.csv("speedAccData4Stata.csv")#

simData     = jonesData[c("pnum",
                          "speedaccuracyYN",
                          "freqCondition",
                          "blocknum",
                          "trialNumWithinBlock")]

addAlphaBeta = function(rho=rho4PPP,dataset)
{
  attach(dataset)
  RHO = 1
  ALPHA = numeric()
  BETA = numeric()
  
  ALPHA[1]=1
  BETA[1]=1
  
  if (rho==999)
  {
    dataset$ALPHA=1
    dataset$BETA=1
    return(dataset)
    }
  
  for (i in 2:length(RTsec))
  {
    if (pnum[i] == pnum[i-1])
    {
      if (blocknum[i] == blocknum[i-1])
      {
        if (freqCondition[i] < 4)
        {
          ALPHA[i] = RHO * (ALPHA[i-1]-1) + 1 + 1  
          BETA[i]  = RHO * (BETA[i-1]-1) + 1
        }
        else
        {
          ALPHA[i] = RHO * (ALPHA[i-1]-1) + 1   
          BETA[i]  = RHO * (BETA[i-1]-1) + 1 + 1
        }
      }
      else
      {
        ALPHA[i]=BETA[i]=1
      }
    }
    else
    {
      ALPHA[i]=BETA[i]=1
    }
  }
  dataset$ALPHA=ALPHA
  dataset$BETA=BETA
  return(dataset)
}

simDataPlus = addAlphaBeta(dataset=simData)
caseSimData = simDataPlus[simDataPlus$pnum == pnum4PPP ,]#

balbaRace = function(v1, v2, A, b, alpha, beta)
  #
{
  v1i = rtruncnorm(n = 1,
                   mean = v1,
                   sd = 1,
                   a = 0)#
  v2i = rtruncnorm(n = 1,
                   mean = v2,
                   sd = 1,
                   a = 0)#
  z1i = rbeta(1, alpha, beta)#
  z2i = rbeta(1, beta, alpha)#
  t1i = (b - z1i) / v1i#
  t2i = (b - z2i) / v2i#
  rt = min(t1i, t2i)#
  corrYN = ifelse(rt == t1i, 1, 0)#
  return(list(rt = rt, corrYN = corrYN))#
}

pppQQplot = function(pnum4PPP,
                     magFactor = 10,
                     caseSimData=caseSimData,
                     jonesData=jonesData)
  {
  n <- 96 # trials per block
  num_blocks <- 20
  
  burnin = 500
  lastiter = 2000
  
  convertFromPNUM = function(x)
  {
    y = ifelse (x == 1, 1, x - 1)
    return(y)
  }
  arrayID = convertFromPNUM(pnum4PPP)#
  
  # posterior means post burnin
  pm_b_acc   = exp(mean(theta[5, arrayID, , (burnin + 1):lastiter]))#
  pm_b_speed = exp(mean(theta[6, arrayID, , (burnin + 1):lastiter]))#
  pm_v1_A    = exp(mean(theta[1, arrayID, , (burnin + 1):lastiter]))#
  pm_v1_B    = exp(mean(theta[2, arrayID, , (burnin + 1):lastiter]))#
  pm_v1_C    = exp(mean(theta[3, arrayID, , (burnin + 1):lastiter]))#
  pm_v2      = exp(mean(theta[4, arrayID, , (burnin + 1):lastiter]))#
  
  ### Not needed for real data
  N = 1#
  cond_vec <- sample(rep(1:6, 96 / 6), size = 96, replace = F)#
  
  
  #jonesData = read.csv("speedAccData4Stata.csv")#
  
  #simData     = jonesData[c("pnum",
  #                          "speedaccuracyYN",
  #                          "freqCondition",
  #                          "blocknum",
  #                          "trialNumWithinBlock")]#
  #caseSimData = simData[simData$pnum == pnum4PPP ,]#
  
  v1vector    = c(pm_v1_A, pm_v1_B, pm_v1_C, pm_v1_A, pm_v1_B, pm_v1_C)#
  
  caseSimData$v1  = v1vector[caseSimData$freqCondition]#
  caseSimData$v2  = pm_v2#
  
  A = 1#
  
  simCorrYN  = numeric()#
  simrt      = numeric()#
  simSpeedyn = numeric()#
  simfreq6   = numeric()#
  
  k = 0#
  
  for (i in 1:length(caseSimData$pnum))
  {
    for (j in 1:magFactor)
    {
      print(paste0("K=", k)) #
      k = k + 1 #
      Fb = ifelse (caseSimData$speedaccuracyYN[i] == 1, pm_b_speed, pm_b_acc) #
      Fv1 = caseSimData$v1[i] #
      Fv2 = caseSimData$v2[i] #
      FA = 1 #
      Falpha = caseSimData$ALPHA[i] #
      Fbeta = caseSimData$BETA[i] #
      Fs1 = 1 #
      Fs2 = 1 #
      pair = balbaRace(Fv1, Fv2, FA, Fb, Falpha, Fbeta) #
      simCorrYN[k] = pair$corrYN #
      simrt[k]     = pair$rt #
      simSpeedyn[k] = caseSimData$speedaccuracyYN[i] #
      simfreq6[k]  = caseSimData$freqCondition[i] #
    }
  }
  
  simDataForPPP = as.data.frame(cbind(simCorrYN, simrt, simSpeedyn, simfreq6))
  
  simRTspeedCorrect = simDataForPPP$simrt[simDataForPPP$simSpeedyn == 1 &
                                            simDataForPPP$simCorrYN == 1]
  simRTspeedWrong   = simDataForPPP$simrt[simDataForPPP$simSpeedyn == 1 &
                                            simDataForPPP$simCorrYN == 0]
  simRTaccCorrect   = simDataForPPP$simrt[simDataForPPP$simSpeedyn == 0 &
                                            simDataForPPP$simCorrYN == 1]
  simRTaccWrong     = simDataForPPP$simrt[simDataForPPP$simSpeedyn == 0 &
                                            simDataForPPP$simCorrYN == 0]
  simRTspeedCorrect = simRTspeedCorrect[simRTspeedCorrect > .180 & simRTspeedCorrect<3]
  simRTspeedWrong   = simRTspeedWrong  [simRTspeedWrong   > .180 & simRTspeedWrong  <3]
  simRTaccCorrect   = simRTaccCorrect  [simRTaccCorrect   > .180 & simRTaccCorrect  <3]
  simRTaccWrong     = simRTaccWrong    [simRTaccWrong     > .180 & simRTaccWrong    <3]
  
  jonesRTspeedCorrect = jonesData$RTsecMinusMin[jonesData$corrYN ==   1 &
                                                  jonesData$speedaccuracyYN == 1 &
                                                  jonesData$censorYN == 0 &
                                                  jonesData$pnum == pnum4PPP]
  jonesRTspeedWrong   = jonesData$RTsecMinusMin[jonesData$corrYN ==   0 &
                                                  jonesData$speedaccuracyYN == 1 &
                                                  jonesData$censorYN == 0 &
                                                  jonesData$pnum == pnum4PPP]
  jonesRTaccCorrect   = jonesData$RTsecMinusMin[jonesData$corrYN ==   1 &
                                                  jonesData$speedaccuracyYN == 0 &
                                                  jonesData$censorYN == 0 &
                                                  jonesData$pnum == pnum4PPP]
  jonesRTaccWrong     = jonesData$RTsecMinusMin[jonesData$corrYN ==   0 &
                                                  jonesData$speedaccuracyYN == 0 &
                                                  jonesData$censorYN == 0 &
                                                  jonesData$pnum == pnum4PPP]
  
  #simRTspeedCorrectQs   = quantile(simRTspeedCorrect,seq(.05,.95,.05),na.rm=TRUE)
  #jonesRTspeedCorrectQs = quantile(jonesRTspeedCorrect,seq(.05,.95,.05),na.rm=TRUE)
  #plot(simRTspeedCorrectQs~jonesRTspeedCorrectQs,type="l")
  
  png(paste0("JHF_end_v006_",pnum4PPP,".png"),
      width=6,
      height=6,
      units = "in",
      res=300
  )
  par(mfrow = c(2, 2))
  
  qqplot(simRTaccCorrect,
         jonesRTaccCorrect,
         xlim = c(0, 3),
         ylim = c(0, 3),
         main = "Accuracy Condition, Corrects",
         xlab = "Predicted RT",
         ylab = "Observed RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
  
  qqplot(simRTaccWrong,
         jonesRTaccWrong,
         xlim = c(0, 3),
         ylim = c(0, 3),
         main = "Accuracy Condition, Errors",
         xlab = "Predicted RT",
         ylab = "Observed RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
  
  qqplot(
    simRTspeedCorrect,
    jonesRTspeedCorrect,
    xlim = c(0, 3),
    ylim = c(0, 3),
    main = "Speed Condition, Corrects",
    xlab = "Predicted RT",
    ylab = "Observed RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
  
  qqplot(
    simRTspeedWrong,
    jonesRTspeedWrong,
    xlim = c(0, 3),
    ylim = c(0, 3),
    main = "Speed Condition, Errors",
    xlab = "Predicted RT",
    ylab = "Observed RT"
  )
  arrows(0, 0,3, 3, len = 0, lty = 3)
  dev.off()
  
}

pppQQplot(pnum4PPP,10,caseSimData,jonesData)


