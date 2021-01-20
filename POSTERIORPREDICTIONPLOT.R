#posteriorPredictionPlot function

##### Posterior Prediction Plots #####
##### WARNING FIRST STEP CLEARS NAMESPACE #####


#### Returns a list of response times, choices, and whether the parameters violate a taboo
#### Taboos generally represent situations that could cause problems for simulation

fillDTCH = function(V1A, V1B, V1C, V2, B_ACC, B_SPEED, A, Sv)
{
  set.seed(get_seed()) # prevents it from doing the same thing every time
  predDT = predCh = rep(NA, length(V1A)) # for each actual response in the data set we will predict one
  #taboo1 = taboo2 = taboo3 = taboo4 = taboo5 = taboo6 = taboo7 = taboo8 =
  #  taboo9 = tabooYN = rep(NA, length(jab$pnum))
  #jonesAugmented$tabooYN = tabooYN
  for (i in 1:length(V1A))
    # loop through cases in jones dataset
  {
    #print(paste0("i=", i))
    if (jab$freqCondition[i] == 1 |
        jab$freqCondition[i] == 4) {
      v1 = V1A[i]
    } else if (jab$freqCondition[i] == 2 |
               jab$freqCondition[i] == 5) {
      v1 = V1B[i]
    } else if (jab$freqCondition[i] == 3 |
               jab$freqCondition[i] == 6) {
      v1 = V1C[i]
    }
    v2 = V2[i]
    if (jab$freqCondition[i] < 4)
    {
      alpha = jab$ALPHA[i]
      beta = jab$BETA[i]
    } else if (jab$freqCondition[i] > 3)
    {
      alpha = jab$BETA[i]
      beta = jab$ALPHA[i]
    }
    if (jab$speedaccuracyYN[i] == 0)
    {
      b = B_ACC[i]
    } else if (jab$speedaccuracyYN[i] == 1)
    {
      b = B_SPEED[i]
    }
    #taboo1[i] = ifelse(jonesAugmented$V1A[i] < jonesAugmented$V1B[i], 1, 0)
    #taboo2[i] = ifelse(jonesAugmented$V1A[i] < jonesAugmented$V1C[i], 1, 0)
    #taboo3[i] = ifelse(jonesAugmented$V1B[i] < jonesAugmented$V1C[i], 1, 0)
    #taboo4[i] = ifelse(jonesAugmented$V1A[i] < jonesAugmented$V2[i], 1, 0)
    #taboo5[i] = ifelse(jonesAugmented$V1B[i] < jonesAugmented$V2[i], 1, 0)
    #taboo6[i] = ifelse(jonesAugmented$V1C[i] < jonesAugmented$V2[i], 1, 0)
    #taboo7[i] = ifelse(jonesAugmented$B_ACC[i] < jonesAugmented$B_SPEED[i], 1, 0)
    #taboo8[i] = ifelse(jonesAugmented$B_SPEED[i] < jonesAugmented$A[i], 1, 0)
    #taboo9[i] = ifelse(jonesAugmented$B_ACC[i] < jonesAugmented$A[i], 1, 0)
    #tabooYN[i] = max(c(
    #  taboo1[i],
    #  taboo2[i],
    #  taboo3[i],
    #  taboo4[i],
    #  taboo5[i],
    #  taboo6[i],
    #  taboo7[i],
    #  taboo8[i],
    #  taboo9[i]
    #))
    #jonesAugmented$tabooYN[i] = tabooYN[i]
    #if (tabooYN[i] == 1)
    #{
    #  predDT[i] = NA
    #  predCh[i] = NA
    #}
    #else
    #{
    print(str(jab))
    print(paste0("v1=",v1))
    print(paste0("v2=",v2))
    print(paste0("A =",jab$A[i]))
    print(paste0("alpha=",alpha))
    print(paste0("beta=",beta))
    print(paste0("b=",b))

        pair = RBALBAinR(
        v1 = v1,
        v2 = v2,
        A = jab$A[i],
        alpha = alpha,
        beta = beta,
        b = b,
        s1 = 1,
        s2 = 1
      )
      predDT[i] = as.numeric(pair$rt)
      predCh[i] = as.numeric(pair$Ch)
    #}
  }
  return(list(
    predDT = predDT,
    predCh = predCh))#,
    #tabooYN = tabooYN
  #))
}

GETPARAMVALUE = function (pnum, paramname, burnin, lastiter)
{
  fitpnum = ifelse(pnum == 1, 1, pnum - 1)
  print(paste0("pnum=", pnum))
  print(paste0("fitpnum=", fitpnum))
  firstiter = burnin + 1
  print(paste0("burnin=", burnin))
  print(paste0("firstiter=", firstiter))
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
    #else
    #{
    #  MX = mean(theta[paramname, fitpnum, , firstiter:lastiter])
    #  return(exp(MX))
    #}
  }
}

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


posteriorPredictionPlot = function (pnum, paramname, burnin, lastiter)
{
  library(tidyverse)
  #bring in dataset of interest
  #load("/Users/OUConline/Desktop/DissertationEnds/jonesHF_v004_fit001_line524.RData")
  # the following has to be changed very carefully depending on the exact fit to be plotted
  # looking at phi_names will help
  # note that because of case 2 fitpnum does not equal pnum
  
  fitpnum2pnum = c(1, 3:17) # index of this indicates which case in target dataset (pnum)
  
  #OXI # parameter vector with c(omega,xi,iota)
  # 1,0,0 nonsequential BALBA
  # 1,0,1 GORE updating rule
  # 0,1,1 VAN-ZANDT updating rule
  
  jab = read.table (
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
  jab = jonesAugmentedBase[c(
    "pnum",
    "blocknum",
    "practiceYN",
    "speedaccuracyYN",
    "freqCondition",
    "RTsec",
    "Ch",
    "censorYN"
  )]
  jab <- filter(jab, pnum != 2)
  jab <- filter(jab, practiceYN == 0)
  summary(jab)
  attach(jab)
  table(pnum)
  
  ### set this each time based on the data file
  OXI = c(1, 0, 0) # nonsequential BALBA
  
  jab$OMEGA <- OXI[1]
  jab$XI <- OXI[2]
  jab$IOTA <- OXI[3]
  
  ### We need to add the estimated value of RHO before getting ALPHA,BETA for each trial
  for (i in 1:length(jab$pnum))
  {
    jab$RHO[i] <- GETPARAMVALUE(jab$pnum[i], "rho", burnin, lastiter)
  }
  summary(jab)
  
  jab$ALPHA[1] <- jab$OMEGA[1] + jab$XI[1]
  jab$BETA[1]  <- jab$OMEGA[1] + jab$XI[1]
  
  for (i in 2:length(jab$pnum))
  {
    if (jab$pnum[i] == jab$pnum[i - 1])
    {
      if (jab$blocknum[i] == jab$blocknum[i - 1])
      {
        if (jab$freqCondition[i] < 4)
        {
          jab$ALPHA[i] <- jab$RHO[i] * (jab$ALPHA[i -1] - jab$OMEGA[i]) + jab$OMEGA[i] + jab$IOTA[i]
          jab$BETA[i]  <- jab$RHO[i] * (jab$BETA[i- 1]  - jab$OMEGA[i]) + jab$OMEGA[i]
        }
        else
        {
          jab$ALPHA[i] <- jab$RHO[i] * (jab$ALPHA[i - 1] - jab$OMEGA[i]) + jab$OMEGA[i]
          jab$BETA[i]  <- jab$RHO[i] * (jab$BETA[i - 1]  - jab$OMEGA[i]) + jab$OMEGA[i] + jab$IOTA[i]
        }
      }
      else
      {
        jab$ALPHA[i] <- jab$OMEGA[i] +jab$XI[i]
        jab$BETA[i]  <- jab$OMEGA[i] +jab$XI[i]
      }
    }
    else
    {
      jab$ALPHA[i] <- jab$OMEGA[i] + jab$XI[i]
      jab$BETA[i]  <- jab$OMEGA[i] + jab$XI[i]
    }
  }
  firstiter = burnin + 1
  # Now we need to add the parameters as fitted or defined initially
  V1A=numeric()
  V1B=numeric()
  V1C=numeric()
  V2=numeric()
  B_ACC=numeric()
  B_SPEED=numeric()
  A=numeric()
  Sv=numeric()
  
  for (i in 1:length(jab$pnum))
  {
    V1A[i]     <- GETPARAMVALUE(jab$pnum[i], "v1_A",    burnin, lastiter)
    V1B[i]     <- GETPARAMVALUE(jab$pnum[i], "v1_B",    burnin, lastiter)
    V1C[i]     <- GETPARAMVALUE(jab$pnum[i], "v1_C",    burnin, lastiter)
    V2[i]      <- GETPARAMVALUE(jab$pnum[i], "v2",      burnin, lastiter)
    B_ACC[i]   <- GETPARAMVALUE(jab$pnum[i], "b_acc",   burnin, lastiter)
    B_SPEED[i] <- GETPARAMVALUE(jab$pnum[i], "b_speed", burnin, lastiter)
    A[i]       <- GETPARAMVALUE(jab$pnum[i], "A",       burnin, lastiter)
    Sv[i]      <- 1
  }
  
  jonesAugmented <<- jonesAugmentedBase
  # Now we can loop through and create a version of jonesAugmentedBase for each element of a sampling distribution of predicted values
  NsimsFromParms = 100
  
  
  outlist = fillDTCH(V1A, V1B, V1C, V2, B_ACC, B_SPEED, A, Sv)
  #jab$tabooYN <- outlist[[3]]
  jab$RTsim1 <- outlist[[1]]
  jab$Chsim1 <- outlist[[2]]
  ## repeat above by hand using cut and paste methods to generate a small sampling distribution of RTsim and Chsim
  
  ### posterior prediction plots ###
  ### pnum = 5 is a good case
  #pnum = 5
  
  quantileList = seq(.05, .95, .05)
  attach(jab)
  par(mfrow = c(2, 2))
  qqplot(
    RTsec[censorYN == 0 &
            speedaccuracyYN == 1 &
            Ch == 0],
    RTsim1[censorYN == 0 &
             speedaccuracyYN == 1 &
             Ch == 0],
    xlim = c(0, 3),
    ylim = c(0, 3),
    pch = 1,
    col = "darkgray",
    main = "Speed Condition Errors",
    xlab = "Observed RT",
    ylab = "Predicted RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
  qqplot(
    RTsec[censorYN == 0 &
            speedaccuracyYN == 1 &
            Ch == 1],
    RTsim1[censorYN == 0 &
             speedaccuracyYN == 1 &
             Ch == 1],
    xlim = c(0, 3),
    ylim = c(0, 3),
    pch = 1,
    col = "darkgray",
    main = "Speed Condition Correct",
    xlab = "Observed RT",
    ylab = "Predicted RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
  qqplot(
    RTsec[censorYN == 0 &
            speedaccuracyYN == 0 &
            Ch == 0],
    RTsim1[censorYN == 0 &
             speedaccuracyYN == 0 &
             Ch == 0],
    xlim = c(0, 3),
    ylim = c(0, 3),
    pch = 1,
    col = "darkgray",
    main = "Accuracy Condition Errors",
    xlab = "Observed RT",
    ylab = "Predicted RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
  qqplot(
    RTsec[censorYN == 0 &
            speedaccuracyYN == 0 &
            Ch == 1],
    RTsim1[censorYN == 0 &
             speedaccuracyYN == 0 &
             Ch == 1],
    xlim = c(0, 3),
    ylim = c(0, 3),
    pch = 1,
    col = "darkgray",
    main = "Accuracy Condition Correct",
    xlab = "Observed RT",
    ylab = "Predicted RT"
  )
  arrows(0, 0, 3, 3, len = 0, lty = 3)
}
