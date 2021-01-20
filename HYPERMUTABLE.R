# hyperMuTable function
# makes a table of posterior
# parameter mean estimates with other descriptive statistics

hyperMuTable = function(burnin,
                        phi,
                        phi_names,
                        Niters,
                        estimation_space,
                        caseNumInArray = 0)
{
  library(kableExtra)
  oddsMeans = (1:(length(phi_names) / 2) * 2) - 1
  estimation_space = c(rep("Log", 6), "Logit")
  inverseLogit = function(x) {
    return(exp(x) / (1 + exp(x)))
  }
  # subject estimates table
  
  #pnum    = 9
  #arraynum = 8
  if (caseNumInArray == 0) {
    hyper_mu_table   = matrix(rep(NA, 7 * length(oddsMeans)), ncol = 7)
    #subject_mu_table = matrix(rep(NA, 7 * length(oddsMeans)), ncol = 7)
    colnames(hyper_mu_table) = c("Parameter",
                                 "Mean",
                                 "Median",
                                 "M-1.96SD",
                                 "M+1.96SD",
                                 "2.5ptile",
                                 "97.5ptile")
    #colnames(subject_mu_table) = c("Parameter",
    # "Mean",
    # "Median",
    # "M-1.96SD",
    # "M+1.96SD",
    # "2.5ptile",
    # "97.5ptile")
    
    roundrule = 3
    for (i in 1:length(oddsMeans))
    {
      newRow  = numeric() # for hyper parameters
      newRow[1]  = phi_names[oddsMeans[i]]
      rawValue  = mean(phi[oddsMeans[i]         , , (burnin + 1):Niters])
      newRow[2] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue),
        inverseLogit(rawValue)
      ),
      roundrule)
      rawValue   = median(phi  [oddsMeans[i]         , , (burnin + 1):Niters])
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
      
      rawValue   = as.numeric(quantile((phi[oddsMeans[i], , (burnin + 1):Niters]), c(.975)))
      newRow[7] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue),
        inverseLogit(rawValue)
      ),
      roundrule)
      hyper_mu_table[i, ] = newRow
    }
    kbl(hyper_mu_table, "latex", booktabs = T)
  }
  else{
    subject_mu_table = matrix(rep(NA, 7 * length(oddsMeans)), ncol = 7)
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
      newRow2 = numeric() # for subject parameters
      newRow2[1] = phi_names[oddsMeans[i]]
      rawValue2 = mean(theta[i , caseNumInArray, , (burnin + 1):Niters])
      
      newRow2[2] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue2),
        inverseLogit(rawValue2)
      ),
      roundrule)
      
      rawValue2  = median(theta[i, caseNumInArray, , (burnin + 1):Niters])
      
      newRow2[3] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue2),
        inverseLogit(rawValue2)
      ),
      roundrule)
      
      
      rawValue2 = mean(theta[i, caseNumInArray, , (burnin + 1):Niters]) - 1.96 * sd(theta[i, caseNumInArray, , (burnin +
                                                                                                      1):Niters])
      
      newRow2[4] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue2),
        inverseLogit(rawValue2)
      ),
      roundrule)
      
      rawValue2 = mean(theta[i, caseNumInArray, , (burnin + 1):Niters]) + 1.96 * sd(theta[i, caseNumInArray, , (burnin +
                                                                                                      1):Niters])
      
      newRow2[5] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue2),
        inverseLogit(rawValue2)
      ),
      roundrule)
      
      rawValue2 = as.numeric(quantile((theta[i, caseNumInArray, , (burnin + 1):Niters]), c(.025)))
      
      newRow2[6] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue2),
        inverseLogit(rawValue2)
      ),
      roundrule)
      
      rawValue2  = as.numeric(quantile((theta[i, caseNumInArray, , (burnin + 1):Niters]), c(.975)))
      newRow2[7] = round(ifelse(
        estimation_space[i] == "Log",
        exp(rawValue2),
        inverseLogit(rawValue2)
      ),
      roundrule)
      
      subject_mu_table[i, ] = newRow2
    }
    kbl(subject_mu_table, "latex", booktabs = T)
  }
}
