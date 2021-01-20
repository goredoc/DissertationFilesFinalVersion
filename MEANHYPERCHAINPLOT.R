# mean hyper chain plot function

# PHI - an array of samples from the posterior
# PHI_NAMES - a vector of parameter names
# BURNIN - number of initial iterations to discard as burnin
# NAMETAG - file name tag (usually version and RH or JHF)
# HYPERPARAMSREFERENCE - true or other reference parameters (such as the mean of the prior distribution)
# estimation_space - a vector of "Log" or "Logit", one for each parameter
# HYPER - TRUE if hyper parameters and FALSE if case parameters
# CASEID - the case to be plotted in case based (in array enumeration, not pnum)

# Still need to figure out
# how to get the right comparison line on a case plot (especially in recovery)
# file names in a case plot are not working right

meanHyperChainPlot = function(PHI,
                              PHI_NAMES,
                              BURNIN,
                              NAMETAG,
                              HYPERPARAMSREFERENCE,
                              ESTIMATION_SPACE,
                              HYPER = TRUE,
                              WITHBURNIN = TRUE,
                              CASEID = 2)
{
  if (CASEID == 1) {pnum=1}
  else (pnum=CASEID-1)
  if (HYPER==FALSE) {NAMETAG = paste0(NAMETAG,"_case_",pnum)}
  oddsMeans = (1:(length(phi_names) / 2) * 2) - 1
  if (HYPER == FALSE)
  {
    PHI = PHI[, CASEID, , ]
  }
  else
  {
    PHI = PHI[oddsMeans, , ]
  }
  
  Nparameters = dim(PHI)[1]
  Nchains     = dim(PHI)[2]
  Niters      = dim(PHI)[3]
  
  type = "hyper"
  if (HYPER == "FALSE") {
    type = "case"
  }
  
  if (WITHBURNIN == TRUE)
  {
    png(
      paste0(IDTAG, "_withburnin_", type, "_chain_means.png"),
      width = 6,
      height = 8,
      units = "in",
      res = 200
    )
  }
  else
  {
    png(
      paste0(IDTAG, "_noburnin_", type, "_chain_means.png"),
      width = 6,
      height = 8,
      units = "in",
      res = 200
    )
  }
  
  if (WITHBURNIN == TRUE)
  {
    startpoint = 1
  }
  else
  {
    startpoint = burnin + 1
  }
  
  par(mfrow = c(4, 2))
  for (i in 1:length(oddsMeans))
  {
    ymin = min(PHI[i, , startpoint:Niters], HYPERPARAMSREFERENCE[i])
    ymax = max(PHI[i, , startpoint:Niters], HYPERPARAMSREFERENCE[i])
    plot(
      PHI[i, 1, startpoint:Niters] ~ c(startpoint:Niters),
      type = "l",
      col = sample(c(
        "red", "grey", "darkgray", "lightgray"
      ), 1),
      xlab = phi_names[oddsMeans[i]],
      ylab = paste0("Estimate in ", ESTIMATION_SPACE[i], " Space"),
      ylim = c(ymin, ymax)
    )
    for (j in 2:Nchains)
    {
      points(PHI[i, j,] ~ c(1:Niters),
             type = "l",
             col = sample(c(
               "red", "grey", "lightgray", "darkgray"
             ),
             1))
    }
    arrows(
      startpoint - 1,
      HYPERPARAMSREFERENCE[i],
      Niters,
      HYPERPARAMSREFERENCE[i],
      len = 0,
      lwd = 3,
      col = "black"
    )
      arrows(
        burnin,
        mean(PHI[i, , burnin:Niters]),
        Niters,
        mean(PHI[i, , burnin:Niters]),
        len = 0,
        lwd = 3,
        col = "black",
        lty = 3
      )

  }
  dev.off()
}