pairsPlot = function(phi,burnin, IDTAG,Nparameters,Niters,Nchains)
  
{
  library(psych)
  
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
  oddsMeans = (1:(length(phi_names) / 2) * 2) - 1
  longchainsMeans = longchains[, oddsMeans]
  colnames(longchainsMeans) = phi_names[oddsMeans]
  #pairs(longchainsMeans)
 
  png(paste0(IDTAG,"PAIRSPLOT.png"),
      width = 6,
      height=6,
      units = "in",
      res = 300
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
}