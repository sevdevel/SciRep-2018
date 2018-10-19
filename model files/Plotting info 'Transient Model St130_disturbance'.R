initialise.plotting.info <- function(plot.info,PL)
{
  #=============================================================================
  # Preparation of the way the simulation output is plotted
  #=============================================================================
  
  plot.info$spec.mfrow <- c(3,6)
  plot.info$spec.mar <- c(2,4,5,2)+0.1
  
  plot.info$reac.mfrow <- c(3,6)
  plot.info$reac.mar <- c(2,4,5,2)+0.1
  
  # conversion from mol C cm-3 of solids to percentage of solids [%]
  
  MW.C <- 12
  MW.CaCO3 <- 100.09
  C.fac <- 100/(1E06*PL$rho.sed) # from [umol cm-3 of solids] to percentage [g C / g total sed %]
  FeS.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  CaCO3.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  MnO2.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  FeOOH.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  X.Fe.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  X.NH4.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  X.Mn.fac <- 1/(PL$rho.sed) # from [umol cm-3 of solids] to [umol g-1 dry wt]
  
  plot.info$SL[["OC.f"]]$value.fac <- C.fac*MW.C 
  plot.info$SL[["OC.f"]]$value.label <- "OC fast [%]"
  
  plot.info$SL[["OC.s"]]$value.fac <- C.fac*MW.C 
  plot.info$SL[["OC.s"]]$value.label <- "OC slow [%]"
  
  plot.info$SL[["OC.r"]]$value.fac <- C.fac*MW.C 
  plot.info$SL[["OC.r"]]$value.label <- "OC refractory [%]"
  
  plot.info$SL[["MnO2"]]$value.fac <- MnO2.fac 
  plot.info$SL[["MnO2"]]$value.label <- "MnO2 [µmol g-1]"
  
  plot.info$SL[["FeOOH"]]$value.fac <- FeOOH.fac 
  plot.info$SL[["FeOOH"]]$value.label <- "FeOOH [µmol g-1]"
  
  plot.info$SL[["CaCO3"]]$value.fac <- CaCO3.fac 
  plot.info$SL[["CaCO3"]]$value.label <- "CaCO3 [µmol g-1]"
  
  plot.info$SL[["FeS"]]$value.fac <- FeS.fac 
  plot.info$SL[["FeS"]]$value.label <- "FeS [µmol g-1]"
  
  plot.info$SL[["X.Fe"]]$value.fac <- X.Fe.fac 
  plot.info$SL[["X.Fe"]]$value.label <- "X.Fe [µmol g-1]"
  
  plot.info$SL[["X.NH4"]]$value.fac <- X.NH4.fac 
  plot.info$SL[["X.NH4"]]$value.label <- "X.NH4 [µmol g-1]"
  
  plot.info$SL[["X.Mn"]]$value.fac <- X.Mn.fac 
  plot.info$SL[["X.Mn"]]$value.label <- "X.Mn [µmol g-1]"
  
  plot.info$SL[["OC"]]$value.fac <- C.fac*MW.C 
  plot.info$SL[["OC"]]$value.label <- "OC [%]"
  
  plot.info$SL[["O2"]]$depth.lim <- c(0,5)
  
  plot.info$SL[["OC.f"]]$visible <- FALSE
  plot.info$SL[["OC.s"]]$visible <- TRUE
  plot.info$SL[["OC.r"]]$visible <- FALSE
  plot.info$SL[["OC"]]$visible <- TRUE
  plot.info$SL[["MnO2"]]$visible <- TRUE
  plot.info$SL[["FeOOH"]]$visible <- TRUE
  plot.info$SL[["MnCO3"]]$visible <- TRUE
  plot.info$SL[["FeS"]]$visible <- FALSE
  plot.info$SL[["X.Fe"]]$visible <- FALSE
  plot.info$SL[["X.NH4"]]$visible <- FALSE
  plot.info$SL[["X.Mn"]]$visible <- FALSE

  plot.info$SL[["H2O"]]$visible <- FALSE
  plot.info$SL[["H"]]$visible <- FALSE
  plot.info$SL[["NH4"]]$visible <- TRUE
  plot.info$SL[["NO3"]]$visible <- FALSE
  plot.info$SL[["N2"]]$visible <- FALSE
  plot.info$SL[["SO4"]]$visible <- TRUE
  plot.info$SL[["CH4"]]$visible <- TRUE
  plot.info$SL[["HCO3"]]$visible <- TRUE
  plot.info$SL[["HS"]]$visible <- TRUE
  plot.info$SL[["Fe"]]$visible <- TRUE
  plot.info$SL[["Mn"]]$visible <- TRUE
  plot.info$SL[["Ca"]]$visible <- FALSE
  
  plot.info$RL[["Cmin.ftot"]]$visible <- FALSE
  plot.info$RL[["Cmin.stot"]]$visible <- FALSE
  plot.info$RL[["Cmin.rtot"]]$visible <- FALSE
  plot.info$RL[["ISP"]]$visible <- TRUE
  plot.info$RL[["ISD"]]$visible <- TRUE
  plot.info$RL[["CFO"]]$visible <- TRUE
  plot.info$RL[["AR.f"]]$visible <- FALSE
  plot.info$RL[["MR.f"]]$visible <- FALSE
  plot.info$RL[["SR.f"]]$visible <- FALSE
  plot.info$RL[["FR.f"]]$visible <- FALSE
  plot.info$RL[["DN.f"]]$visible <- FALSE
  plot.info$RL[["MG.f"]]$visible <- FALSE
  plot.info$RL[["AR.s"]]$visible <- FALSE
  plot.info$RL[["MR.s"]]$visible <- TRUE
  plot.info$RL[["SR.s"]]$visible <- TRUE
  plot.info$RL[["FR.s"]]$visible <- TRUE
  plot.info$RL[["DN.s"]]$visible <- FALSE
  plot.info$RL[["MG.s"]]$visible <- FALSE
  plot.info$RL[["AR.r"]]$visible <- FALSE
  plot.info$RL[["MR.r"]]$visible <- FALSE
  plot.info$RL[["SR.r"]]$visible <- FALSE
  plot.info$RL[["FR.r"]]$visible <- FALSE
  plot.info$RL[["DN.r"]]$visible <- FALSE
  plot.info$RL[["MG.r"]]$visible <- FALSE
  plot.info$RL[["FIO"]]$visible <- FALSE
  plot.info$RL[["CSO"]]$visible <- FALSE
  plot.info$RL[["MnO"]]$visible <- FALSE
  plot.info$RL[["SMn"]]$visible <- FALSE
  plot.info$RL[["CSFO"]]$visible <- FALSE
  plot.info$RL[["FMR"]]$visible <- FALSE
  plot.info$RL[["IOA"]]$visible <- FALSE
  plot.info$RL[["SIO"]]$visible <- FALSE
  plot.info$RL[["SAO"]]$visible <- FALSE
  plot.info$RL[["AMO"]]$visible <- FALSE
  plot.info$RL[["SMO"]]$visible <- FALSE
  plot.info$RL[["FMO"]]$visible <- FALSE
  plot.info$RL[["CCP"]]$visible <- TRUE
  plot.info$RL[["CCD"]]$visible <- FALSE
  
  
  
  
  return(plot.info)
}  
