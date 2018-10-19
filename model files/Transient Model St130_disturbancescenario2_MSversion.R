###############################################################################
# Model of early diagenesis
# Type: Carbon,sulfur and iron cycling  
# Author: Filip Meysman (filip.meysman@nioz.nl)
#         Sebastiaan van de Velde (sevdevel@vub.ac.be) 
# Affiliation: NIOZ, Korringaweg 7, 4401 NT, Yerseke 
#              VUB, Pleinlaan 2, 1050 BE, Brussel
###############################################################################

# public packages
require(ReacTran)
require(marelac)
#require(AquaEnv)

# non-public packages
require(diagenesis)

#=============================================================================
# Definition of specific functions
#=============================================================================

#-----------------------------------------------------------------------------
# Saturation function
#-----------------------------------------------------------------------------

FSAT <- function(C,K,n) (C/K)^n/((C/K)^n + 1)

#-----------------------------------------------------------------------------
# Function: tortuosity 
#-----------------------------------------------------------------------------

tortuosity <- function(por) 1-2*log(por)

#-----------------------------------------------------------------------------
# Function: init.diffusion.coefficient 
#-----------------------------------------------------------------------------

init.diffusion.coefficient <- function (species,S,TC,P,conv.fac,grid,tort.grid)
{
  Dmol <- conv.fac*diffcoeff(S=S,t=TC,P=P,species=species)[[species]]
  D.grid <- setup.prop.1D(value=Dmol,grid=grid)
  D.grid$mid <- D.grid$mid/tort.grid$mid
  D.grid$int <- D.grid$int/tort.grid$int
  return(D.grid)
}

#=============================================================================
# Definition of objects: species, reactions, elements, simulations
#=============================================================================

# SL = list of species incorporated in the model 
SL <- list()

# chemical species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="OC.f",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="OC.s",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="OC.r",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="MnO2",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="MnCO3",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Fe",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Mn",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.NH4",phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="HCO3",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="NH4",phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="O2",phase="solute",type="chemical")
#SL <- initialise.list.species(SL,"1D",name="N2",phase="solute",type="chemical")
#SL <- initialise.list.species(SL,"1D",name="NO3",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="SO4",phase="solute",type="chemical")

#SL <- initialise.list.species(SL,"1D",name="Ca",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Mn",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="HS",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="CH4",phase="solute",type="chemical")

#SL <- initialise.list.species(SL,"1D",name="H",phase="solute",type="chemical")
#SL <- initialise.list.species(SL,"1D",name="H2O",phase="solute",type="chemical")

# Composite species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="OC",phase="solid",type="composite")

# Artificial species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="Omega.FeS",phase="solute",type="artificial")
SL <- initialise.list.species(SL,"1D",name="Omega.MnCO3",phase="solute",type="artificial")

# Initialise species matrix 

spec.mat <- initialise.species.matrix(SL)
spec.mat["OC",c("OC.f","OC.s","OC.r")] <- c(1,1,1)

spec.table <- as.data.frame(spec.mat)
remove(spec.mat)


# List of reactions included 

RL <- list()

#---------------------
# mineralization terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="Cmin.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.r",type="kinetic")

RL <- initialise.list.reaction(RL,"1D",name="Cmin.ftot",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.stot",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.rtot",type="kinetic")

#---------------------
# individual mineralization terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="AR.f",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="DN.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MR.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FR.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG.f",type="kinetic")

RL <- initialise.list.reaction(RL,"1D",name="AR.s",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="DN.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MR.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FR.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG.s",type="kinetic")

RL <- initialise.list.reaction(RL,"1D",name="AR.r",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="DN.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MR.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FR.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG.r",type="kinetic")

#---------------------
# aerobic oxidation terms
# excluding reactions with sulfide minerals
#---------------------

#RL <- initialise.list.reaction(RL,"1D",name="AAO",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="SAO",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="MnO",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="SMnO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="AMO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSO",type="kinetic")

#---------------------
# anaerobic oxdiation terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="FMR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSFO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SMO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FMO",type="kinetic")

#---------------------
# sulfide mineral reactions
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="ISP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISD",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="ISO",type="kinetic")

#---------------------
# carbonate mineral reactions
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="CCP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CCD",type="kinetic")

#---------------------
# Iron, Manganese and Ammonium adsorption
#---------------------

#RL <- initialise.list.reaction(RL,"1D",name="AMS",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="FIS",type="kinetic")
#RL <- initialise.list.reaction(RL,"1D",name="MnS",type="kinetic")


# List of elements included 

EL <- list() # create empty element object list

EL <- initialise.list.element(EL,"1D",name="C")
EL <- initialise.list.element(EL,"1D",name="O")
EL <- initialise.list.element(EL,"1D",name="H")
EL <- initialise.list.element(EL,"1D",name="S")
#EL <- initialise.list.element(EL,"1D",name="Ca")
EL <- initialise.list.element(EL,"1D",name="Fe")
EL <- initialise.list.element(EL,"1D",name="Mn")
EL <- initialise.list.element(EL,"1D",name="N")
EL <- initialise.list.element(EL,"1D",name="e")

# Initialise element matrix 

elt.mat <- initialise.element.matrix(SL, EL)

elt.mat["OC.f",c("C","H","O","N")] <- c(1,25/11,1,1/11)
elt.mat["OC.s",c("C","H","O","N")] <- c(1,25/11,1,1/11)
elt.mat["OC.r",c("C","H","O","N")] <- c(1,25/11,1,1/11)

elt.mat["MnO2",c("Mn","O")] <- c(1,2)
elt.mat["FeOOH",c("Fe","O","H")] <- c(1,2,1)
elt.mat["FeS",c("S","Fe")] <- c(1,1)
elt.mat["MnCO3",c("Mn","C","O")] <- c(1,1,3)

elt.mat["HCO3",c("H","C","O","e")] <- c(1,1,3,1)
elt.mat["NH4",c("N","H","e")] <- c(1,4,-1)
elt.mat["X.NH4",c("N","H","e")] <- c(1,4,-1)

elt.mat["O2",c("O")] <- c(2)
#elt.mat["N2",c("N")] <- c(2)
#elt.mat["NO3",c("N","O","e")] <- c(1,3,1)
elt.mat["SO4",c("S","O","e")] <- c(1,4,2)

#elt.mat["Ca",c("Ca","e")] <- c(1,-2)
elt.mat["Mn",c("Mn","e")] <- c(1,-2)
elt.mat["X.Mn",c("Mn","e")] <- c(1,-2)
elt.mat["Fe",c("Fe","e")] <- c(1,-2)
elt.mat["X.Fe",c("Fe","e")] <- c(1,-2)
elt.mat["HS",c("H","S","e")] <- c(1,1,1)
elt.mat["CH4",c("C","H")] <- c(1,4)

#elt.mat["H",c("H","e")] <- c(1,-1)
#elt.mat["H2O",c("H","O")] <- c(2,1)

elt.table <- as.data.frame(elt.mat)
remove(elt.mat)

################################################################################
# initialise.parameters function
# This function initialises a default set of parameters (baseline simulation),
# but this function is also used when updating the parameter set.
################################################################################

initialise.parameters <- function(PL=NULL)
  
  # Model parameters are intialised on the first call (PL = NULL), and updated 
  # when a parmeter list is supplied (PL not NULL)
  
  # begin initialise.parameters
{
  
  #=============================================================================
  # Physical units 
  #=============================================================================
  
  if (is.null(PL$units$M)) PL$units$M <- "umol" # Mass unit
  if (is.null(PL$units$T)) PL$units$T <- "yr"   # Time unit
  if (is.null(PL$units$L)) PL$units$L <- "cm"   # Length unit
  
  # Limiting concentration: numerical limit on consumption in kinetic reactions
  PL$C.lim <- 1E-3
  
  #=============================================================================
  # Model domain and grid definition
  #=============================================================================
  
  if (is.null (PL$L)) PL$L <- 30   # depth of sediment domain [cm]
  if (is.null (PL$N)) PL$N <- 400  # number of grid layers
  
  # Call "setup.grid.1D" from ReacTran
  if (is.null (PL$grid.type)) PL$grid.type <- "uniform"  # number of grid layers
  if (is.null (PL$dx.1)) PL$dx.1 <- PL$L/(10*PL$N)
  #if (PL$grid.type == "uniform")
  #{
    # even grid - all depth layers have same thickness 
    #PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N)
  #} else {
    # uneven grid - higher resolution near the sediment water interface
    PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N,dx.1=PL$dx.1)
  #}
  
  #=============================================================================
  # Sediment parameters: 
  #=============================================================================
  
  # Environmental parameters
  
  if (is.null(PL$S))  PL$S  <- 33    # salinity
  if (is.null(PL$TC)) PL$TC <- 20    # temperature [deg C]
  if (is.null(PL$P))  PL$P  <- 1.013 # pressure [bar]
  
  # Porosity profile 
  
  if (is.null(PL$por.0)) PL$por.0 <- 0.75     # porosity at the sediment-water interface
    
  PL$por.grid <- setup.prop.1D(value=PL$por.0,grid=PL$grid)
  PL$svf.grid <- setup.prop.1D(value=(1-PL$por.0),grid=PL$grid)
  
  
  # Initialisation tortuosity 

  PL$tort.grid <- setup.prop.1D(value=0,grid=PL$grid)
  PL$tort.grid$mid <- tortuosity(PL$por.grid$mid)
  PL$tort.grid$int <- tortuosity(PL$por.grid$int)

  # Fixed profiles: pH
  
  PL$pH      <- 7.5
  PL$pH.grid <-  setup.prop.1D(value=PL$pH,grid=PL$grid)
  
  #=============================================================================
  # Diffusion coefficients 
  # Uses routine 'diffcoeff' from 'marelac' package
  #=============================================================================
  
  c.fac <- 10000*(3600*24*365.25) # conversion from [m2 s-1] to [cm2 yr-1]
  
  # Diffusion coefficients from marelac package
  
  PL$D.HCO3.grid <- init.diffusion.coefficient(species="HCO3",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.NH4.grid <- init.diffusion.coefficient(species="NH4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.O2.grid <- init.diffusion.coefficient(species="O2",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  #PL$D.N2.grid <- init.diffusion.coefficient(species="N2",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  #PL$D.NO3.grid <- init.diffusion.coefficient(species="NO3",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.SO4.grid <- init.diffusion.coefficient(species="SO4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  #PL$D.Ca.grid <- init.diffusion.coefficient(species="Ca",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.Mn.grid <- init.diffusion.coefficient(species="Mn",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.Fe.grid <- init.diffusion.coefficient(species="Fe",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.HS.grid <- init.diffusion.coefficient(species="HS",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.CH4.grid <- init.diffusion.coefficient(species="CH4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  #PL$D.H2O.grid <- init.diffusion.coefficient(species="H2O",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  #PL$D.H.grid <- init.diffusion.coefficient(species="H",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  #=============================================================================
  # Transport parameters
  #=============================================================================
  
  # Diffusive boundary layer
  
  if (is.null (PL$x.DBL)) PL$x.DBL <- 0 # thickness of DBL [cm]
  
  # Advective velocities 
  
  if (is.null(PL$rho.sed)) PL$rho.sed <- 2.6  # density solid sediment [g cm-3]
  if (is.null(PL$u.0)) PL$u.0 <- 0  # sedimentation velocity pore water [cm yr-1]
  if (is.null(PL$v.0)) PL$v.0 <- 0  # sedimentation velocity solids [cm yr-1]
  
  PL$SedFlux <- PL$rho.sed*(1-PL$por.0)*PL$v.0  # sedimentation flux [g cm-2 yr-1]
  
  PL$v.grid <- setup.compaction.1D(v.0=PL$v.0,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v # advection velocity solid
  PL$u.grid <- setup.compaction.1D(v.0=PL$v.0,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u # advection velocity pore water
  
  # Bioturbation profile
  
  if (is.null(PL$x.L))  PL$x.L    <- 0
  if (is.null(PL$x.att))PL$x.att  <- 0
  if (is.null(PL$y.inf))PL$y.inf  <- 0
  if (is.null(PL$Db.0)) PL$Db.0   <- 0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
  PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
  
  # Irrigation profile
  
  if (is.null(PL$irr.0)) PL$irr.0 <- 0     # irrigation rate at SWI [yr-1]
  PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=0,x.L=0,x.att=3,grid=PL$grid)
  
  #=============================================================================
  # Reaction parameters
  #=============================================================================
  
  # Decay constants organic matter
  
  if (is.null(PL$k.fast))      PL$k.f     <- 10     # decay constant of organic matter [yr-1]
  if (is.null(PL$k.s))         PL$k.s     <- 0.1    # decay constant of organic matter [yr-1]
  if (is.null(PL$k.r))         PL$k.r     <- 0      # decay constant of organic matter [yr-1] (Meysman et al., 2003)
     
  if (is.null(PL$K_O2))   PL$K_O2     <- 0.001         # Monod constant O2 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  #if (is.null(PL$K_NO3))  PL$K_NO3    <- 0.001         # Monod constant NO3 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_MnO2)) PL$K_MnO2   <- 1*PL$rho.sed  # Monod constant MnIV reduction [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_FeOOH))PL$K_FeOOH <- 4*PL$rho.sed  # Monod constant FeIII reduction [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_SO4))  PL$K_SO4   <- 0.9           # Monod constant SO4 reduction [umol cm-3 or mM] (Meysman et al., 2003)
  
  if (is.null(PL$R_CN))   PL$R_CN    <- 5        # from data
  #if (is.null(PL$R_MC))   PL$R_MC    <- 50        # from data
  
  # aerobic oxidation reactions
  # excluding sulfide mineral reactions
  
  #if (is.null(PL$k.SAO))  PL$k.SAO  <- 0.5      # kinetic constant ammonium oxidation [umol-1 cm3 yr-1] (Meysman et al., 2003)
  #if (is.null(PL$k.AAO))  PL$k.AAO  <- 0.5      # kinetic constant ammonium oxidation [umol-1 cm3 yr-1] (Meysman et al., 2003)
  #if (is.null(PL$k.MnO))  PL$k.MnO  <- 3.4E-4     # kinetic constant manganese oxidation [umol-1 cm3 yr-1] (Berg et al., 2003)
  #if (is.null(PL$k.SMnO))  PL$k.SMnO  <- 3.4E-4     # kinetic constant manganese oxidation [umol-1 cm3 yr-1] (Berg et al., 2003)
  
  if (is.null(PL$k.FIO)) PL$k.FIO <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.SIO)) PL$k.SIO <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.CSO)) PL$k.CSO  <- 1E+07    # kinetic constant sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.AMO)) PL$k.AMO <- 1E+04  # kinetic constant aerobic methane oxidation [µmol-1 cm3 yr-1] (Meysman et al., 2003)

  # anaerobic oxidation reactions
  
  if (is.null(PL$k.FMR)) PL$k.FMR <- 2.4E-6     # kinetic constant iron oxidation by MnIV [umol-1 cm3 yr-1] (Berg et al., 2003)
  if (is.null(PL$k.CSFO)) PL$k.CSFO <- 494     # kinetic constant sulfide oxidation by fresh FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.SMO)) PL$k.SMO <- 1E+04     # kinetic constant Sulfate methane oxidation [µmol-1 cm3 yr-1] (-)
  if (is.null(PL$k.FMO)) PL$k.FMO <- 1E+04     # kinetic constant FeOOH methane oxidation [µmol-1 cm3 yr-1] (-)
  
  # reactions involving sulfide minerals
  
  #if (is.null(PL$k.ISO))  PL$k.ISO  <- 1E+07    # kinetic constant iron sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996)
    
  if (is.null(PL$k.ISP)) PL$k.ISP <- 1E+04   # kinetic constant FeS precipitation [umol-1 cm3 yr-1] (Meysman et al., 2015)
  if (is.null(PL$k.ISD)) PL$k.ISD <- 3       # kinetic constant FeS dissolution [yr-1] (Meysman et al., 2015)
  if (is.null(PL$n.ISP)) PL$n.ISP <- 1       # kinetic exponent FeS precipitation [] (Meysman et al., 2003)
  if (is.null(PL$n.ISD)) PL$n.ISD <- 1       # kinetic exponent FeS dissolution [] (Meysman et al., 2003)
  
  if (is.null(PL$K_IS))  PL$K_IS  <- (10^-PL$pH)*1000*3.16    # Saturation constant [umol-1 cm3] (Rickard, 2006)
  
  # reactions involving carbonate minerals
  
  if (is.null(PL$k.CCP)) PL$k.CCP <- 1E-1   # kinetic constant CaCO3 precipitation [umol cm-3 yr-1] (Meysman et al., 2003)
  if (is.null(PL$k.CCD)) PL$k.CCD <- 0       # kinetic constant CaCO3 dissolution [yr-1] (Meysman et al., 2003)
  if (is.null(PL$n.CCP)) PL$n.CCP <- 1       # kinetic exponent CaCO3 precipitation [] (Meysman et al., 2003)
  if (is.null(PL$n.CCD)) PL$n.CCD <- 1       # kinetic exponent CaCO3 dissolution [] (Meysman et al., 2003)
  
  if (is.null(PL$K_CC))  PL$K_CC  <- 1e-06    # Saturation constant [umol-1 cm3] (Meysman, 2003)
  
  
  #if (is.null(PL$k.FIS)) PL$k.FIS <- 1    # kinetic constant Ferrous iron sorption  [yr-1] (Berg et al., 2003)
  if (is.null(PL$K.FIS)) PL$K.FIS <- PL$rho.sed*268   # Equilibrium constant Ferrous iron sorption [] (Berg et al., 2003)
  if (is.null(PL$K.AMS)) PL$K.AMS <- PL$rho.sed*1.8  # Equilibrium constant ammonium sorption [] (Berg et al., 2003; Mackin and Aller, 1984)
  if (is.null(PL$K.MnS)) PL$K.MnS <- PL$rho.sed*13   # Equilibrium constant Manganese sorption [] (Berg et al., 2003)
  
  #if (is.null(PL$k.FIS)) PL$k.FIS <- 1   # Kinetic constant Ferrous iron sorption [] (Berg et al., 2003)??
  #if (is.null(PL$k.AMS)) PL$k.AMS <- 1   # Kinetic constant ammonium sorption [] (Berg et al., 2003)??
  #if (is.null(PL$k.MnS)) PL$k.MnS <- 1   # Kinetic constant Manganese sorption [] (Berg et al., 2003)??
  
  #=============================================================================
  # Boundary conditions
  #=============================================================================
  
  # Upper boundary fluxes

  F.OC.f <- 0 # mmol m-2 d-1
  if (is.null(PL$F.OC.f)) PL$F.OC.f <- F.OC.f*365.25/10 # conversion to umol cm-2 yr-1 
  
  #F.OC.s <- 10 # mmol m-2 d-1
  #if (is.null(PL$F.OC.s)) PL$F.OC.s <- F.OC.s*365.25/10 # conversion to umol cm-2 yr-1 
  
  #F.OC.r <- 0 # mmol m-2 d-1
  #if (is.null(PL$F.OC.r)) PL$F.OC.r <- F.OC.r*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.FeS <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS)) PL$F.FeS <- F.FeS*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.MnCO3 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.MnCO3)) PL$F.MnCO3 <- F.MnCO3*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.MnO2 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.MnO2)) PL$F.MnO2 <- F.MnO2*365.25/10 # conversion to umol cm-2 yr-1 
  
  # 1% Fe content of mineral deposition 
  F.FeOOH <- 8/36.525#1E06*0.01*PL$SedFlux/(55.847) # mmol m-2 d-1
  if (is.null(PL$F.FeOOH)) PL$F.FeOOH <- F.FeOOH*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.X.Mn <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.Mn)) PL$F.X.Mn <- F.X.Mn*365.25/10 # conversion to umol cm-2 yr-1 
  F.X.Fe <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.Fe)) PL$F.X.Fe <- F.X.Fe*365.25/10 # conversion to umol cm-2 yr-1 
  F.X.NH4 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.NH4)) PL$F.X.NH4 <- F.X.NH4*365.25/10 # conversion to umol cm-2 yr-1 
  
  # Upper boundary concentrations

  if (is.null(PL$HCO3.ow)) PL$HCO3.ow <- 2.2 # SumCO2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$NH4.ow))  PL$NH4.ow  <- 0
  
  if (is.null(PL$O2.ow))  PL$O2.ow   <- 0.28 # O2 concentration bottom water [umol cm-3 or mM]
  #if (is.null(PL$N2.ow))  PL$N2.ow   <- 0
  #if (is.null(PL$NO3.ow)) PL$NO3.ow  <- 0
  if (is.null(PL$SO4.ow)) PL$SO4.ow  <- 140.0*PL$S/((32.065+4*15.999)*1.80655)
  
  #if (is.null(PL$Ca.ow))  PL$Ca.ow   <- 8
  if (is.null(PL$Mn.ow))  PL$Mn.ow   <- 0
  if (is.null(PL$Fe.ow))  PL$Fe.ow   <- 0
  if (is.null(PL$HS.ow))  PL$HS.ow   <- 0   # SumH2S concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$CH4.ow)) PL$CH4.ow  <- 0
  
  #if (is.null(PL$H2O.ow)) PL$H2O.ow   <- 0
  #if (is.null(PL$H.ow))   PL$H.ow     <- 0 # Proton concentration bottom water, without equilibrium with H2O
    
  # Boundary condition lower boundary 
  # All species: no gradient
  
  if (is.null(PL$HCO3.ds)) PL$HCO3.ds <- NA
  if (is.null(PL$NH4.ds))  PL$NH4.ds  <- NA
  
  if (is.null(PL$O2.ds)) PL$O2.ds   <- NA
  #if (is.null(PL$N2.ds)) PL$N2.ds   <- NA
  #if (is.null(PL$NO3.ds)) PL$NO3.ds <- NA
  if (is.null(PL$SO4.ds)) PL$SO4.ds <- NA
    
  #if (is.null(PL$Ca.ds)) PL$Ca.ds   <- NA
  if (is.null(PL$Mn.ds)) PL$Mn.ds   <- NA
  if (is.null(PL$Fe.ds)) PL$Fe.ds   <- NA
  if (is.null(PL$HS.ds)) PL$HS.ds   <- NA
  if (is.null(PL$CH4.ds)) PL$CH4.ds <- NA
  
  #if (is.null(PL$H2O.ds)) PL$H2O.ds <- NA
  #if (is.null(PL$H.ds)) PL$H.ds     <- NA
  
  #=============================================================================
  # Flags 
  #=============================================================================
  
  if (is.null(PL$simulation.type)) PL$simulation.type <- "direct.steady.state"
  #simulation.type <- "dynamic.steady.state"
  #simulation.type <- "time.dependent"
  
  ##############################################################################
  return(PL)
  # end initialise.parameters
}
##############################################################################

#=============================================================================
# Model formulation
#=============================================================================

CSFe.model <- function (t,state,parameters,summary.call=FALSE) 
{
  
  with(as.list(c(parameters)),{
    
    #---------------------------------------------------------------------------
    # Initialisation of 
    # SV : matrix of state variables
    # OL : Object List
    #---------------------------------------------------------------------------
    
    SV <- matrix(nrow=N,ncol=N.var,data=state)
    for (i in 1:N.var) OL[[var.names[i]]]$C <- SV[,i]
    
    #---------------------------------------------------------------------------
    OL <- within(OL,{
    #---------------------------------------------------------------------------
      
    #=========================================================================
    # Initialisation of concentration depth profiles
    #=========================================================================
    #-------------------------------------------------------------------------
    # Initialisation of species with fixed concentration depth profiles 
    #-------------------------------------------------------------------------
    if (is.null(PL$OC.f.fix)) OC.f.fix <- rep(0,grid$N)
    if (is.null(PL$OC.s.fix)) OC.s.fix <- rep(0,grid$N)
    if (is.null(PL$OC.r.fix)) OC.r.fix <- rep(0,grid$N)
    
    if (!("OC.f" %in% var.names)) OC.f$C <- OC.f.fix
    if (!("OC.s" %in% var.names)) OC.s$C <- OC.s.fix
    if (!("OC.r" %in% var.names)) OC.r$C <- OC.r.fix
    
    X.Fe$C  <- K.FIS*Fe$C
    X.NH4$C <- K.AMS*NH4$C
    X.Mn$C  <- K.MnS*Mn$C
    
    #-------------------------------------------------------------------------
    # Upper and lower boundary conditions
    #-------------------------------------------------------------------------

    HCO3$C.up <- HCO3.ow
    NH4$C.up  <- NH4.ow
    
    O2$C.up   <- O2.ow
    #N2$C.up   <- N2.ow
    #NO3$C.up  <- NO3.ow
    SO4$C.up  <- SO4.ow
    
    #Ca$C.up   <- Ca.ow
    Fe$C.up   <- Fe.ow
    Mn$C.up   <- Mn.ow
    HS$C.up   <- HS.ow
    CH4$C.up  <- CH4.ow
    
    #H2O$C.up  <- H2O.ow
    #H$C.up    <- H.ow
    
    if (is.na(HCO3.ds)) HCO3$C.down <- HCO3$C[N] else HCO3$C.down <- HCO3.ds 
    if (is.na(NH4.ds))  NH4$C.down <- NH4$C[N] else NH4$C.down <- NH4.ds 
    
    if (is.na(O2.ds))   O2$C.down  <- O2$C[N] else O2$C.down <- O2.ds 
    #if (is.na(N2.ds))   N2$C.down  <- N2$C[N] else N2$C.down <- N2.ds 
    #if (is.na(NO3.ds))  NO3$C.down <- NO3$C[N] else NO3$C.down <- NO3.ds 
    if (is.na(SO4.ds))  SO4$C.down <- SO4$C[N] else SO4$C.down <- SO4.ds 
    
    #if (is.na(Ca.ds))   Ca$C.down  <- Ca$C[N] else Ca$C.down <- Ca.ds 
    if (is.na(Mn.ds))   Mn$C.down  <- Mn$C[N] else Mn$C.down <- Mn.ds 
    if (is.na(Fe.ds))   Fe$C.down  <- Fe$C[N] else Fe$C.down <- Fe.ds 
    if (is.na(HS.ds))   HS$C.down  <- HS$C[N] else HS$C.down <- HS.ds 
    if (is.na(CH4.ds))  CH4$C.down <- CH4$C[N] else CH4$C.down <- CH4.ds 
    
    #if (is.na(H2O.ds))  H2O$C.down <- H2O$C[N] else H2O$C.down <- H2O.ds 
    #if (is.na(H.ds))    H$C.down   <- H$C[N] else H$C.down <- H.ds 
    
    
    #=========================================================================
    # TRANSPORT terms using tran.1D from ReacTran
    #=========================================================================
    
    if (summary.call) {ind <- c("tran","C.up","C.down","dif.flux","adv.flux","flux","flux.up","flux.down")} else {ind <- c("tran","flux.up","flux.down")}
    
    # Transport terms: Solids [umol cm-3 solid yr-1]
    
    OC.f[ind]  <- tran.1D(C=OC.f$C,flux.up=F.OC.f,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    #OC.s[ind]  <- tran.1D(C=OC.s$C,flux.up=F.OC.s,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    #OC.r[ind]  <- tran.1D(C=OC.r$C,flux.up=F.OC.r,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    MnO2[ind]  <- tran.1D(C=MnO2$C,flux.up=F.MnO2,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH[ind]<- tran.1D(C=FeOOH$C,flux.up=F.FeOOH,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS[ind]   <- tran.1D(C=FeS$C,flux.up=F.FeS,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    MnCO3[ind]  <- tran.1D(C=MnCO3$C,flux.up=F.MnCO3,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    X.Fe[ind]  <- tran.1D(C=X.Fe$C,flux.up=F.X.Fe,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    X.NH4[ind] <- tran.1D(C=X.NH4$C,flux.up=F.X.NH4,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    X.Mn[ind]  <- tran.1D(C=X.Mn$C,flux.up=F.X.Mn,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    # Transport terms: Solutes [umol cm-3 pore water yr-1]
    
    HCO3[ind] <- tran.1D(C=HCO3$C,C.up=HCO3$C.up,C.down=HCO3$C.down,v=u.grid,D=D.HCO3.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    NH4[ind]  <- tran.1D(C=NH4$C,C.up=NH4$C.up,C.down=NH4$C.down,v=u.grid,D=D.NH4.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    O2[ind]  <- tran.1D(C=O2$C,C.up=O2$C.up,C.down=O2$C.down,v=u.grid,D=D.O2.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    #N2[ind]  <- tran.1D(C=N2$C,C.up=N2$C.up,C.down=N2$C.down,v=u.grid,D=D.N2.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    #NO3[ind] <- tran.1D(C=NO3$C,C.up=NO3$C.up,C.down=NO3$C.down,v=u.grid,D=D.NO3.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    SO4[ind] <- tran.1D(C=SO4$C,C.up=SO4$C.up,C.down=SO4$C.down,v=u.grid,D=D.SO4.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    #Ca[ind]  <- tran.1D(C=Ca$C,C.up=Ca$C.up,C.down=Ca$C.down,v=u.grid,D=D.Fe.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    Mn[ind]  <- tran.1D(C=Mn$C,C.up=Mn$C.up,C.down=Mn$C.down,v=u.grid,D=D.Fe.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    Fe[ind]  <- tran.1D(C=Fe$C,C.up=Fe$C.up,C.down=Fe$C.down,v=u.grid,D=D.Fe.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    HS[ind]  <- tran.1D(C=HS$C,C.up=HS$C.up,C.down=HS$C.down,v=u.grid,D=D.HS.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    CH4[ind] <- tran.1D(C=CH4$C,C.up=CH4$C.up,C.down=CH4$C.down,v=u.grid,D=D.CH4.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    #H2O[ind] <- tran.1D(C=H2O$C,C.up=H2O$C.up,C.down=H2O$C.down,v=u.grid,D=D.H2O.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    #H[ind]   <- tran.1D(C=H$C,C.up=H$C.up,C.down=H$C.down,v=u.grid,D=D.H.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
        
    #-----------------------------------------------------------------------------
    # Saturation state calculation
    #-----------------------------------------------------------------------------
    
    # Saturation state FeS 
    
    Omega.FeS$C      <- Fe$C*HS$C/(K_IS)#*FSAT(Fe$C,C.lim/100,5)*FSAT(HS$C,C.lim/100,5)
    Omega.FeS$C.up   <- Fe.ow*HS$C.up/K_IS
    Omega.FeS$C.down <- Fe$C.down*HS$C.down/K_IS
    
    # Saturation state CaCO3 
    
    Omega.MnCO3$C      <- Mn$C*HCO3$C/(K_CC)#*FSAT(Fe$C,C.lim/100,5)*FSAT(HS$C,C.lim/100,5)
    Omega.MnCO3$C.up   <- Mn.ow*HCO3$C.up/K_CC
    Omega.MnCO3$C.down <- Mn$C.down*HCO3$C.down/K_CC
    
    #=========================================================================
    # Rates of Kinetic reactions 
    # Expressed per volume of bulk sediment [umol cm-3 yr-1]
    #=========================================================================
    
    # Mineralisation reactions
    
    Cmin.f$R <- svf.grid$mid*k.f*OC.f$C*FSAT(OC.f$C,C.lim,5)
    Cmin.s$R <- svf.grid$mid*k.s*OC.s$C*FSAT(OC.s$C,C.lim,5)
    Cmin.r$R <- svf.grid$mid*k.r*OC.r$C*FSAT(OC.r$C,C.lim,5)

    O2.lim <- O2$C/(O2$C+K_O2)*FSAT(O2$C,C.lim,5) #*(O2$C>sqrt(.Machine$double.eps))
    O2.inh <- K_O2/(O2$C+K_O2)
    MnO2.lim <- MnO2$C/(MnO2$C+K_MnO2)*FSAT(MnO2$C,C.lim,5)
    MnO2.inh <- K_MnO2/(MnO2$C+K_MnO2)
    FeOOH.lim <- FeOOH$C/(FeOOH$C+K_FeOOH)*FSAT(FeOOH$C,C.lim,5)
    FeOOH.inh <- K_FeOOH/(FeOOH$C+K_FeOOH)
    SO4.lim <- SO4$C/(SO4$C+K_SO4)*FSAT(SO4$C,C.lim,5)
    SO4.inh <- K_SO4/(SO4$C+K_SO4)
    
    O2.f <- O2.lim
    Mn.f <- MnO2.lim*O2.inh
    Fe.f <- FeOOH.lim*MnO2.inh*O2.inh
    S.f  <- SO4.lim*MnO2.inh*O2.inh
    CH4.f <- SO4.inh*FeOOH.inh*MnO2.inh*O2.inh
    sum.inh <- O2.f + Mn.f + Fe.f +S.f + CH4.f
    
     AR.f$R <- O2.f/sum.inh*Cmin.f$R  # Aerobic respiration
     MR.f$R <- Mn.f/sum.inh*Cmin.f$R # Denitrifiction
     FR.f$R <- Fe.f/sum.inh*Cmin.f$R # FeIII reduction
     SR.f$R <- S.f/sum.inh*Cmin.f$R # Sulfate reduction
     MG.f$R <- CH4.f/sum.inh*Cmin.f$R # methanogenesis

     AR.s$R <- O2.f/sum.inh*Cmin.s$R  # Aerobic respiration
     MR.s$R <- Mn.f/sum.inh*Cmin.s$R # Denitrifiction
     FR.s$R <- Fe.f/sum.inh*Cmin.s$R # FeIII reduction
     SR.s$R <- S.f/sum.inh*Cmin.s$R # Sulfate reduction
     MG.s$R <- CH4.f/sum.inh*Cmin.s$R # methanogesis

     AR.r$R <- O2.f/sum.inh*Cmin.r$R  # Aerobic respiration
     MR.r$R <- Mn.f/sum.inh*Cmin.r$R # Denitrifiction
     FR.r$R <- Fe.f/sum.inh*Cmin.r$R # FeIII reduction
     SR.r$R <- S.f/sum.inh*Cmin.r$R # Sulfate reduction
     MG.r$R <- CH4.f/sum.inh*Cmin.r$R # methanogesis

    # Aerobic oxidation reactions
    
    #AAO$R  <- por.grid$mid*k.AAO*O2$C*FSAT(O2$C,C.lim/100,5)*NH4$C*FSAT(NH4$C,C.lim/100,5) # ammonium oxidation
    #SAO$R  <- svf.grid$mid*k.SAO*O2$C*FSAT(O2$C,C.lim/100,5)*X.NH4$C*FSAT(X.NH4$C,C.lim/100,5) # Sorbed ammonium oxidation
    #MnO$R  <- por.grid$mid*k.MnO*O2$C*FSAT(O2$C,C.lim/100,5)*Mn$C*FSAT(Mn$C,C.lim/100,5) # manganese oxidation
    #SMnO$R  <- svf.grid$mid*k.SMnO*O2$C*FSAT(O2$C,C.lim/100,5)*X.Mn$C*FSAT(X.Mn$C,C.lim/100,5) # Sorbed manganese oxidation
    FIO$R  <- por.grid$mid*k.FIO*O2$C*FSAT(O2$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5) # iron re-oxidation
    SIO$R  <- svf.grid$mid*k.SIO*O2$C*FSAT(O2$C,C.lim,5)*X.Fe$C*FSAT(X.Fe$C,C.lim,5) # sorbed iron re-oxidation
    CSO$R  <- por.grid$mid*k.CSO*O2$C*FSAT(O2$C,C.lim/100,5)*HS$C*FSAT(HS$C,C.lim/100,5) # Sulfide re-oxidation
    AMO$R  <- por.grid$mid*k.AMO*O2$C*FSAT(O2$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # aerobic methane oxidation
 
    # Anaerobic oxidation reactions
    
    FMR$R <- svf.grid$mid*k.FMR*MnO2$C*FSAT(MnO2$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5)# iron re-oxidation by MnIV
    CSFO$R <- svf.grid$mid*k.CSFO*FeOOH$C*FSAT(FeOOH$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    SMO$R  <- por.grid$mid*k.SMO*SO4$C*FSAT(SO4$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # Sulfate methane oxidation
    FMO$R  <- svf.grid$mid*k.FMO*FeOOH$C*FSAT(FeOOH$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # iron oxide methane oxidation

    # Reactions involving sulfide and carbonate minerals
   
    #ISO$R  <- svf.grid$mid*k.ISO*O2$C*FSAT(O2$C,C.lim,5)*FeS$C*FSAT(FeS$C,C.lim,5) # iron sulfide re-oxidation
    ISP$R  <- svf.grid$mid*k.ISP*((Omega.FeS$C-1)^n.ISP)*(Omega.FeS$C>1)               #iron sulfide precipitation
    #ISD$R  <- svf.grid$mid*k.ISD*FeS$C*FSAT(FeS$C,C.lim,5)*(1-Omega.FeS$C)^n.ISD*(Omega.FeS$C<1) #iron sulfide dissolution
    CCP$R  <- svf.grid$mid*k.CCP*((Omega.MnCO3$C-1)^n.CCP)*(Omega.MnCO3$C>1)               #calcium carbonate precipitation
    CCD$R  <- svf.grid$mid*k.CCD*MnCO3$C*FSAT(MnCO3$C,C.lim,5)*(1-Omega.MnCO3$C)^n.CCD*(Omega.MnCO3$C<1) #calcium carbonate dissolution

    # Reactions involving mineral sorption

    #AMS$R  <- svf.grid$mid*k.AMS*(K.AMS*NH4$C - X.NH4$C)*(NH4$C>0)*(X.NH4$C>0) # Ammonium sorption
    #FIS$R  <- svf.grid$mid*k.FIS*(K.FIS*Fe$C - X.Fe$C)*(Fe$C>0)*(X.Fe$C>0)     # Ferric Iron sorption
    #MnS$R  <- svf.grid$mid*k.MnS*(K.MnS*Mn$C - X.Mn$C)*(Mn$C>0)*(X.Mn$C>0)     # Manganese sorption

    #=========================================================================
    # Consumption terms: kinetic reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================
    
    OC.f$kin.reac  <- - (AR.f$R + MR.f$R+ FR.f$R + SR.f$R + MG.f$R)/svf.grid$mid
    #OC.s$kin.reac  <- - (AR.s$R + DN.s$R + MR.s$R+ FR.s$R + SR.s$R + MG.s$R)/svf.grid$mid
    #OC.r$kin.reac  <- - (AR.r$R + DN.r$R + MR.r$R+ FR.r$R + SR.r$R + MG.r$R)/svf.grid$mid
    MnO2$kin.reac <- (- 2*MR.f$R - 2*MR.s$R - 2*MR.r$R - FMR$R)/svf.grid$mid  
    FeOOH$kin.reac <- (- 4*FR.f$R - 4*FR.s$R - 4*FR.r$R + 4*FIO$R + 4*SIO$R  - 8*CSFO$R - 8*FMO$R)/svf.grid$mid #+ ISO$R   
    #FeS$kin.reac   <- (- ISO$R)/svf.grid$mid
    X.Fe$kin.reac  <- (- 4*SIO$R)/svf.grid$mid
    #X.NH4$kin.reac <- (- SAO$R)/svf.grid$mid
    #X.Mn$kin.reac  <- (- 2*SMnO$R)/svf.grid$mid
   
    HCO3$kin.reac<- (AR.f$R + MR.f$R+ FR.f$R + SR.f$R + AR.s$R + MR.s$R+ FR.s$R + SR.s$R + AR.r$R + MR.r$R+ FR.r$R + SR.r$R + 1/2*MG.f$R + 1/2*MG.s$R + 1/2*MG.r$R + AMO$R + SMO$R + FMO$R)/por.grid$mid
    NH4$kin.reac <-  (+ 1/R_CN*AR.f$R + 1/R_CN*MR.f$R+ 1/R_CN*FR.f$R + 1/R_CN*SR.f$R + 1/R_CN*MG.f$R + 1/R_CN*AR.s$R + 1/R_CN*MR.s$R+ 1/R_CN*FR.s$R + 1/R_CN*SR.s$R + 1/R_CN*MG.s$R + 1/R_CN*AR.r$R + 1/R_CN*MR.r$R+ 1/R_CN*FR.r$R + 1/R_CN*SR.r$R + 1/R_CN*MG.r$R)/por.grid$mid

    O2$kin.reac  <- (- AR.f$R - AR.s$R - AR.r$R - FIO$R - SIO$R - 2*CSO$R  - 2*AMO$R)/por.grid$mid #- 9/4*ISO$R
    #N2$kin.reac  <-  (+ 2/5*DN.f$R + 2/5*DN.s$R + 2/5*DN.r$R )/por.grid$mid
    #NO3$kin.reac <-  (- 4/5*DN.f$R - 4/5*DN.s$R - 4/5*DN.r$R + AAO$R + SAO$R)/por.grid$mid
    SO4$kin.reac <- (- 0.5*SR.f$R - 0.5*SR.s$R - 0.5*SR.r$R+ CSFO$R + CSO$R  - SMO$R)/por.grid$mid #+ ISO$R

    Mn$kin.reac  <- (2*MR.f$R + 2*MR.s$R + 2*MR.r$R + FMR$R)/por.grid$mid

    Fe$kin.reac  <- (4*FR.f$R + 4*FR.s$R + 4*FR.r$R - 4*FIO$R  + 8*CSFO$R - 2*FMR$R + 8*FMO$R)/por.grid$mid
    HS$kin.reac  <- (0.5*SR.f$R + 0.5*SR.s$R + 0.5*SR.r$R - CSFO$R - CSO$R + SMO$R)/por.grid$mid
    CH4$kin.reac <- (+ 1/2*MG.f$R + 1/2*MG.s$R + 1/2*MG.r$R - AMO$R - SMO$R - FMO$R)/por.grid$mid
    
    #H2O$kin.reac <- (+ 2/5*DN.f$R + 2/5*DN.s$R + 2/5*DN.r$R + 2*MR.f$R + 2*MR.s$R + 2*MR.r$R + 6*FR.f$R + 6*FR.s$R + 6*FR.r$R + 1/2*MG.f$R + 1/2*MG.s$R + 1/2*MG.r$R + AAO$R + SAO$R - 2*MnO$R - 2*SMnO$R - 2*FMR$R - 6*FIO$R - 6*SIO$R + 12*CSFO$R  + 2*AMO$R + 2*SMO$R + 14*FMO$R)/por.grid$mid #- 3/2*ISO$R
    #H$kin.reac   <- (+ (R_CN - 1)/R_CN*AR.f$R + (R_CN - 1)/R_CN*AR.s$R + (R_CN - 1)/R_CN*AR.r$R + (R_CN-5)/(5*R_CN)*DN.f$R + (R_CN-5)/(5*R_CN)*DN.s$R + (R_CN-5)/(5*R_CN)*DN.r$R - (3*R_CN + 1)/R_CN*MR.f$R - (3*R_CN + 1)/R_CN*MR.s$R - (3*R_CN + 1)/R_CN*MR.r$R - (7*R_CN + 1)/R_CN*FR.f$R - (7*R_CN + 1)/R_CN*FR.s$R - (7*R_CN + 1)/R_CN*FR.r$R + (R_CN-2)/(2*R_CN)*SR.f$R + (R_CN-2)/(2*R_CN)*SR.s$R + (R_CN-2)/(2*R_CN)*SR.r$R - (3*R_CN + 2)/(2*R_CN)*MG.f$R - (3*R_CN + 2)/(2*R_CN)*MG.s$R - (3*R_CN + 2)/(2*R_CN)*MG.r$R + 2*AAO$R + 2*SAO$R + 4*MnO$R + 8*SMnO$R + 2*FMR$R + 8*FIO$R + 8*SIO$R + CSO$R - 15*CSFO$R  - SMO$R - 16*FMO$R)/por.grid$mid #+ 2*ISO$R

    
    #=========================================================================
    # Equilibrium reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    FeS$eq.reac  <- (ISP$R - ISD$R)/svf.grid$mid
    #X.Fe$eq.reac <- (FIS$R)/svf.grid$mid
    #Fe$eq.reac   <- (- ISP$R + ISD$R - FIS$R)/por.grid$mid
    Fe$eq.reac   <- (- ISP$R + ISD$R)/por.grid$mid
    HS$eq.reac   <- (- ISP$R + ISD$R)/por.grid$mid
    #H$eq.reac   <- (ISP$R - ISD$R - CCP$R + CCD$R)/por.grid$mid

    MnCO3$eq.reac  <- (CCP$R - CCD$R)/svf.grid$mid
    HCO3$eq.reac  <- (- CCP$R + CCD$R)/por.grid$mid
    #Mn$eq.reac  <- (- CCP$R + CCD$R - MnS$R)/por.grid$mid
    Mn$eq.reac  <- (- CCP$R + CCD$R)/por.grid$mid
    #X.Mn$eq.reac  <- (MnS$R)/svf.grid$mid

    #NH4$eq.reac  <- (- AMS$R)/por.grid$mid
    #X.NH4$eq.reac  <- (AMS$R)/svf.grid$mid


    #=========================================================================
    # Irrigation 
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    HCO3$irr <- irr.grid$mid*(HCO3.ow - HCO3$C)
    NH4$irr  <- irr.grid$mid*(NH4.ow - NH4$C)    

    O2$irr   <- irr.grid$mid*(O2.ow - O2$C)
    #N2$irr   <- irr.grid$mid*(N2.ow - N2$C)    
    #NO3$irr  <- irr.grid$mid*(NO3.ow - NO3$C)    
    SO4$irr  <- irr.grid$mid*(SO4.ow - SO4$C)

    Fe$irr   <- irr.grid$mid*(Fe.ow - Fe$C)
    HS$irr   <- irr.grid$mid*(HS.ow - HS$C)
    CH4$irr  <- irr.grid$mid*(SO4.ow - SO4$C)
    
    #=========================================================================
    # Rate of changes composite species 
    #=========================================================================
    
    #OC$tran     <- OC.f$tran + OC.s$tran + OC.r$tran
    #OC$kin.reac <- OC.f$kin.reac + OC.s$kin.reac + OC.r$kin.reac
    #OC$irr      <- OC.f$irr + OC.s$irr + OC.r$irr
    #OC$eq.reac  <- OC.f$eq.reac + OC.s$eq.reac + OC.r$eq.reac

    #=========================================================================
    # Rate of changes of all species
    #=========================================================================
    
  # Fe$irr <- alfa*(PL$Fe.ow - Fe$C)

    OC.f$ddt  <- OC.f$tran  + OC.f$kin.reac  + OC.f$eq.reac  + OC.f$irr
    #OC.s$ddt  <- OC.s$tran  + OC.s$kin.reac  + OC.s$eq.reac  + OC.s$irr
    #OC.r$ddt  <- OC.r$tran  + OC.r$kin.reac  + OC.r$eq.reac  + OC.r$irr
    MnO2$ddt  <- MnO2$tran + MnO2$kin.reac + MnO2$eq.reac + MnO2$irr
    FeOOH$ddt <- FeOOH$tran + FeOOH$kin.reac + FeOOH$eq.reac + FeOOH$irr
    FeS$ddt   <- FeS$tran   + FeS$kin.reac   + FeS$eq.reac   + FeS$irr
    MnCO3$ddt <- MnCO3$tran   + MnCO3$kin.reac   + MnCO3$eq.reac   + MnCO3$irr

    HCO3$ddt  <- HCO3$tran  + HCO3$kin.reac  + HCO3$eq.reac  + HCO3$irr
    NH4$ddt   <- 1/(por.grid$mid+svf.grid$mid*K.AMS)*(por.grid$mid*(NH4$tran   + NH4$kin.reac   + NH4$eq.reac   + NH4$irr) + svf.grid$mid*(X.NH4$tran   + X.NH4$kin.reac   + X.NH4$eq.reac   + X.NH4$irr))
    #NH4$ddt   <- NH4$tran   + NH4$kin.reac   + NH4$eq.reac   + NH4$irr
    #X.NH4$ddt <- X.NH4$tran   + X.NH4$kin.reac   + X.NH4$eq.reac   + X.NH4$irr


    O2$ddt    <- O2$tran    + O2$kin.reac    + O2$eq.reac    + O2$irr
    #N2$ddt    <- N2$tran    + N2$kin.reac    + N2$eq.reac    + N2$irr
    #NO3$ddt   <- NO3$tran   + NO3$kin.reac   + NO3$eq.reac   + NO3$irr
    SO4$ddt   <- SO4$tran   + SO4$kin.reac   + SO4$eq.reac   + SO4$irr

    #Ca$ddt    <- Ca$tran    + Ca$kin.reac    + Ca$eq.reac    + Ca$irr

    Mn$ddt    <- 1/(por.grid$mid+svf.grid$mid*K.MnS)*(por.grid$mid*(Mn$tran + Mn$kin.reac + Mn$eq.reac + Mn$irr) + svf.grid$mid*(X.Mn$tran + X.Mn$kin.reac + X.Mn$eq.reac + X.Mn$irr))
    #Mn$ddt    <- Mn$tran + Mn$kin.reac + Mn$eq.reac + Mn$irr
    #X.Mn$ddt  <- X.Mn$tran + X.Mn$kin.reac + X.Mn$eq.reac + X.Mn$irr
    Fe$ddt    <- 1/(por.grid$mid+svf.grid$mid*K.FIS)*(por.grid$mid*(Fe$tran + Fe$kin.reac + Fe$eq.reac + Fe$irr) + svf.grid$mid*(X.Fe$tran + X.Fe$kin.reac + X.Fe$eq.reac + X.Fe$irr))
    #X.Fe$ddt  <- X.Fe$tran + X.Fe$kin.reac + X.Fe$eq.reac + X.Fe$irr
    #Fe$ddt    <- Fe$tran + Fe$kin.reac + Fe$eq.reac + Fe$irr
    HS$ddt    <- HS$tran    + HS$kin.reac    + HS$eq.reac    + HS$irr
    CH4$ddt   <- CH4$tran    + CH4$kin.reac    + CH4$eq.reac    + CH4$irr
    
    #H2O$ddt   <- H2O$tran + H2O$kin.reac + H2O$eq.reac + H2O$irr
    #H$ddt     <- H$tran + H$kin.reac + H$eq.reac + H$irr
    
    OC$ddt <- OC$tran + OC$kin.reac + OC$eq.reac + OC$irr
    
    #=========================================================================
    }) # end within call
  #=========================================================================
  
  #-----------------------------------------------------------------------------
  # Assemble matrix with rate of changes (RC) 
  #-----------------------------------------------------------------------------
  
  RC <- 0*SV
  for (i in 1:N.var) RC[,i] <- OL[[var.names[i]]]$ddt
  
  #-----------------------------------------------------------------------------
  # Monitoring of simulation progession during dynamic simulation 
  #-----------------------------------------------------------------------------
  
  if ((simulation.type %in% c("dynamic.steady.state","time.dependent")) && (t > sim.info$tc)){
    
    # sim.info$fc = keeps track of number of function calls
    # sim.info$tc = keeps track of incremental simulation time
    
    sim.info$fc <<- sim.info$fc + 1  
    sim.info$tc <<- t*1
    
    # Calculate for each state variable 
    #  ACC = absolute concentration change 
    #  RCC = relative concentration change 
    
    SV.ref <- rep(1E-03,N.var)
    ACC.t <- ACC(RC)
    RCC.t <- RCC(RC,SV,SV.ref)
    
    # Output of simulation progress to the screen 
    screenprint.simulation.progress(sim.info$fc,t,SV,ACC.t,RCC.t,var.names)
    
    # Store simulation progress results 
    sim.info$progress <<- rbind(sim.info$progress,c(t,ACC.t,RCC.t))
    points(log10(t),log10(sqrt(sum(RCC.t^2)/length(RCC.t))),pch=16,col="red")
    
  }
  
  #===========================================================================
  # return statement 
  #===========================================================================
  
  if (!summary.call) {return(list(RC))}
  else {
    
    # Initialisation of lists
    
    SL <- OL[names(SL)]
    RL <- OL[names(RL)]
    EL <- EL
    
    # Updating species information
    
    SL <- updating.chemical.species(SL,PL) 
    SL <- updating.composite.species(spec.table,SL) 
    
    # Updating reaction information
    
    for (i in 1:length(RL)) RL[[i]]$R.int <- sum(RL[[i]]$R*grid$dx) 
    
    # Updating element information
    
    EL <- updating.elements(elt.table,EL,SL) 
    
    return(list(SL=SL,RL=RL,EL=EL,PL=parameters))
    
   }
   
   #-----------------------------------------------------------------------------
 }) # end with call
}
#-----------------------------------------------------------------------------

################################################################################
# Initialisation of state variable matrix
################################################################################

initialise.state.variables <- function(PL)
{
  
  PL$var.names <- c("OC.f",#"OC.s","OC.r",
                    "MnO2","FeOOH","FeS","MnCO3","HCO3","NH4",#"X.NH4",
                    "O2",#"N2","NO3",
                    "SO4","Mn",#"X.Mn",
                    "Fe",#"X.Fe","Ca",
                    "HS","CH4")#,"H2O","H")
  PL$non.var.names <- NULL
  
  PL$N <- PL$grid$N
  
  PL$N.var <- length(PL$var.names)
  PL$N.non.var <- length(PL$non.var.names)
  
  return(PL)
}

initialise.species.arrays <- function(PL,SL)
{
  
  for (i in 1:length(SL))
  {
    SL[[i]]$C <- rep(0,PL$N)
    SL[[i]]$tran <- rep(0,PL$N)
    SL[[i]]$kin.reac <- rep(0,PL$N)
    SL[[i]]$eq.reac <- rep(0,PL$N)
    SL[[i]]$irr <- rep(0,PL$N)
  }
  return(SL)
}

# Function that creates a rough guess for the initial state

crude.initial.state <- function(PL)
{
  with(PL,{
    
    # Initialisation State Variable matrix (SV) 
    
    SV.guess <- matrix(nrow=N,ncol=N.var,dimnames=list(1:N,var.names))
    
    if ("OC.f" %in% PL$var.names) SV.guess[,"OC.f"] <- rep(0,length.out=N)
    #if ("OC.s" %in% PL$var.names) SV.guess[,"OC.s"] <- rep(0,length.out=N)
    #if ("OC.r" %in% PL$var.names) SV.guess[,"OC.r"] <- rep(0,length.out=N)
    
    if ("MnO2" %in% PL$var.names)  SV.guess[,"MnO2"] <- rep(0,length.out=N)
    if ("FeOOH" %in% PL$var.names)  SV.guess[,"FeOOH"] <- rep(0,length.out=N)
    if ("FeS" %in% PL$var.names)  SV.guess[,"FeS"] <- rep(0,length.out=N)
    if ("MnCO3" %in% PL$var.names)  SV.guess[,"MnCO3"] <- rep(0,length.out=N)
    
    if ("HCO3" %in% PL$var.names) SV.guess[,"HCO3"] <- rep(HCO3.ow,length.out=N)
    if ("NH4" %in% PL$var.names) SV.guess[,"NH4"] <- rep(NH4.ow,length.out=N)
    #if ("X.NH4" %in% PL$var.names) SV.guess[,"X.NH4"] <- rep(0,length.out=N)
    
    if ("O2" %in% PL$var.names) SV.guess[,"O2"] <- rep(O2.ow,length.out=N)
    #if ("N2" %in% PL$var.names) SV.guess[,"N2"] <- rep(N2.ow,length.out=N)
    #if ("NO3" %in% PL$var.names) SV.guess[,"NO3"] <- rep(NO3.ow,length.out=N)
    if ("SO4" %in% PL$var.names) SV.guess[,"SO4"] <- rep(SO4.ow,length.out=N)
    
    if ("Mn" %in% PL$var.names) SV.guess[,"Mn"] <- rep(Mn.ow,length.out=N)
    #if ("X.Mn" %in% PL$var.names) SV.guess[,"X.Mn"] <- rep(0,length.out=N)
    #if ("Ca" %in% PL$var.names) SV.guess[,"Ca"] <- rep(Ca.ow,length.out=N)
    if ("Fe" %in% PL$var.names) SV.guess[,"Fe"] <- rep(Fe.ow,length.out=N)
    #if ("X.Fe" %in% PL$var.names) SV.guess[,"X.Fe"] <- rep(0,length.out=N)
    if ("HS" %in% PL$var.names) SV.guess[,"HS"] <- rep(HS.ow,length.out=N)
    if ("CH4" %in% PL$var.names) SV.guess[,"CH4"] <- rep(CH4.ow,length.out=N)
    
    #if ("H2O" %in% PL$var.names) SV.guess[,"H2O"] <- rep(H2O.ow,length.out=N)
    #if ("H" %in% PL$var.names) SV.guess[,"H"] <- rep(H.ow,length.out=N)
    
    return(SV.guess)
    
  })
}


################################################################################
# Test trial model runs
################################################################################

model <- CSFe.model
sim.info <- list(tc= 0,fc=0)

# Initialisation of object list (OL)  

PL <- initialise.parameters(PL=list())
PL$simulation.type <- "direct.steady.state"
SL <- initialise.species.arrays(PL,SL) # default action
OL <- c(SL,RL) # default action

# Test trial 

PL <- initialise.state.variables(PL)
SV <- crude.initial.state(PL)

system.time({
  result <- model(t=0,state=SV,parameters=PL,summary.call=FALSE)
  result <- model(t=0,state=SV,parameters=PL,summary.call=TRUE)
})

remove(result)
remove(SV)

