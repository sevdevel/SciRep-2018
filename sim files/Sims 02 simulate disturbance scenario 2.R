# ##############################################################################
# Simulation file 
# Transient diagenesis model 
# Author: Sebastiaan van de Velde
# ##############################################################################

require(signal)
load("measured data.Rdata")

# Source file containing model function
source("Transient Model St130_disturbancescenario2_MSversion.R")

# Source file containing plotting info function
source("Plotting info 'Transient Model St130_disturbance'.R")

# Number of simulations

model <- CSFe.model

# =============================================================================
# Dynamic simulation (start from steady state)
# Find recovery of DIC and NH4
# disturb + follow evolution over time
# =============================================================================

load("00 Transient find start steady state try 01 3.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "02 Transient disturbance_scenario2_MSrun"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$F.OC.f <- 0
PL$k.CCP <- 5E-5
PL$k.ISP <- 1E+4
PL$K_CC  <- 1E-5
PL$K_MnO2 <- 1*PL$rho.sed
PL$K_FeOOH <- 4*PL$rho.sed 
PL$K.AMS <- 1.75   #4.32
PL$K.MnS <- 5   #31.2 according to Berg
PL$K.FIS <- 10 #643.2 according to Berg #100 
PL$k.FMO <- 1E+04
#PL$k.FIS <- 1

PL$u.0 <- 0
PL$v.0 <- 0

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

#SV[,"Fe"] <- approx(y=c(data.transientdiagenesis$Fe$may[[1]]$C,data.transientdiagenesis$Fe$may[[2]]$C)/1000,x=c(data.transientdiagenesis$Fe$may[[1]]$x,data.transientdiagenesis$Fe$may[[2]]$x)*100,xout=PL$grid$x.mid,rule=2)$y
#SV[,"Mn"] <- approx(y=c(data.transientdiagenesis$Mn$may[[1]]$C,data.transientdiagenesis$Mn$may[[2]]$C)/1000,x=c(data.transientdiagenesis$Fe$may[[1]]$x,data.transientdiagenesis$Fe$may[[2]]$x)*100,xout=PL$grid$x.mid,rule=2)$y

l.dist <- 15
#g <- length(which(PL$grid$x.mid<l.dist))
start <- max(which(PL$grid$x.mid<l.dist))
end <- length(PL$grid$x.mid)

SV[,"OC.f"][start:end]  <- approx(y=c(SV[,"OC.f"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"O2"][start:end]    <- approx(y=c(SV[,"O2"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"MnO2"][start:end]  <-  approx(y=c(SV[,"MnO2"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"FeOOH"][start:end] <-  approx(y=c(SV[,"FeOOH"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"Fe"][start:end]    <-  approx(y=c(data.transientdiagenesis$Fe$may[[1]]$C,data.transientdiagenesis$Fe$may[[2]]$C)/1000, x = c(data.transientdiagenesis$Fe$may[[1]]$x,data.transientdiagenesis$Fe$may[[2]]$x)*100,xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h, rule=2)$y
#SV[,"X.Fe"][start:end]  <-  PL$K.FIS*SV[,"Fe"][start:end]
SV[,"Mn"][start:end]    <-  approx(y=c(data.transientdiagenesis$Mn$may[[1]]$C,data.transientdiagenesis$Mn$may[[2]]$C)/1000, x = c(data.transientdiagenesis$Mn$may[[1]]$x,data.transientdiagenesis$Mn$may[[2]]$x)*100,xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h, rule=2)$y
#SV[,"X.Mn"][start:end]  <-  PL$K.MnS*SV[,"Mn"][start:end]
SV[,"NH4"][start:end]   <-  approx(y=c(SV[,"NH4"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
#SV[,"X.NH4"][start:end] <- PL$K.AMS*SV[,"NH4"][start:end]
SV[,"MnCO3"][start:end]<-  approx(y=c(SV[,"MnCO3"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"HCO3"][start:end] <-  approx(y=c(SV[,"HCO3"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"SO4"][start:end]  <-  approx(y=c(SV[,"SO4"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"HS"][start:end]   <-  approx(y=c(0,0), x = c(PL$grid$x.mid[1],PL$grid$x.mid[(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
SV[,"FeS"][start:end]  <-  approx(y=c(SV[,"FeS"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
#SV[,"H2O"][start:end]  <-  approx(y=c(SV[,"H2O"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y
#SV[,"H"][start:end]    <-  approx(y=c(SV[,"H"][1:(start-1)]), x = c(PL$grid$x.mid[1:(start-1)]),xout=(PL$grid$x.mid[start:end]-15), method = "linear", n = h,rule=2)$y

g <- length(PL$grid$x.mid[1:(start-1)])

SV[,"OC.f"][1:(start-1)]  <- rep(65,length.out=g)
SV[,"O2"][1:(start-1)]    <- rep(PL$O2.ow,length.out=g)
SV[,"MnO2"][1:(start-1)]  <-  rep(1*PL$rho.sed,length.out=g)
SV[,"FeOOH"][1:(start-1)] <-  rep(11*PL$rho.sed,length.out=g)
SV[,"Fe"][1:(start-1)]    <-  rep(PL$Fe.ow,length.out=g)
#SV[,"X.Fe"][1:(start-1)]  <-  PL$K.FIS*SV[,"Fe"][1:(start-1)]
SV[,"Mn"][1:(start-1)]    <-  rep(PL$Mn.ow,length.out=g)
#SV[,"X.Mn"][1:(start-1)]  <-  PL$K.MnS*SV[,"Mn"][1:(start-1)]
SV[,"NH4"][1:(start-1)]   <-  rep(PL$NH4.ow,length.out=g)
#SV[,"X.NH4"][1:(start-1)] <- PL$K.AMS*SV[,"NH4"][1:(start-1)]
SV[,"MnCO3"][1:(start-1)] <-  rep(0,length.out=g)
SV[,"HCO3"][1:(start-1)]  <-  rep(PL$HCO3.ow,length.out=g)
SV[,"SO4"][1:(start-1)]   <-  rep(PL$SO4.ow,length.out=g)
SV[,"HS"][1:(start-1)]    <-  rep(PL$HS.ow,length.out=g)
SV[,"FeS"][1:(start-1)]   <-  rep(0,length.out=g)
#SV[,"H2O"][1:(start-1)]  <-  rep(PL$H2O.ow,length.out=g)
#SV[,"H"][1:(start-1)]    <-  rep(PL$H.ow,length.out=g)


# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

# -------------------------------------------------------------------------------
# Sequence of time points where output is needed
# -------------------------------------------------------------------------------

sim.info$time.seq <- c(0,5,6,7,8,9,10,11,12,13,14,15,1*24,1.5*24,2*24,2.5*24,3*24,5*24,7*24,10*24,13*24,15*24,20*24,25*24,30*24,24*32,24*35,24*40,24*45,24*50,24*55,24*60,24*65,24*70,24*75,24*78)/(24*365.25)#,24*150,24*300)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","5 hr","6 hr","7 hr","8 hr","9 hr","10 hr","11 hr","12 hr","13 hr","14 hr","15 hr","1 day","1.5 days","2 days","2.5 days","3 days","5 days","7 days","10 days","13 days","15 days","20 days","25 days","30 days","32 days","35 days","40 days","45 days","50 days","55 days","60 days","65 days","70 days","75 days","78 days")#,"150 days","300 days")
sim.info$N.out <- length(sim.info$time.seq)

# ------------------------------------------------------------------------------
# Dynamic simulation
# -------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

# -------------------------------------------------------------------------------
# Plotting preparation
# -------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11(16,8)
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11(16,8)
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

# -------------------------------------------------------------------------------
# Save the simulation output list
# -------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

# ============================================================================
# Follow disturbance over longer time
# 
# 
# =============================================================================

load("02 Transient disturbance_scenario2_MSrun 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "02 Transient disturbance_scenario2_MSrun"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

# -------------------------------------------------------------------------------
# Sequence of time points where output is needed
# -------------------------------------------------------------------------------

sim.info$time.seq <- c(0,24*2,24*72,24*222,24*287.25,24*(287.25+365.25/2),24*(287.25+365.25))/(24*365.25)#,24*150,24*300)/(24*365.25) 
sim.info$time.schedule <- c("78 days","80 days","150 days","300 days","1 year","1.5 years","2 years")#,"150 days","300 days")
sim.info$N.out <- length(sim.info$time.seq)

# ------------------------------------------------------------------------------
# Dynamic simulation
# -------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

# -------------------------------------------------------------------------------
# Plotting preparation
# -------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11(16,8)
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11(16,8)
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

# -------------------------------------------------------------------------------
# Save the simulation output list
# -------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


