###############################################################################
# Simulation file 
# Transient diagenesis model 
# Author: Sebastiaan van de Velde
###############################################################################

require(signal)
load("measured data.Rdata")

# Source file containing model function
#source("Transient Model St130_disturbance_instantadsorption.R")
source("Transient Model St130_disturbance_instantadsorption_OC not dynamic simpl.R")

# Source file containing plotting info function
source("Plotting info 'Transient Model St130_disturbance'.R")

# Number of simulations

model <- CSFe.model

#=============================================================================
# Dynamic simulation (start from steady state)
# Find recovery of DIC and NH4
# disturb + follow evolution over time
#=============================================================================

load("00 Transient find start steady state try 01 3.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "01 Transient disturbance try 01"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$v.0 <- 0
PL$u.0 <- 0

PL <- initialise.parameters(PL)

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

l.dist <- 20
g <- length(which(PL$grid$x.mid<l.dist))

SV[,"O2"][which(PL$grid$x.mid<l.dist)]   <- rep(PL$O2.ow,length.out=g)
SV[,"MnO2"][which(PL$grid$x.mid<l.dist)] <-  rep(5*PL$rho.sed,length.out=g)
SV[,"FeOOH"][which(PL$grid$x.mid<l.dist)]<-  rep(50*PL$rho.sed,length.out=g)
SV[,"Fe"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Fe.ow,length.out=g)
SV[,"Mn"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Mn.ow,length.out=g)
SV[,"NH4"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$NH4.ow,length.out=g)
SV[,"MnCO3"][which(PL$grid$x.mid<l.dist)]<-  rep(0,length.out=g)
SV[,"HCO3"][which(PL$grid$x.mid<l.dist)] <-  rep(PL$HCO3.ow,length.out=g)
SV[,"SO4"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$SO4.ow,length.out=g)
SV[,"HS"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$HS.ow,length.out=g)
SV[,"FeS"][which(PL$grid$x.mid<l.dist)]  <-  rep(0,length.out=g)
SV[,"H2O"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$H2O.ow,length.out=g)
SV[,"H"][which(PL$grid$x.mid<l.dist)]    <-  rep(PL$H.ow,length.out=g)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,2*24,24*30,24*59)/(24*365.25)#,3*24*30,4*24*30)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","2 days","30 days","59 days")#,"3 months","4 months")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation (start from steady state)
# Find recovery of DIC and NH4
# add OC.f after disturbance + follow evolution over time
#=============================================================================

load("00 Transient find start steady state try 01 3.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "01 Transient disturbance try 02"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL <- initialise.parameters(PL)

PL$v.0 <- 0
PL$u.0 <- 0

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

l.dist <- 20
g <- length(which(PL$grid$x.mid<l.dist))

SV[,"OC.f"][which(PL$grid$x.mid<l.dist)]   <- rep(100,length.out=g)
SV[,"O2"][which(PL$grid$x.mid<l.dist)]   <- rep(PL$O2.ow,length.out=g)
SV[,"MnO2"][which(PL$grid$x.mid<l.dist)] <-  rep(5*PL$rho.sed,length.out=g)
SV[,"FeOOH"][which(PL$grid$x.mid<l.dist)]<-  rep(50*PL$rho.sed,length.out=g)
SV[,"Fe"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Fe.ow,length.out=g)
SV[,"Mn"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Mn.ow,length.out=g)
SV[,"NH4"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$NH4.ow,length.out=g)
SV[,"X.NH4"][which(PL$grid$x.mid<l.dist)]   <-  rep(0,length.out=g)
SV[,"MnCO3"][which(PL$grid$x.mid<l.dist)]<-  rep(0,length.out=g)
SV[,"HCO3"][which(PL$grid$x.mid<l.dist)] <-  rep(PL$HCO3.ow,length.out=g)
SV[,"SO4"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$SO4.ow,length.out=g)
SV[,"HS"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$HS.ow,length.out=g)
SV[,"FeS"][which(PL$grid$x.mid<l.dist)]  <-  rep(0,length.out=g)
SV[,"H2O"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$H2O.ow,length.out=g)
SV[,"H"][which(PL$grid$x.mid<l.dist)]    <-  rep(PL$H.ow,length.out=g)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,2*24,7*24,14*24,24*30,24*45,24*59,24*90)/(24*365.25)#,24*120,3*24*30,4*24*30)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","2 days","7 days","14 days","30 days","45 days","59 days","90 days")#,"120 days")#,"3 months","4 months")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation (start from steady state)
# tune Mn dynamics
# add 'mixed region' between 15 and 20 cm
# add OC.f after disturbance + follow evolution over time
#=============================================================================

load("00 Transient find start steady state try 01 3.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "01 Transient disturbance try 02"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL <- initialise.parameters(PL)

#PL$F.OC.s <- 36.525
PL$k.CCP <- 5E-5
PL$k.ISP <- 1E+4
PL$K_CC  <- 1E-5
PL$K_MnO2 <- 1*PL$rho.sed
PL$K_FeOOH <- 4*PL$rho.sed 
PL$K.AMS <- 1.5   #4.32
PL$K.MnS <- 5   #31.2 according to Berg
PL$K.FIS <- 100 #643.2 according to Berg
  
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
g <- length(which(PL$grid$x.mid<l.dist))

SV[,"OC.f"][which(PL$grid$x.mid<l.dist)]   <- rep(55,length.out=g)
SV[,"O2"][which(PL$grid$x.mid<l.dist)]   <- rep(PL$O2.ow,length.out=g)
SV[,"MnO2"][which(PL$grid$x.mid<l.dist)] <-  rep(1*PL$rho.sed,length.out=g)
SV[,"FeOOH"][which(PL$grid$x.mid<l.dist)]<-  rep(30*PL$rho.sed,length.out=g)
SV[,"Fe"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Fe.ow,length.out=g)
SV[,"Mn"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Mn.ow,length.out=g)
SV[,"NH4"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$NH4.ow,length.out=g)
SV[,"MnCO3"][which(PL$grid$x.mid<l.dist)]<-  rep(0,length.out=g)
SV[,"HCO3"][which(PL$grid$x.mid<l.dist)] <-  rep(PL$HCO3.ow,length.out=g)
SV[,"SO4"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$SO4.ow,length.out=g)
SV[,"HS"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$HS.ow,length.out=g)
SV[,"FeS"][which(PL$grid$x.mid<l.dist)]  <-  rep(0,length.out=g)
SV[,"H2O"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$H2O.ow,length.out=g)
SV[,"H"][which(PL$grid$x.mid<l.dist)]    <-  rep(PL$H.ow,length.out=g)

#l.mix <- 20
#start <- max(which(PL$grid$x.mid<l.dist))
#end <- min(which(PL$grid$x.mid>l.mix))
#h <- length(PL$grid$x.mid[start:end])

#SV[,"OC.f"][start:end] <- approx(y=c(SV[,"OC.f"][start],SV[,"OC.f"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"O2"][start:end]   <- approx(y=c(SV[,"O2"][start],SV[,"O2"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"MnO2"][start:end] <- approx(y=c(SV[,"MnO2"][start],SV[,"MnO2"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"FeOOH"][start:end]<- approx(y=c(SV[,"FeOOH"][start],SV[,"FeOOH"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"Fe"][start:end]   <- approx(y=c(SV[,"Fe"][start],SV[,"Fe"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"Mn"][start:end]   <- approx(y=c(SV[,"Mn"][start],SV[,"Mn"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"NH4"][start:end]  <- approx(y=c(SV[,"NH4"][start],SV[,"NH4"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"MnCO3"][start:end]<- approx(y=c(SV[,"MnCO3"][start],SV[,"MnCO3"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"HCO3"][start:end] <- approx(y=c(SV[,"HCO3"][start],SV[,"HCO3"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"SO4"][start:end]  <- approx(y=c(SV[,"SO4"][start],SV[,"SO4"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"HS"][start:end]   <- approx(y=c(SV[,"HS"][start],SV[,"HS"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"FeS"][start:end]  <- approx(y=c(SV[,"FeS"][start],SV[,"FeS"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"H2O"][start:end]  <- approx(y=c(SV[,"H2O"][start],SV[,"H2O"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y
#SV[,"H"][start:end]    <- approx(y=c(SV[,"H"][start],SV[,"H"][end]), x = c(PL$grid$x.mid[start],PL$grid$x.mid[end]),xout=PL$grid$x.mid[start:end], method = "linear", n = h)$y

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,5*24,34*24,80*24)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","5 days","34 days","80 days")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation (start from steady state)
# Find recovery of DIC and NH4
# add OC.f after disturbance + follow evolution over time
#=============================================================================

load("00 Transient find start steady state try 01 3.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "01 Transient disturbance try 03"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$v.0 <- 0.1
PL$u.0 <- 0.1

PL <- initialise.parameters(PL)

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

l.dist <- 15
g <- length(which(PL$grid$x.mid<l.dist))

SV[,"OC.f"][which(PL$grid$x.mid<l.dist)]   <- rep(100,length.out=g)
SV[,"O2"][which(PL$grid$x.mid<l.dist)]   <- rep(PL$O2.ow,length.out=g)
SV[,"MnO2"][which(PL$grid$x.mid<l.dist)] <-  rep(5*PL$rho.sed,length.out=g)
SV[,"FeOOH"][which(PL$grid$x.mid<l.dist)]<-  rep(50*PL$rho.sed,length.out=g)
SV[,"Fe"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Fe.ow,length.out=g)
SV[,"Mn"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$Mn.ow,length.out=g)
SV[,"NH4"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$NH4.ow,length.out=g)
SV[,"X.NH4"][which(PL$grid$x.mid<l.dist)]   <-  rep(0,length.out=g)
SV[,"MnCO3"][which(PL$grid$x.mid<l.dist)]<-  rep(0,length.out=g)
SV[,"HCO3"][which(PL$grid$x.mid<l.dist)] <-  rep(PL$HCO3.ow,length.out=g)
SV[,"SO4"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$SO4.ow,length.out=g)
SV[,"HS"][which(PL$grid$x.mid<l.dist)]   <-  rep(PL$HS.ow,length.out=g)
SV[,"FeS"][which(PL$grid$x.mid<l.dist)]  <-  rep(0,length.out=g)
SV[,"H2O"][which(PL$grid$x.mid<l.dist)]  <-  rep(PL$H2O.ow,length.out=g)
SV[,"H"][which(PL$grid$x.mid<l.dist)]    <-  rep(PL$H.ow,length.out=g)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,2*24,7*24,14*24,24*30,24*45,24*59,24*90)/(24*365.25)#,24*120,3*24*30,4*24*30)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","2 days","7 days","14 days","30 days","45 days","59 days","90 days")#,"120 days")#,"3 months","4 months")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
