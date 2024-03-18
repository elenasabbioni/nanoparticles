####### Comparison trajetories of monomers in the exact process (7) and in the approximated process (13)

# Clear the workspace
rm(list=ls())

########################################
#           PACKAGES
########################################
# Load necessary libraries
require(data.table)


########################################
#           INITIALIZATION 
#           OF PARAMETERS
########################################
# Set initial parameters
mon <- 1e+5                             # initial number of monomers for process (7)
growth <- 5                             # growth rate constant 
mins <- 3                               # number of monomers used during nucleation n in the exact process (7)
nuc <- growth*(mon^(1-mins))            # nucleation rate constant 


chrMon <- format(mon, scientific = TRUE)
chrMon <- substr(chrMon, nchar(chrMon), nchar(chrMon))



########################################
#              LOAD THE RESULTS
#              OF THE SIMULATIONS
########################################
# Set the path where the simulation results are
pathInput <- "nanoparticlesProject/code/output"     # path where the results are stores


# load the trajectories and the jump times of the exact process (7)
traj <- fread(file=paste(pathInput, "/traj_", chrMon, ".csv", sep = ""))
 # fread(file=paste("C:/Users/elena/Dropbox (Politecnico Di Torino Studenti)/rebeka-turin/Code Enrico/output/traj/traj_", chrMon, ".csv", sep = "")) 

# load the trajectories and the jump times of the approximated process (13) 
timesAPP <- fread(file=paste(pathInput, "/times_trajAPPROX_", chrMon, ".csv", sep = ""))# fread(file=paste("C:/Users/elena/Dropbox (Politecnico Di Torino Studenti)/rebeka-turin/Code Enrico/output/traj/tempi_trajAPPROX_", chrMon, ".csv", sep =""))
stateAPP <- fread(file=paste(pathInput, "/mon_trajAPPROX_", chrMon, ".csv", sep ="")) # fread(file=paste("C:/Users/elena/Dropbox (Politecnico Di Torino Studenti)/rebeka-turin/Code Elena/output/traj/mon_trajAPPROX_", chrMon, ".csv", sep =""))



########################################
#              PLOTS
########################################
pathOutput <- "nanoparticlesProject/graphics"     # path where the plots will be stored


# --- exact process
traj <- as.matrix(traj)
k <- 1                                      # repetition that we analyze
times <- c(0, traj[k, 2:(mon+1)])           # extract the jumps'times 
states <- traj[k, (mon + 1):(2*mon)]        # extract the number of monomers at each jump step
rm(traj)

last.time <- min(which(times[-1] == 0))
last.state <- min(which(states[-1] == 0))
last <- max(last.time, last.state)
ff <- stepfun(times[2:last], states[1:last])                # interpolate the (times, # monomers) to obtain a continuos function describing monomers'evolution
x <- seq(0, times[last], length.out = 1000)

# --- approximatd process
timesAPP <- as.matrix(timesAPP)
stateAPP <- as.matrix(stateAPP)
k <- 1                                                      # repetition that we analyze
timesAPP <- timesAPP[k, ]                                   # extract the jumps'times 
statesAPP <- stateAPP[k,]                                   # extract the number of monomers at each jump step

last.timeAPP <- length(timesAPP)
last.stateAPP <-min(which(statesAPP == 0))
lastAPP <- max(last.timeAPP, last.stateAPP)
ffAPP <-stepfun(timesAPP[2:lastAPP], statesAPP[1:lastAPP])    # interpolate the (times, # monomers) to obtain a continuos function describing monomers'evolution
xAPP <- seq(0,timesAPP[lastAPP], length.out = 1000)


pdf(paste(pathOutput, "/cfrTraj_", chrMon, ".pdf", sep = ""))
plot(x = x, y = ff(x), type = "l", col = "red", main = paste("10^", chrMon, sep = ""))
lines(x = xAPP, y = ffAPP(xAPP), type="l", col="black")
dev.off()


summary(times)
summary(times2)















