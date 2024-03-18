################################
# Code to produce trajectories of the exact process in (7) using Gillespie's algorithm, following the times at which the jumps occur and the number of available monomers at each jump time.

# When the initial number of monomers is bigger than 10^5, we suggest not to run 
# the following code on your local machine, but to use a cluster where it can be 
# run in parallel.
################################
using Pkg


pathProject = "nanoparticlesProject"                            # path of the project
pathOutput =  "nanoparticlesProject/code/output/traj"           # path where the results will be saved
Pkg.activate(pathProject)

### LOAD JULIA PACKAGES
using Revise, Distributions, Random, LinearAlgebra, PDMats, StatsBase
using RCall
using npformation
using DelimitedFiles
using BenchmarkTools


growth::Float64 = 5.0                               # growth rate constant \gamma 
mon::Int32 = 1e+05                                  # initial number of monomers N 
mins::Int32 = 3                                     # number of monomers used during nucleation n
nuc::Float64 = Float64(mon)^(1-mins)*growth
repetition::Int64 = 3                               # number of simulated trajectories
threshold::Int64 = ceil(Int64, 5 * sqrt(mon))       # maximum typical size of nanoparticles in the original process


rngseed = 123456;                                  # seed
Random.seed!(rngseed);

### RUN GILLESPIE'S ALGORITHM TO SIMULATE THE TRAJECTORIES OF (7)
timeGillespie_traj = @elapsed traj = gillespie_traj_mon(nuc, growth, mon, mins, repetition)


### SAVE THE RESULTS 
# Send the results and the parameters to R
@rput traj;
@rput timeGillespie_traj;
monR::Float64 = Float64(mon)
@rput monR;
@rput pathOutput;

R"""
if(monR == 1e+04){
    chrMon <-4
}else if(monR == 1e+05){
    chrMon <- 5
}else if(monR == 1e+06){
    chrMon <- 6
}else if(monR == 1e+07){
    chrMon <- 7
}else if(monR == 1e+08){
    chrMon <- 8
}else if(monR == 1e+09){
    chrMon <- 9
}else if(monR == 1e+10){
    chrMon <- 10
}else if(monR == 1e+11){
    chrMon <- 11
}else if(monR == 1e+12){
    chrMon <- 12
}


write.csv(traj, file = paste(pathOutput, "/traj_",chrMon, ".csv", sep = ""))
write.csv(timeGillespie_traj, file = paste(pathOutput, "/traj_",chrMon, "_TIME.csv", sep = ""))
""" 




