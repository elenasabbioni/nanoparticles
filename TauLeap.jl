################################
# Code to produce exact trajectories of (7) using Tau-Leap with Post-Leap checks algorithm
# and riproduce results described in Section 3.2 of 
# Sabbioni E., Szabò R., Siri P., Cappelletti D., Lente G., Bibbona E. 
# "Final nanoparticle size distribution under unusual parameter regimes", 
# arXiv...

# When the initial number of monomers is bigger than 10^5, we suggest not to run 
# the following code on your local machine, but to use a cluster where it can be 
# run in parallel.
################################

using Pkg

pathProject = "nanoparticlesProject"                                # path of the project
pathOutput =  "nanoparticlesProject/code/output/TauLeap"            # path where the results will be saved
Pkg.activate(pathProject)

### LOAD JULIA PACKAGES
using Revise, Distributions, Random, LinearAlgebra, PDMats, StatsBase
using RCall
using npformation
using DelimitedFiles
using BenchmarkTools

### SET THE VALUES OF THE PARAMETERS
growth::Float64 = 5.0                          # growth rate constant \gamma
mon::Int64 = 1e+04                             # initial number of monomers N
mins::Int64 = 3                                # number of monomers used during nucleation n
nuc::Float64 = Float64(mon)^(1-mins)*growth    # nucleation rate constant \nu, derived such that the new alternative scaling holds (see Section 2.1.2)
repetition::Int64 = 6                          # number of simulated trajectories
threshold::Int64 = ceil(Int64, 5 * sqrt(mon))  # maximum typical size of nanoparticles

### SET ALGORITHM-SPECIFIC PARAMETERS 
epsilon::Float64 = 0.7
p::Float64 = 0.5
pStar::Float64 = 0.9
q::Float64 = 0.9


rngseed = 123456;                              # seed
Random.seed!(rngseed);

### RUN TAU-LEAP'S ALGORITHM TO SIMULATE THE TRAJECTORIES OF (7)
timeTau = @elapsed resTau=TauLeap_PostLeap(nuc, growth, mon, mins, repetition, epsilon, p, pStar, q; threshold)

### SAVE THE RESULTS 
# Send the results and the parameters to R
@rput resTau;
@rput timeTau;

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


write.csv(resTau,  file = paste(pathOutput, "/tau-leap_",chrMon, ".csv", sep = ""))
write.csv(timeTau, file = paste(pathOutput, "/tau-leap_",chrMon, "TIME.csv", sep = ""))
""" 

