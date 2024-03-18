################################
# Code to produce trajectories of the approximated process in (13) using Gillespie's algorithm
# and riproduce results described in Section 3.3 of 
# Sabbioni E., Szabò R., Siri P., Cappelletti D., Lente G., Bibbona E. 
# "Final nanoparticle size distribution under unusual parameter regimes", 
# arXiv...

# When the initial number of monomers is bigger than 10^8, we suggest not to run 
# the following code on your local machine, but to use a cluster where it can be 
# run in parallel.
################################
using Pkg

pathProject = "nanoparticlesProject"                                        # path of the project
pathOutput =  "nanoparticlesProject/code/output/ApproxGillespie"            # path where the results will be saved
Pkg.activate(pathProject)

### LOAD JULIA PACKAGES
using Revise, Distributions, Random, LinearAlgebra, PDMats, StatsBase
using RCall
using npformation
using DelimitedFiles
using BenchmarkTools

### SET THE VALUES OF PARAMETERS FOR THE ASSOCIATED EXACT PROCESS (7)
growth::Float64 = 5.0                           # growth rate constant \gamma 
mon::Int64 = 1e+04                              # initial number of monomers N 
mins::Int32 = 3                                 # number of monomers used during nucleation n
nuc::Float64 = Float64(mon)^(1-mins)*growth     # nucleation rate constant \nu, derived such that the new alternative scaling holds (see Section 2.1.2)
repetition::Int64 = 6                           # number of simulated trajectories
threshold::Int32 = ceil(Int32, 5 * sqrt(mon))   # maximum typical size of nanoparticles in the original process

### DERIVE THE VALUES OF PARAMETERS FOR THE APPROXIMATED PROCESS (13)
a::Int64 = ceil(Int64, mon^(1/4))                  # bin amplitude
b::Int64 = ceil(Int64, mon^(3/4))                  # number of bin
thresholdBin::Int64 = ceil(Int64, 5 * sqrt(mon)/a) # maximum typical size of nanoparticles in the approximated process
 
growthBin = growth                                 # growth rate constant of the approximated process
nucBin = nuc                                       # nucleation rate constant of the approximated process

rngseed = 123456;                                  # seed
Random.seed!(rngseed);

### RUN GILLESPIE'S ALGORITHM TO SIMULATE THE TRAJECTORIES OF (13)
timeGillespieApprox = @elapsed resGillApprox = gillespieAPPROX(nucBin, growthBin, mon, b, mins, repetition; thresholdBin = thresholdBin, monNucleation = Int32(1),parNucleation = Int16(1), monGrowth = Int32(1), parGrowth = Int16(1))


### SAVE THE RESULTS 
# Send the results and the parameters to R
@rput resGillApprox;
@rput timeGillespieApprox;

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

write.csv(resGillApprox, file = paste(pathOutput, "/approxGill_",chrMon, ".csv", sep = ""))
write.csv(timeGillespieApprox, file = paste(pathOutput, "/approxGill_",chrMon, "_TIME.csv", sep = ""))
""" 