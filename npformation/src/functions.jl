# ------------------------------------------------------------------------
# Auxiliary functions necessary to run the simulation described in 
# Bibbona E., Cappelletti D., Siri P., Sabbioni E., Szabò R., Lente G. 
# "Final nanoparticle size distribution under unusual parameter regimes", 
# arXiv...
# ------------------------------------------------------------------------


####################################################
####################################################
#      GILLESPIE algorithm - EXACT PROCESS
####################################################
####################################################

# Function that simulates different trajectories of the process (7) of the paper using the Gillespie algorithm.
# The code returns the final size counts for each trajectory. 
function gillespie(
    nuc::Float64,               # nucleation rate constant \nu
    growth::Float64,            # growth rate constant \gamma
    mon::Int32,                 # initial number of monomers N
    mins::Int32,                # number of monomers used during nucleation n
    repetition::Int64;          # number of simulated trajectories
    threshold::Int64 = -10      # maximum typical size of nanoparticles
)


    parNuc::Int16 = 1           # number of nanoparticles created during nucleation
    parGrowth::Int16 = 1        # number of nanoparticles created during growth
    monGrowth::Int32 = 1        # number of monomerc used during growth

    sequenza = Int32.(0:(mins-1))

    # if not given as positive number as input, set the maximum typical size of nanoparticles to 5*N^0.5
    if threshold < 0
        threshold = ceil(Int64, 5 * sqrt(mon))
    end

    out::Matrix{Int64} = zeros(Int64, repetition, threshold)    # output matix (repetitions x threshold) with final numer of nanoparticles for each trajectories

    
    Threads.@threads for i = 1:repetition
        println("Repetition ", i)


        state = zeros(Int16, threshold);    # number of nanoparticles for each size: state[k+1] is the number of nanoparticles of size k
        stateMON::Int32 = mon;              # number of monomers
        largest::Int64 = mins + 1;          # largest partcle that can be created at the current time

        ratesGR = zeros(Int32, threshold);  # rates of growth of nanoparticles: rate[k+1] is the rate of growth of nanoparticles of size k (see eq. (5) and (6) of the paper)
        rateNUC::Float64 = 0;               # rate of nucleation (see eq. (4) of the paper)
        
        prob = zeros(Float64, threshold);   # probabilities of reactions to fire
        
        
        newpar::Float64 = nuc/growth;

        while stateMON > 0  # check if there are still available monomers
            
            # NOTE: all the rates are divided by the growth rate constant 
            # to increase the computational speed
            
            # update nucleation rate
            if stateMON >= mins     # nucleation can still fire
                rateNUC = newpar*exp(sum(log.(stateMON .- sequenza))); 
            else                    # nucleation can not fire
                rateNUC = 0.0;
            end 
            
            # update growth rate
            for i = mins:(largest-1)  
                ratesGR[i+1] = stateMON * state[i+1]; 
            end

            sumratesGR::Int64 = sum(ratesGR);
            sumrates = sumratesGR + rateNUC;
            
            prob[1] = exp(log(rateNUC) - log(sumrates));    # update probability that nucleation occurs

            # update probabilities of growth reactions
            for i = (mins+1):largest
                if ratesGR[i] > 0;
                    prob[i]= exp(log(ratesGR[i]) - log(sumrates)); 
                else
                    prob[i] = 0.0
                end
            end

            fires::Int64 = rand(Categorical(prob));     # index of the reaction that fires
            
            if fires == threshold
                error("Dimension of nanoparticles is exceding threshold")
            end

            if fires == 1 # nucleation occurs
                stateMON = stateMON - mins;                     # update number of monomers
                state[mins+1] = state[mins+1] + parNuc;         # update number of nanoparticles of minimal size
            else # growth of a particle of dimension "fires - 1" to a particle of dimension "fires"
                if fires == largest    
                    largest=largest+1; # update the maximum size of nanoparticle that we can obtain
                end
                stateMON = stateMON - monGrowth;                # update number of monomers
                state[fires] = state[fires] - parGrowth;        # update number of nanoparticles of dimension "fires - 1"
                state[fires+1] = state[fires+1] + parGrowth;    # update number of nanoparticles of dimension "fires"
            end
        end
        out[i, :] = Int64.(state);
        out[i, 1] = stateMON
    end

    return out
end


# Function that simulates different trajectories of the process (7) of the paper using the Gillespie algorithm and allows to follow the number of available monomers at each step of the algorithm.
# The code returns, for each trajectory, the times at which each reation occurs, the number of available monomers at each step of the algorithm and the total final number of nanoparticles. 
function gillespie_traj_mon(
    nuc::Float64,               # nucleation rate constant \nu
    growth::Float64,            # growth rate constant \gamma
    mon::Int32,                 # initial number of monomers N
    mins::Int32,                # number of monomers used during nucleation n
    repetition::Int64,          # number of simulated trajectories
)


    parNuc::Int16 = 1           # number of nanoparticles created during nucleation
    monGrowth::Int32 = 1        # number of monomerc used during growth

    sequenza = Int32.(0:(mins-1))

    out::Matrix{Float64} = zeros(Float64, repetition, 2*mon+1)  # output matix: each row store the results of a different repetition, columns 1:mon store the times at which reactions occur, columns (mon+1):(2*mon) store the number of available monomers at each step of the algorithm, column (2*mon + 1) store the total final number of nanoparticles


    Threads.@threads for i = 1:repetition
        println("Repetition ", i)


        statePAR::Int64 = 0;            # total number of nanoparticles (of each size)
        stateMON::Vector{Int64} = zeros(Int64, mon); # number of monomers at each step of the algorithm
        stateMON[1] = mon;       

        rateGR::Float64 = 0;            # rates of growth of nanoparticles (the nanoparticles of the different sizes are considered all together)
        rateNUC::Float64 = 0;           # rate of nucleation (see eq. (4) of the paper)
        
       
        prob = zeros(Float64, 2);       # probabilities of reactions to fire: the first element corresponds to the probability of nucleation, the second element to the probability of a growth

        times::Vector{Float64} = zeros(Float64, mon);  # times at which the reactions occur
        t::Float64 = 0;                 # time   
    
        
        newpar::Float64 = nuc/growth;
        count::Int64 = 1                # current step of the algorithm
        
        while stateMON[count] > 0       # check if there are still monomers available
            
            # NOTE: all the rates are divided by the growth rate constant 
            # to increase the computational speed

            # update nucleation rate
            if stateMON[count] >= mins  # nucleation can still fire
                rateNUC = newpar*exp(sum(log.(stateMON[count] .- sequenza)));
            else                        # nucleation can not fire
                rateNUC=0.0;
            end 

            # update growth rate
            rateGR = stateMON[count] * statePAR; 

            sumrates= rateGR+rateNUC;

            # sample dt such that t + dt is the time at which the next reaction fires
            dt = rand(Exponential(1/(growth*sumrates)));
            t = t + dt;
            times[count] = t
            
             # update probability that nucleation occurs
            if rateNUC > 0
                prob[1] = exp(log(rateNUC) - log(sumrates)); 
            else
                prob[1] = 0; 
            end
             # update probability that a growth occues
            prob[2] = 1 - prob[1];

            fires::Int64 = rand(Categorical(prob)); # index of the reaction that fires
            
            if fires == 1   # nucleation occurs
                stateMON[count+1] = stateMON[count] - mins;     # update number of monomers
                statePAR = statePAR + parNuc;                   # update number of nanoparticles
            else
                stateMON[count+1] = stateMON[count] - monGrowth;   # update number of monomers
            end
            
            count = count + 1
        end

        out[i, 1:mon] = times;
        out[i, (mon + 1):(2*mon)] = stateMON;
        out[i, 2*mon + 1] = statePAR;
    end

    return out
end




####################################################
####################################################
#      GILLESPIE algorithm - APPROXIMATED PROCESS
####################################################
####################################################

# Function that simulates different trajectories of the approximated process (13) of the paper using the Gillespie algorithm.
# The code returns the final size counts for each trajectory. 
function gillespieAPPROX(
    nucBin::Float64,                # nucleation rate constant 
    growthBin::Float64,             # growth rate constant 
    mon::Int64,                     # inital number of monomers N 
    b::Int64,                       # total numer of bins
    mins::Int32,                    # number of monomers used during nucleation n in the associated exact process (7)
    repetition::Int64;              # number of simulated trajectories
    thresholdBin::Int64 = -10,      # typical maximum number of filled bins
    monNucleation::Int32 = 1,       # number of a-mers consumed during nucleation
    parNucleation::Int16 = 1,       # number of nanoparticles removed from each bin during nucleation
    monGrowth::Int32 = 1,           # number of a-mers consumed during growth
    parGrowth::Int16 = 1            # number of nanoparticles removed from each bin during growth
)

    a::Int64 = ceil(Int64, mon^(1/4))  # bin amplitude

    # if not given as positive number as input, set the typical maximum number of filled bins to 5*N^0.5/a
    if thresholdBin < 0
        thresholdBin = ceil(Int64, 5 * sqrt(mon)/a)
    end

    out::Matrix{Int64} = zeros(Int64, repetition, thresholdBin)     # output matix (repetitions x thresholdBin) with final numer of nanoparticles in each bin for each trajectories


    Threads.@threads for i = 1:repetition
        println("Repetition ", i)


        state = zeros(Int16, thresholdBin);      # number of nanoparticles for each bin: state[k+1] is the number of nanoparticles of in the bin k
        stateMON::Int64 = ceil(Int64, mon/a);    # number of a-mers
        largest::Int64= Int64(2)                 # largest bin that can be filled at the current time

        ratesGR = zeros(Int64, thresholdBin);    # rates of growth of nanoparticles in the different bins: rate[k+1] is the rate of growth of nanoparticles in bin k (see eq. (15))
        rateNUC::Float64=0;                      # rate of nucleation (see eq. (14))
        prob = zeros(Float64, thresholdBin);     # probabilities of reactions to fire
        
    

        newpar::Float64 = nucBin/growthBin; 

        while stateMON >= monGrowth              # check if there are still available a-mers
            
            # NOTE: all the rates are divided by the growth rate constant growthBin to increase the computational speed

            # update nucleation rate
            rateNUC = newpar*exp(sum(log.(a*stateMON .-Int32.(0:(mins-1)))));
            
            # update growth rate
            for i = 1:(largest-1)  
                ratesGR[i+1] = stateMON #=a-mers=# *state[i+1] 
            end
            
            # update probabilities of reactions to fire
            test2 = @view prob[1:largest]
            prob[1] = log(rateNUC[1]);
            for i = 2:largest
                prob[i] = (log(ratesGR[i])); 
            end
            # due to computational approximization, use the methods describes in 
            # https://timvieira.github.io/blog/post/2014/07/31/gumbel-max-trick/
            # to sample from the multinomial distribution and to obtain the index of the reaction that fires
            fires = argmax(test2 .+ rand(Gumbel(0.0,1.0),largest))
            
            
            if fires == thresholdBin
                error("Dimension of nanoparticles is exceding thresholdBin")
            end

            if fires == 1                                       # nucleation occurs
                stateMON = stateMON - monNucleation;            # update number of a-mers
                state[2] = state[2] + parNucleation;            # update number of nanoparticles in the first bin
            else                                                # growth of a particle in bin "fires - 1" to particle in bin "fires"
                if fires == largest
                    largest = largest + 1;                      # update the maximum bin size in which we can create particles
                end
                stateMON = stateMON - monGrowth;                # update number of a-mers
                state[fires] = state[fires] - parGrowth;        # update number of nanoparticles in bin "fires - 1"
                state[fires+1] = state[fires+1] + parGrowth;    # update number of nanoparticles in bin "fires"
            end
        end

        out[i, :] .= Int64.(state);
        out[i, 1] = stateMON
    end
    return out
end



# Function that simulates different trajectories of the approximated process (13) of the paper using the Gillespie algorithm and allows to follow the number of available monomers at each step of the algorithm.
# The code returns, for each trajectory, the final size counts of nanoparticles in each bin, the times at which each reation occurs, the number of available a-mers at each step of the algorithm. 
function gillespieAPPROX_traj_mon(
    nucBin::Float64,                # nucleation rate constant 
    growthBin::Float64,             # growth rate constant 
    mon::Int64,                     # inital number of monomers N 
    b::Int64,                       # total numer of bins
    mins::Int32,                    # number of monomers used during nucleation n in the associated exact process (7)
    repetition::Int64;              # number of simulated trajectories
    thresholdBin::Int64 = -10,      # typical maximum number of filled bins
    monNucleation::Int32 = 1,       # number of a-mers consumed during nucleation
    parNucleation::Int16 = 1,       # number of nanoparticles removed from each bin during nucleation
    monGrowth::Int32 = 1,           # number of a-mers consumed during growth
    parGrowth::Int16                # number of nanoparticles removed from each 
)

    a::Int64 = ceil(Int64, mon^(1/4))                                                   # bin amplitude

    
    times::Matrix{Float64} = Matrix{Float64}(undef, repetition, Int64(ceil(mon/a)));    # times at which the reactions occur                    
    times[:, :] .= 0.0
    traj_mon::Matrix{Int64} = zeros(Int64, repetition, Int64(ceil(mon/a)));             # number of availbale a-mers at each step of the algorithm

    # if not given as positive number as input, set the typical maximum number of filled bins to 5*N^0.5/a
    if thresholdBin < 0
        thresholdBin = ceil(Int64, 5 * sqrt(mon)/a)
    end

    out::Matrix{Int64} = zeros(Int64, repetition, thresholdBin)  # output matix (repetitions x thresholdBin) with final numer of nanoparticles in each bin for each trajectories

    Threads.@threads for i = 1:repetition
        println("Repetition ", i)

        state = zeros(Int16, thresholdBin);         # number of nanoparticles for each bin: state[k+1] is the number of nanoparticles of in the bin k
        stateMON::Int64 = ceil(Int64, mon/a);       # number of a-mers
        largest::Int64= Int64(2)                    # largest bin that can be filled at the current time

        ratesGR = zeros(Int64, thresholdBin);       # rates of growth of nanoparticles in the different bins: rate[k+1] is the rate of growth of nanoparticles in bin k (see eq. (15))
        rateNUC::Float64 = 0;                       # rate of nucleation (see eq. (14))
        prob = zeros(Float64, thresholdBin);        # probabilities of reactions to fire
 
        
        t::Float64 = 0;                             # time

        newpar::Float64 = nucBin/growthBin; 
        
        count::Int64 = 1                            # current step of the algorithm

        while stateMON >= monGrowth                 # check if there are still a-mers available 
            
            # NOTE: all the rates are divided by the growth rate constant growthBin to increase the computational speed

            # update nucleation rate
            rateNUC = newpar*exp(sum(log.(a*stateMON #=a-mers=# .- Int32.(0:(mins-1)))));
 
            # update growth rate
            for i = 1:(largest-1)  
                ratesGR[i+1] = stateMON #=a-mers=# *state[i+1] 
            end
            

            sumratesGR::Int64 = sum(ratesGR);
            sumrates = sumratesGR+rateNUC;

            # sample dt such that t + dt is the time at which the next reaction fires
            dt = rand(Exponential(1/(growthBin*sumrates)));
            t = t + dt;
            times[i,count] = t

            # update probabilities of reactions to fire
            test2 = @view prob[1:largest]
            prob[1] = log(rateNUC[1]);
            for i = 2:largest
                prob[i] = (log(ratesGR[i])); 
            end
            # due to computational approximization, use the methods describes in 
            # https://timvieira.github.io/blog/post/2014/07/31/gumbel-max-trick/
            # to sample from the multinomial distribution and to obtain the index of the reaction that fires
            fires = argmax(test2 .+ rand(Gumbel(0.0,1.0), largest))
            
            
            if fires == thresholdBin
                error("Dimension of nanoparticles is exceding thresholdBin")
            end

            if fires == 1                                       # nucleation occurs
                stateMON = stateMON - monNucleation;            # update number of a-mers
                state[2] = state[2] + parNucleation;            # update number of nanoparticles in the first bin
                traj_mon[i, count] = stateMON
            else                                                # growth of a particle in bin "fires - 1" to particle in bin "fires"
                if fires == largest
                    largest = largest + 1;                      # update the maximum bin size in which we can create particles
                end
                stateMON = stateMON - monGrowth;                # update number of a-mers
                state[fires] = state[fires] - parGrowth;        # update number of nanoparticles in bin "fires - 1"
                state[fires+1] = state[fires+1] + parGrowth;    # update number of nanoparticles in bin "fires"
                traj_mon[i, count] = stateMON
            end
            
            count = count + 1

        end

        out[i, :] .= Int64.(state);
        out[i, 1] = stateMON
    end
    return out, times, traj_mon
end


##########################################################################
##########################################################################
#      TAU-LEAP with POST LEAP CHECKS algorithm - EXACT PROCESS
##########################################################################
##########################################################################

# Function that simulates different trajectories of the process (7) of the paper using the Tau-Leap algorithm.
# The code returns the final size counts for each trajectory. 
function TauLeap_PostLeap(
    nuc::Float64,                       # nucleation rate constant \nu
    growth::Float64,                    # growth rate constant \gamma
    mon::Int64,                         # initial number of monomers N
    mins::Int64,                        # number of monomers used during nucleation n
    repetition::Int64,                  # number of simulated trajectories
    epsilon::Float64,                   # epsilon parameter used in the tau-leap condition
    p::Float64,                         # p parameter to update of tau
    pStar::Float64,                     # pStar parameter to update of tau
    q::Float64;                         # q parameter to update of tau
    threshold::Int64 = -10              # maximum typical size of nanoparticles
)


    # if not given as positive number as input, set the maximum typical size of nanoparticles to 5*N^0.5
    if threshold < 0
        threshold = ceil(Int64, 5 * sqrt(mon))
    end


    # check if the conditions on p, pStar and q are satisfied
    if (p > pStar )
        error("Error, p should be smaller than pStar")
    end

    if (p > 1.0) | (p < 0.0)
        error("Error, p should be in (0, 1)")
    end

    if (pStar > 1.0) | (pStar < 0.0)
        error("Error, pStar should be in (0, 1)")
    end

    if (q > 1.0) | (q < 0.0)
        error("Error, q should be in (0, 1)")
    end


    out::Matrix{Int64} = zeros(Int64, repetition, threshold)    # output matix (repetitions x threshold) with final numer of nanoparticles for each trajectories

    
    Threads.@threads for j = 1:repetition
        println("repetition ", j)
        
        state::Vector{Int64} = zeros(Int64, threshold);         # number of nanoparticles for each size: state[k+1] is the number of nanoparticles of size k
        stateMON::Vector{Int64} = [mon];                        # number of monomers
        largest::Int64 = mins + 1;                              # largest partcle that can be created at the current time

        ratesGR::Vector{Float64} = zeros(Float64, threshold);   # rates of growth of nanoparticles: rate[k+1] is the rate of growth of nanoparticles of size k (see eq. (5) and (6) of the paper)
        rateNUC::Vector{Float64} = [0.0];                       # rate of nucleation (see eq. (4) of the paper)
        
        
        # vectors that store the current internal times and values of the Poisson processes and the interla time and values proposed for the future time
        T::Vector{Float64} = zeros(Float64, threshold)          
        C::Vector{Int64} = zeros(Int64, threshold)              
        ST::Matrix{Float64} = zeros(Float64, 200, threshold)    
        SC::Matrix{Int64} = zeros(Int64, 200, threshold)        

        nPois::Vector{Int64} = zeros(Int64, threshold)          # independent random variable (Poisson/Binomial distriuted)

        # auxiliary variables used in algorithm (see Algorithms (2) and (3))
        B::Vector{Int64} = zeros(Int64, threshold)    
        rows::Vector{Int64} = zeros(Int64, threshold)  
        dimS::Vector{Int64} = ones(Int64, threshold)
        update::Vector{Int64} = zeros(Int64, threshold)         


        sequenza = Int64.(0:(mins-1))
        
        # update the rates
        rates!(stateMON, state, rateNUC, ratesGR, mins, nuc, growth, sequenza, largest)

        # absolute time of the system
        t::Float64 = 0.0;
        
        # compute tau and g_i
        output::NamedTuple = computeTAU(stateMON, state, rateNUC, ratesGR, mon, mins, epsilon, threshold, largest)  
        tau::Float64 = output.tau       
        g::Vector{Float64} = output.g
        newT::Float64 = 0.0                         # new proposed internal time

        r::Float64 = 0.0
        i::Int64 = 0

        epsilonStar::Float64 = 3/4*epsilon

        while stateMON[1] > 0 # check if there are still available monomers
            
            # nucleation
            B[1] = dimS[1]
            
            newT = rateNUC[1]*tau + T[1] 
            if newT >= ST[B[1], 1]                                                          # proposed internal time is largest than all the previously proposed internal times
                nPois[1] = rand(Poisson(newT - ST[B[1], 1])) + SC[B[1], 1]- C[1]            # sample from Poisson 
                rows[1] = B[1]                      
            else                                                                            # proposed internal time is between two previously proposed times
                i = searchsortedlast(ST[1:B[1], 1], newT)  
                r = (newT - ST[i, 1])/(ST[i + 1, 1] - ST[i, 1])                             # success probability of the Binomial distribution
                nPois[1] = rand(Binomial(SC[i + 1, 1] - SC[i, 1], r)) + SC[i, 1] - C[1]     # sample from Binomial 
                rows[1] = i 
            end


            # growth from particle of dimension k-1 to dimension k, with k = 2,..., largest
            for k = 2:largest
                
                B[k] = dimS[k]
                newT = ratesGR[k]*tau + T[k]

                if newT >= ST[B[k], k]                                                      # proposed internal time is largest than all the previously proposed internal times
                    nPois[k] = rand(Poisson(newT - ST[B[k], k])) + SC[B[k], k]- C[k]        # sample from Poisson
                    rows[k] = B[k]
                else                                                                        # proposed internal time is between two previously proposed times
                    i = searchsortedlast(ST[1:B[k], k], newT)                    
                    r = (newT - ST[i, k])/(ST[i + 1, k] - ST[i, k])                         # success probability of the Binomial 
                    nPois[k] = rand(Binomial(SC[i + 1, k] - SC[i, k], r)) + SC[i, k] - C[k] # sample from Binomial 

                    rows[k] = i 
                end
            end          


            accept::Bool = true

            # compute the proposed update of system
            updateState!(update, nPois, mins, largest)

            # check if monomers do not satisty leap-condition
            if (abs(update[1]) > max(1.0, epsilon*stateMON[1]/g[1]))
                accept = false
            end
            
            # check if nanoparticles of any dimension satisfy the leap-condition
            i = mins 
            while ((i <= largest) & (accept == true))
                if (abs(update[i + 1]) > max(1.0, epsilon*state[i + 1]/g[i + 1]))
                    i = largest + 1 # nanoparticles of dimension i do not satisty leap-condition and we go out of the cycle
                    accept = false
                else
                    i = i + 1
                end
            end
            

            if accept                       # leap condition satisfied 
                ST[2:(dimS[1] - rows[1] + 1), 1] .= ST[(rows[1]+ 1):dimS[1], 1]
                SC[2:(dimS[1] - rows[1] + 1), 1] .= SC[(rows[1]+ 1):dimS[1], 1]
                dimS[1] = dimS[1] - rows[1] + 1
                ST[1, 1] = T[1] + rateNUC[1]*tau
                SC[1, 1] = C[1] + nPois[1]

                T[1] = ST[1, 1]             # update current internal time of the Poisson process for nucleation
                C[1] = SC[1, 1]             # update current value of the Poisson process for nucleation
                
                for k = 2:largest                  
                    ST[2:(dimS[k] - rows[k] + 1), k] .= ST[(rows[k]+ 1):dimS[k], k]
                    SC[2:(dimS[k] - rows[k] + 1), k] .= SC[(rows[k]+ 1):dimS[k], k]
                    dimS[k] = dimS[k] - rows[k] + 1 
                    ST[1, k] = T[k] + ratesGR[k]*tau
                    SC[1, k] = C[k] + nPois[k]

                    T[k] = ST[1, k]         # update current internal time of the Poisson process for growth
                    C[k] = SC[1, k]         # update current value of the Poisson process for growth
                    
                end

                t = t + tau                 # update absolute time of the system

                # check if the leap condition is satisfied also with epsilonStar = 3/4*epsilon
                accept = true
                
                if (abs(update[1]) > max(1.0, epsilonStar*stateMON[1]/g[1]))
                    accept = false
                end
                
                i = mins 
                while ((i <= largest) & (accept == true)) 
                    if (abs(update[i + 1]) > max(1.0, epsilonStar*state[i + 1]/g[i + 1]))
                        i = largest + 1
                        accept = false
                    else
                        i = i + 1
                    end
                end

                # update tau
                if accept
                    tau = tau^q
                else 
                    tau = tau*pStar
                end

                # update the state of the system
                stateMON[1] = stateMON[1] + update[1]
                state[2:(largest+1)] = state[2:(largest+1)] .+ update[2:(largest+1)]
                if state[largest + 1] > 0
                    largest = largest + 1
                end
           
                # compute the new rates
                rates!(stateMON, state, rateNUC, ratesGR, mins, nuc, growth, sequenza, largest)
               
                # update g_1
                if (stateMON[1] < 100)
                    g[1] = 3 + 1/(stateMON[1] - 1) + 1/(stateMON[1] - 2)
                end 
            else                            # leap condition not satisfied 

                # store the new proposed internal times and values of the Poisson process
                if  ((dimS[1] + 1) > size(dimS)[1])
                    ST = [ST;transpose(zeros(Float64, threshold))]
                    SC = [SC;transpose(zeros(Int64, threshold))]
                end

                ST[(rows[1] + 2):(dimS[1] + 1),  1] = ST[(rows[1] + 1):(dimS[1]), 1]
                SC[(rows[1] + 2):(dimS[1] + 1),  1] = SC[(rows[1] + 1):(dimS[1]), 1]
                ST[rows[1] + 1, 1] = T[1] + rateNUC[1]*tau
                SC[rows[1] + 1, 1] = C[1] + nPois[1]
                dimS[1] = dimS[1] + 1   

                for k = 2:largest
                    ST[(rows[k] + 2):(dimS[k] + 1),  k] = ST[(rows[k] + 1):(dimS[k]), k]
                    SC[(rows[k] + 2):(dimS[k] + 1),  k] = SC[(rows[k] + 1):(dimS[k]), k]
                    ST[rows[k] + 1, k] = T[k] + ratesGR[k]*tau
                    SC[rows[k] + 1, k] = C[k] + nPois[k]
                    dimS[k] = dimS[k] + 1   
                end

                tau = tau*p                 # update tau
            end
        end

        out[j, :] = state
        out[j, 1] = stateMON[1]
    end
    
    return out    
end

# It computes the rates for nucleation and growth for process (7)
function rates!(
    stateMON::Vector{Int64},                # current number of monomers
    state::Vector{Int64},                   # current number of nanoparticles of different sizes
    rateNUC::Vector{Float64},               # nucleation rate
    ratesGR::Vector{Float64},               # growth rates
    mins::Int64,                            # number of monomers used during nucleation n
    nuc::Float64,                           # nucleation rate constant \nu
    growth::Float64,                        # growth rate constant \gamma
    sequenza::Vector{Int64},                # auxialiar paramer with 0:(mins-1)
    largest::Int64,                         # largest partcle that can be created at the current time
)

    # nucleation rate
    if stateMON[1] >= mins
        rateNUC[1] = nuc*exp(sum(log.(stateMON[1] .- sequenza)));
    else 
        rateNUC[1] =0.0;
    end 
    
    # growth rate
    for i = mins:(largest-1)  
        ratesGR[i+1] = growth*stateMON[1] * state[i+1]; 
    end
  
end


# It returns the functions g_i(X_i(t))
function computeG(
    stateMON::Vector{Int64},
    mins::Int64,
    threshold::Int64
)::Vector{Float64}

    # compute HOR
    hor::Vector{Int64} = zeros(Int64, threshold)
    
    # all the growth reactions are bimolecular reactions and requires only one molecule of Xi
    hor[(mins + 1):threshold] .= 2 
    hor[1] = mins

    # compute the g_i
    g::Vector{Float64} = zeros(Float64, threshold)
    g[(mins + 1):threshold] .= 2

    if (hor[1] == 1)
        g[1] = 1
    elseif (hor[1] == 2)
        g[1] = 2 + 1/(stateMON[1] - 1)
    elseif (hor[1] == 3)
        g[1] = 3 + 1/(stateMON[1] - 1) + 2/(stateMON[1] - 2)
    end

    g[2:mins] .= 10000
    
    return g
end


# It returns the initial tau and the g_i
function computeTAU(
    stateMON::Vector{Int64},
    state::Vector{Int64},
    rateNUC::Vector{Float64},
    ratesGR::Vector{Float64},
    mon::Int64,
    mins::Int64,
    epsilon::Float64,  
    threshold::Int64, 
    largest::Int64,
)

    # compute g_i
    g::Vector{Float64} = computeG(stateMON, mins, threshold)

    # compute mu and sigma2 to obtain the intial tau
    num1::Vector{Float64} = max.(epsilon.*state ./ g, 1.0)
    num1[1] = max(epsilon*stateMON[1]./g[1], 1.0)
    num2::Vector{Float64} = num1.^2
    mu::Vector{Float64} = zeros(Float64, threshold)
    sigma2::Vector{Float64} = zeros(Float64, threshold)

    mu[1] = -mins*rateNUC[1] - sum(ratesGR[(mins + 1): largest])
    sigma2[1] = (mins)^2*rateNUC[1] + sum(ratesGR[(mins + 1): threshold])
    mu[mins + 1] = rateNUC[1] - ratesGR[mins + 1]
    sigma2[mins + 1] = rateNUC[1] + ratesGR[mins + 1]
    for i = (mins + 1):(largest-1)
        mu[i + 1] = ratesGR[i] - ratesGR[i + 1]
        sigma2[i + 1] = ratesGR[i] + ratesGR[i + 1]
    end
    
    min1::Vector{Float64} = num1./abs.(mu)
    min2::Vector{Float64} = num2./sigma2
    
    # compute tau
    tau::Float64 = min(minimum(min1), minimum(min2))

    return (tau = tau, g = g)

end

# It computes the updated value of the system 
function updateState!(
    update::Vector{Int64},
    nPois::Vector{Int64},
    mins::Int64,
    largest::Int64,
)

    update[1] = - nPois[1]*mins - sum(nPois[(mins + 1):largest])
    update[mins + 1] = nPois[1] - nPois[mins + 1]

    for i = (mins + 1):(largest - 1)
        update[i + 1] = nPois[i] - nPois[i + 1]
    end

    update[largest + 1] = nPois[largest] # first formation of a particle of dimension largest 
end












#######################################################################################
#######################################################################################
#      TAU-LEAP with POST LEAP CHECKS algorithm algorithm - APPROXIMATED PROCESS
#######################################################################################
#######################################################################################

# Function that simulates different trajectories of the approximated process (13) of the paper using the Gillespie algorithm.
# The code returns the final size counts for each trajectory. 
function TauLeap_PostLeapAPPROX(
    nucBin::Float64,                    # nucleation rate constant 
    growthBin::Float64,                 # growth rate constant
    mon::Int64,                         # inital number of monomers N 
    b::Int64,                           # total numer of bins
    mins::Int64,                        # number of monomers used during nucleation n in the associated exact process (7)
    repetition::Int64,                  # number of simulated trajectories
    epsilon::Float64,                   # epsilon parameter used in the tau-leap condition
    p::Float64,                         # p parameter to update of tau
    pStar::Float64,                     # pStar parameter to update of tau
    q::Float64;                         # q parameter to update of tau
    thresholdBin::Int64 = -10,          # typical maximum number of filled bins
    monNucleation::Int32 = 1,           # number of a-mers consumed during nucleation
    parNucleation::Int16 = 1,           # number of nanoparticles removed from each bin during nucleation
    monGrowth::Int32 = 1,               # number of a-mers consumed during growth
    parGrowth::Int16                    # number of nanoparticles removed from each bin during growth
)

    
    # ----- le reazioni che consideriamo sono 
    # ----- nuc: M_d --> C1 con rateConstant = \nu --- per fare come il gillespie, il rateNUC sarà nuc*(d*stateMON)^(mins-1)
    # ----- gro: M + C_i --> C_{i+1} con rateConstant = \growth (viene data già in input come growthDelNonApprox/d) --- per fare come il gillespie il ratesGR sarà growth*d*stateMON*C_i
    # ----- Dobbiamo ricordarci che in entrambi i casi vengono consumati d monomeri

    a::Int64 = ceil(Int64, mon^(1/4)) # # bin amplitude

    # if not given as positive number as input, set the typical maximum number of filled bins to 5*N^0.5/a
    if thresholdBin < 0
        thresholdBin = ceil(Int64, 5 * sqrt(mon)/a)
    end

    # check if the conditions on p, pStar and q are satisfied
    if (p > pStar )
        error("Error, p should be smaller than pStar")
    end

    if (p > 1.0) | (p < 0.0)
        error("Error, p should be in (0, 1)")
    end

    if (pStar > 1.0) | (pStar < 0.0)
        error("Error, pStar should be in (0, 1)")
    end

    if (q > 1.0) | (q < 0.0)
        error("Error, q should be in (0, 1)")
    end


    out::Matrix{Int64} = zeros(Int64, repetition, thresholdBin) # output matix (repetitions x threshold) with final numer of nanoparticles for each trajectories


    Threads.@threads for j = 1:repetition
        println("repetition ", j)
        
        state::Vector{Int16} = zeros(Int16, thresholdBin);      # number of nanoparticles for each bin: state[k+1] is the number of nanoparticles of in the bin k
        stateMON::Vector{Int64} = [ceil(Int64, mon/a)];         # number of a-mers
        largest::Int64=2;                                       # largest bin that can be filled at the current time
        
        ratesGR::Vector{Float64} = zeros(Int64, thresholdBin);  # rates of growth of nanoparticles in the different bins: rate[k+1] is the rate of growth of nanoparticles in bin k (see eq. (15))
        rateNUC::Vector{Float64} = [0.0];                       # rate of nucleation (see eq. (14))

       
        
        # vectors that store the current internal times and values of the Poisson processes and the interla time and values proposed for the future time
        T::Vector{Float64} = zeros(Float64, thresholdBin) 
        C::Vector{Int64} = zeros(Int64, thresholdBin)
        ST::Matrix{Float64} = zeros(Float64, 100, thresholdBin)
        SC::Matrix{Int64} = zeros(Int64, 100, thresholdBin)

        nPois::Vector{Int64} = zeros(Int64, thresholdBin)          # independent random variable (Poisson/Binomial distriuted)

        # auxiliary variables used in algorithm (see Algorithms (2) and (3))
        B::Vector{Int64} = zeros(Int64, thresholdBin)
        rows::Vector{Int64} = zeros(Int64, thresholdBin)
        dimS::Vector{Int64} = ones(Int64, thresholdBin)
        update::Vector{Int64} = zeros(Int64, thresholdBin)

        
        sequenza = Int32.(0:(mins-1))  
        
        
        

        # update the rates
        rates_APPROX!(stateMON, state, rateNUC, ratesGR, nucBin, growthBin, sequenza, largest, monGrowth, a)

        # absolute time of the system
        t::Float64 = 0.0;
        
        # compute tau and g_i
        output::NamedTuple = computeTAU_APPROX(stateMON, state, rateNUC, ratesGR, mon, mins, epsilon, thresholdBin, largest)  
        tau::Float64 = output.tau 
        g::Vector{Float64} = output.g

        newT::Float64 = 0.0                                     # new proposed internal time
        r::Float64 = 0.0
        i::Int64 = 0

        epsilonStar::Float64 = 3/4*epsilon

        while (stateMON[1] >= monGrowth)                        # check if there are still available a-mers
           
            # nucleation
            B[1] = dimS[1]
            
            newT = rateNUC[1]*tau + T[1]
            if newT >= ST[B[1], 1]                                                          # proposed internal time is largest than all the previously proposed internal times
                nPois[1] = rand(Poisson(newT - ST[B[1], 1])) + SC[B[1], 1]- C[1]            # sample from Poisson 
                rows[1] = B[1]                      
            else                                                                            # proposed internal time is between two previously proposed times
                
                i = searchsortedlast(ST[1:B[1], 1], newT)                  
                r = (newT - ST[i, 1])/(ST[i + 1, 1] - ST[i, 1])                             # success probability of the Binomial distribution
                nPois[1] = rand(Binomial(SC[i + 1, 1] - SC[i, 1], r)) + SC[i, 1] - C[1]     # sample from Binomial 

                rows[1] = i 
            end


            # growth from particle of dimension k-1 to dimension k, with k = 2,..., largest
            for k = 2:largest
                
                B[k] = dimS[k]
                newT = ratesGR[k]*tau + T[k]

                if newT >= ST[B[k], k]                                                      # proposed internal time is largest than all the previously proposed internal times
                    nPois[k] = rand(Poisson(newT - ST[B[k], k])) + SC[B[k], k]- C[k]        # sample from Poisson
                    rows[k] = B[k]
                else                                                                        # proposed internal time is between two previously proposed times
                    
                    i = searchsortedlast(ST[1:B[k], k], newT) 
                    r = (newT - ST[i, k])/(ST[i + 1, k] - ST[i, k])                         # success probability of the Binomial 
                    nPois[k] = rand(Binomial(SC[i + 1, k] - SC[i, k], r)) + SC[i, k] - C[k] # sample from Binomial 

                    rows[k] = i 
                end
            end          


            accept::Bool = true

            # compute the proposed update of system
            updateState_APPROX!(update, nPois, mins, largest)
            
            # check if monomers do not satisty leap-condition
            if (abs(update[1]) > max(1.0, epsilon*stateMON[1]/g[1]))
                accept = false
            end
            
            # check if nanoparticles in any bin satisfy the leap-condition
            i = 1
            while ((i <= largest) & (accept == true))
                if (abs(update[i + 1]) > max(1.0, epsilon*state[i + 1]/g[i + 1]))
                    i = largest + 1                 # nanoparticles in bin i do not satisty leap-condition and we go out of the cycle
                    accept = false
                else
                    i = i + 1
                end
            end
            
            
            if accept                       # leap condition satisfied 
                ST[2:(dimS[1] - rows[1] + 1), 1] .= ST[(rows[1]+ 1):dimS[1], 1]
                SC[2:(dimS[1] - rows[1] + 1), 1] .= SC[(rows[1]+ 1):dimS[1], 1]
                dimS[1] = dimS[1] - rows[1] + 1
                ST[1, 1] = T[1] + rateNUC[1]*tau
                SC[1, 1] = C[1] + nPois[1]

                T[1] = ST[1, 1]             # update current internal time of the Poisson process for nucleation
                C[1] = SC[1, 1]             # update current value of the Poisson process for nucleation
                
                for k = 2:largest                  
                    ST[2:(dimS[k] - rows[k] + 1), k] .= ST[(rows[k]+ 1):dimS[k], k]
                    SC[2:(dimS[k] - rows[k] + 1), k] .= SC[(rows[k]+ 1):dimS[k], k]
                    dimS[k] = dimS[k] - rows[k] + 1
                    ST[1, k] = T[k] + ratesGR[k]*tau
                    SC[1, k] = C[k] + nPois[k]

                    T[k] = ST[1, k]         # update current internal time of the Poisson process for growth
                    C[k] = SC[1, k]         # update current value of the Poisson process for growth                    
                end

                t = t + tau                 # update absolute time of the system

                accept = true
                i = 1

                # monomers do not satisty leap-condition
                if (abs(update[1]) > max(1.0 #=monomeroni=#, epsilonStar*stateMON[1]/g[1]))
                    accept = false
                end
                
                # check if the leap condition is satisfied also with epsilonStar = 3/4*epsilon
                i = 1
                while ((i <= largest) & (accept == true))
                    if (abs(update[i + 1]) > max(1.0, epsilonStar*state[i + 1]/g[i + 1]))
                        i = largest + 1
                        accept = false
                    else
                        i = i + 1
                    end
                end

                # update tau
                if accept
                    tau = tau^q
                else 
                    tau = tau*pStar
                end

                # update the state of the system
                stateMON[1] = stateMON[1] + update[1]
                state[2:(largest+1)] = state[2:(largest+1)] .+ update[2:(largest+1)]
                if state[largest + 1] > 0
                    largest = largest + 1
                end
           
                # compute the new rates
                rates_APPROX!(stateMON, state, rateNUC, ratesGR, nucBin, growthBin, sequenza, largest, monGrowth, a)

            else                            # leap condition not satisfied 
                
                # store the new proposed internal times and values of the Poisson process
                ST[(rows[1] + 2):(dimS[1] + 1),  1] = ST[(rows[1] + 1):(dimS[1]), 1]
                SC[(rows[1] + 2):(dimS[1] + 1),  1] = SC[(rows[1] + 1):(dimS[1]), 1]
                ST[rows[1] + 1, 1] = T[1] + rateNUC[1]*tau
                SC[rows[1] + 1, 1] = C[1] + nPois[1]
                dimS[1] = dimS[1] + 1   

                for k = 2:largest
                    ST[(rows[k] + 2):(dimS[k] + 1),  k] = ST[(rows[k] + 1):(dimS[k]), k]
                    SC[(rows[k] + 2):(dimS[k] + 1),  k] = SC[(rows[k] + 1):(dimS[k]), k]

                    ST[rows[k] + 1, k] = T[k] + ratesGR[k]*tau
                    SC[rows[k] + 1, k] = C[k] + nPois[k]

                    dimS[k] = dimS[k] + 1   
                end

                tau = tau*p                 # update tau
            end
        end

        out[j, :] = state
        out[j, 1] = stateMON[1]
    end
    
    return out    
end



# It computes the rates for nucleation and growth for process (13)
function rates_APPROX!(
    stateMON::Vector{Int64},            # current number of a-mers
    state::Vector{Int16},               # current number of nanoparticles in each bin
    rateNUC::Vector{Float64},           # nucleation rate
    ratesGR::Vector{Float64},           # growth rates
    nucBin::Float64,                    # nucleation rate constant 
    growthBin::Float64,               # growth rate constant 
    sequenza::Vector{Int32},            # auxialiar paramer with 0:(mins-1)
    largest::Int64,                     # largest bin that can be filled at the current time
    monGrowth::Int32,                   # number of a-mers consumed during growth
    a::Int64,                           # bin amplitude
)

    # nucleation rate
    if stateMON[1] >= monGrowth
        rateNUC[1] = nucBin*exp(sum(log.(a*stateMON[1] #=a-mers=# .- sequenza)));
    else 
        rateNUC[1] = 0.0;
    end 
    
    # growth rate
    for i = 1:(largest-1)  
        ratesGR[i+1] = growthBin* stateMON[1] #=monomeroni=# * state[i+1] 
    end
  
end

# It returns the functions g_i(X_i(t)) for the approximated process (13)
function computeG_APPROX(
    stateMON::Vector{Int64},
    thresholdBin::Int64
)::Vector{Float64}

    # compute HOR
    hor::Vector{Int64} = zeros(Int64, thresholdBin)
    
    # all the growth reactions are bimolecular reactions and requires only one molecule of Xi
    hor[2:thresholdBin] .= 2 
    hor[1] = 1

    # compute the g_i
    g::Vector{Float64} = zeros(Float64, thresholdBin)
    g[2:thresholdBin] .= 2

    if (hor[1] == 1)
        g[1] = 1
    elseif (hor[1] == 2)
        g[1] = 2 + 1/(stateMON[1] - 1)
    elseif (hor[1] == 3)
        g[1] = 3 + 1/(stateMON[1] - 1) + 2/(stateMON[1] - 2)
    end

    return g
end


# It returns the initial tau and the g_i for the approximated process (13)
function computeTAU_APPROX(
    stateMON::Vector{Int64},
    state::Vector{Int16},
    rateNUC::Vector{Float64},
    ratesGR::Vector{Float64},
    mon::Int64,
    mins::Int64,
    epsilon::Float64,  
    thresholdBin::Int64, 
    largest::Int64,
)

    # compute g_i
    g::Vector{Float64} = computeG_APPROX(stateMON, thresholdBin)

    # compute mu and sigma2 to obtain the initial tau
    num1::Vector{Float64} = max.(epsilon.*state ./ g, 1.0)
    num1[1] = max(epsilon*stateMON[1]./g[1], 1.0)
    num2::Vector{Float64} = num1.^2
    mu::Vector{Float64} = zeros(Float64, thresholdBin)
    sigma2::Vector{Float64} = zeros(Float64, thresholdBin)

    mu[1] = -rateNUC[1] - sum(ratesGR[2: largest])
    sigma2[1] = rateNUC[1] + sum(ratesGR[2: thresholdBin])
    mu[2] = rateNUC[1] - ratesGR[2]
    sigma2[2] = rateNUC[1] + ratesGR[2]
    for i = 2:(largest-1)
        mu[i + 1] = ratesGR[i] - ratesGR[i + 1]
        sigma2[i + 1] = ratesGR[i] + ratesGR[i + 1]
    end
    
    min1::Vector{Float64} = num1./abs.(mu)
    min2::Vector{Float64} = num2./sigma2
    
    # compute tau
    tau::Float64 = min(minimum(min1), minimum(min2))

    return (tau = tau, g = g)
end

# It computes the updated value of the system for process in (13)
function updateState_APPROX!(
    update::Vector{Int64},
    nPois::Vector{Int64},
    mins::Int64,
    largest::Int64
)

    update[1] = - nPois[1] - sum(nPois[2:largest]) 
    update[2] = nPois[1] - nPois[2]

    for i = 2:(largest - 1)
        update[i + 1] = nPois[i] - nPois[i + 1]
    end

    update[largest + 1] = nPois[largest] # first formation of a nanparticle in the largest-th bin (contained in position largest + 1)
end







