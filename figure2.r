################################
# Code to produce Figure 2 of 
# Sabbioni E., Szabò R., Siri P., Cappelletti D., Lente G., Bibbona E. 
# "Final nanoparticle size distribution under unusual parameter regimes", 
# arXiv...
################################

# Clear the workspace
rm(list=ls())

########################################
#           PACKAGES
########################################
# Load necessary libraries
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(dplyr)
library(scales)

########################################
#           GRAPHICAL PARAMETERS
########################################
textThemeOriginal <- theme(
    plot.title = element_text(
    size = 30,          
    color = "black",     
    face = "bold",     
    hjust = 0.5
  ), 
    axis.title.x = element_text(
    size = 30,          
    color = "black",     
    face = "bold",     
  ), 
  axis.title.y = element_text(
    size = 30,          
    color = "black",     
  ), 
  axis.text.x = element_text(
    size = 23,      
    color = "black"  
  ), 
  axis.text.y = element_text(
    size = 23,      
    color = "black" 
  )
  )



########################################
#           INITIALIZATION 
#           OF PARAMETERS
########################################
# Set initial parameters
mon <- 1e+8         # initial number of monomers for process (7)
growth <- 5         # growth rate constant 
mins <- 3           # number of monomers used during nucleation n in the exact process (7)


nu <- growth * (mon^(1-mins))              # nucleation rate constant 

threshold <- ceiling(growth * sqrt(mon))   # typical maximum size of nanoparticles



chrMon <- format(mon, scientific = TRUE)
chrMon <- substr(chrMon, nchar(chrMon), nchar(chrMon))
chrMonInt <- as.integer(chrMon)


########################################
#              LOAD THE RESULTS
#              OF THE SIMULATIONS
########################################
# Set the path where the simulation results are
pathInput <- "nanoparticlesProject/code/output"     

# Read the Gillespie simulation data from a CSV file
setGill <- fread(file=paste(pathInput, "/Gillespie/gillespie_", chrMon, ".csv", sep = ""))
# Convert data table to a matrix
matGill <- data.frame(rep = rep(seq(1, dim(setGill)[1]), times = dim(setGill[, -1])[2]), counts = as.vector(unname(as.matrix(setGill[,-1]))), nbin = rep(seq(1, dim(setGill[,-1])[2]), each = dim(setGill)[1]))

# Perform a sanity check to ensure the total mass at the end is consistent with the initial mass
N <- mon
for (i in unique(matGill$rep)){
  if (sum(setGill[i,-1] * c(1, 1:(N-1))) != N) {
    stop("Error in the mass conservation for Gillespie simulation!")  
  } 
} 



########################################
#         BUILD BINNED RESULTS
########################################
# Set initial values for repetitions and bins
rep <- dim(setGill)[1]
a <- ceiling(N^(1/4))
bins <- length(setGill[1,-1])%/% a

# Initialize a numeric vector for binned state in Gillespie simulations
binned_stateGill <- data.frame(rep = rep(seq(1, rep), times = bins), counts = 0, nbin = rep(seq(1, bins), each = rep))

# Bin together different states in order to visualize better the results of the Gillespie algorithm 
for (r in 1:rep){
    for (b in 1:(bins-1)){
        # Calculate the sum of particles within the b-th bin
        binned_stateGill$counts[which(binned_stateGill$rep == r & binned_stateGill$nbin == b)] <- sum(matGill$counts[matGill$rep == r][1+(1:a) + (b-1)*a])
    }
        
    # Calculate the sum of monomers for the last bin 
    binned_stateGill$counts[which(binned_stateGill$rep == r & binned_stateGill$nbin == bins)] <- sum(matGill$counts[matGill$rep == r][1+(1:a-1) + (bins-1)*a])
}






########################################
#                PLOT 
########################################
# set the path where the plots will be saved
pathOutput <- "nanoparticlesProject/graphics"     # path where the plots will be stored



# ----------------------------------------------------------------
# ---------- original Gillespie's results ---------
# ----------------------------------------------------------------

# Initialize variables to find relevant states in the original simulation
minnGill <- N
maxxGill <- 0


# Loop across the different repetitions to find the minimum and the maximum states  with particles'concentration greater than 0 for original Gillespie simulations
for (r in 1:rep){ 
    posintGill <- which(matGill$counts[matGill$rep == r] > 0)
    minnGill <- min(c(posintGill, minnGill))
    maxxGill <- max(c(posintGill, maxxGill))
}

# Create two vectors of relevant states in the original Gillespie simulation respectively
relevant_statesGill <- minnGill:maxxGill

# Scale the states we have found by sqrt(N), such that we obtain a comparable scale when N grows
scaled_relevant_statesGill <- round(relevant_statesGill / sqrt(N), 1)
saveGill <- which(scaled_relevant_statesGill <= 2.05)
relevant_statesGill <- relevant_statesGill[saveGill]
scaled_relevant_statesGill <- scaled_relevant_statesGill[saveGill]



# set the graphical elements
limX <- seq(0, max(relevant_statesGill))
scaled_limX <- round(limX/sqrt(N), 1)
labels <- round(seq(0, max(scaled_limX), 0.2), 1)
index <- rep(NA, length(labels))
for(i in 1:length(labels)){
  index[i] <- which(round(scaled_limX, 1) == labels[i])[1]
}


# ----------------------------------------------------------------
# ---------- binned Gillespie's results -----------
# ----------------------------------------------------------------

# Initialize variables to find relevant bins in the binned simulation
minnGill_binned <- bins
maxxGill_binned <- 0


# Loop across the different repetitions to find the minimum and the maximum states  with particles'concentration greater than 0 for binned Gillespie simulations
for (r in 1:rep){ 
    posintGill_binned <- which(binned_stateGill$counts[binned_stateGill$rep == r]>0)
    minnGill_binned <- min(c(posintGill_binned, minnGill_binned))
    maxxGill_binned <- max(c(posintGill_binned, maxxGill_binned))
}

# Create two vectors of relevant states in the original Gillespie simulation respectively
relevant_statesGill_binned <- minnGill_binned:maxxGill_binned

# Scale the states we have found by sqrt(N)/d, such that we obtain a comparable scale when N grows
scaled_relevant_statesGill_binned <- round(relevant_statesGill_binned*a / sqrt(N), 1)
saveGill <- which(scaled_relevant_statesGill_binned <= 2.05)
relevant_statesGill_binned <- relevant_statesGill_binned[saveGill]
scaled_relevant_statesGill_binned <- scaled_relevant_statesGill_binned[saveGill]

# graphical parameters
limX_binned <- seq(min(relevant_statesGill_binned), max(relevant_statesGill_binned))
scaled_limX_binned <- round(limX_binned*a/sqrt(N), 1)
labels_binned <- round(seq(0, max(scaled_limX_binned), 0.2), 1)
index_binned <- rep(NA, length(labels_binned))
for(i in 1:length(labels_binned)){
  index_binned[i] <- which(round(scaled_limX_binned, 1) == labels_binned[i])[1]
}

#----- FIGURE 2 OF THE PAPER -------
# plots of the results of the Gillespie's algorithm after the binning
gillespiePLOT_binned <- list()
for(r in 1:rep){
    pdf(file=paste(pathOutput, "GillespieBinned10_", chrMon,"-",r , ".pdf", sep = ""), width = 9, height = 5)

    if(r == 1){
      labY <- "Counts"
    }else{
      labY <- ""
    }

    if(chrMon == "8"){
      labX <- expression(paste("Scaled nanoparticle size (in the unit of ", sqrt(N), ")"), sep = "")
    }else{
      labX <- ""
    }

    if(chrMon =="6"){
      xText <- 40 
      yText <- max(binned_stateGill$counts) - 0.5
    }else if(chrMon == "4"){
      xText <- 7
      yText <- max(binned_stateGill$counts) - 0.5
    }else if(chrMon =="8"){
      xText <- 70
      yText <- 115
    }

    lab <- sprintf("Binned process, $N = 10^{%d}, rep = %d$",chrMonInt, r)
    labels <- sprintf("%.1f",labels_binned)
    labels[1] <- ""

    gillespiePLOT_binned[[r]] <- ggplot(subset(binned_stateGill, rep == r & nbin %in% limX_binned), aes(x = factor(nbin), y = counts)) + geom_bar(stat="identity",  color="black", fill="black") + ylim(c(0,max(binned_stateGill$counts)))+ labs(x =  labX,  y = labY) + scale_x_discrete(breaks = limX_binned[index_binned], labels =  labels)  + textThemeOriginal + theme(plot.margin = margin(0, 0.5, 0, 0.5, unit = "cm")) +  annotate(geom="text", x=xText, y=yText, label=TeX(lab), color="black", size = 8)
    print(gillespiePLOT_binned[[r]])
    dev.off()
}
# ------------------------------------------
# ------------------------------------------


