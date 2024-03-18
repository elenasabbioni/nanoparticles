################################
# Code to produce Figure 1 of 
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

for(mon in c(1e+4, 1e+6, 1e+8)){
  # Set initial parameters
  mon <- 1e+4         # initial number of monomers for process (7)
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

  # Loop across the different repetitions to find the minimum and the maximum states  with particles'concentration greater than 0 for original Gillespie  simulations
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


  #----- FIGURE 1 OF THE PAPER -------
  # plots for the Gillespie's algorithm
  gillespiePLOT <- list()
  for(r in 1:rep){
      pdf(file=paste(pathOutput, "Gillespie10_", chrMon,"-",r , ".pdf", sep = ""), width = 9, height = 5)

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
        xText <- 400 
      }else if(chrMon == "4"){
        xText <- 40
      }else if(chrMon =="8"){
        xText <- 4000
      }

      lab <- sprintf("$N = 10^{%d}, rep = %d$",chrMonInt, r)

      gillespiePLOT[[r]] <- ggplot(subset(matGill, rep == r & nbin %in% limX), aes(x = factor(nbin), y = counts)) + geom_bar(stat="identity",  color="black", fill="black") + ylim(c(0,6)) +  labs(x = labX,  y = labY) + scale_x_discrete(breaks = limX[index], labels =  sprintf("%.1f",labels))  + textThemeOriginal + theme(plot.margin = margin(0, 0.5, 0, 0.5, unit = "cm")) + annotate(geom="text", x=xText, y=5.5, label=TeX(lab), color="black", size = 8)


      print(gillespiePLOT[[r]])
      dev.off()
  }
  # ------------------------------------------
  # ------------------------------------------






}
