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
textThemeApprox <- theme(
    plot.title = element_text(
    size = 30,          
    color = "black",     
    face = "bold",     
    hjust = 0.5
  ), 
    axis.title.x = element_text(
    size = 20,          
    color = "black",     
    face = "bold",     
  ), 
  axis.title.y = element_text(
    size = 20,          
    color = "black",     
  ), 
  axis.text.x = element_text(
    size = 20,       
    color = "black"  
  ), 
  axis.text.y = element_text(
    size = 20,      
    color = "black" 
  )
  )

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
mon <- 1e+4         # initial number of monomers for process (7)
monApp <- 1e+04     # initial number of monomers for the approximated process (13)
growth <- 5         # growth rate constant 
mins <- 3           # number of monomers used during nucleation n in the exact process (7)


nu <- growth * (mon^(1-mins))              # nucleation rate constant 

threshold <- ceiling(growth * sqrt(mon))   # typical maximum size of nanoparticles



chrMon <- format(mon, scientific = TRUE)
chrMon <- substr(chrMon, nchar(chrMon), nchar(chrMon))
chrMonInt <- as.integer(chrMon)
chrMonApp <- format(monApp, scientific = TRUE)
chrMonApp <- substr(chrMonApp, nchar(chrMonApp)-1, nchar(chrMonApp))
if(substr(chrMonApp, nchar(chrMonApp)-1, nchar(chrMonApp)-1) == "0"){
  chrMonApp <- substr(chrMonApp, nchar(chrMonApp), nchar(chrMonApp))
}

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

# Read the Tau-Leap simulation data from a CSV file
setTau <- fread(file=paste(pathInput, "/TauLeap/tau-leap_", chrMon,  ".csv", sep = ""))
# Convert data table to a matrix
matTau <-  data.frame(rep = rep(seq(1, dim(setTau)[1]), times = dim(setTau[, -1])[2]), counts = as.vector(unname(as.matrix(setTau[,-1]))), nbin = rep(seq(1, dim(setTau[,-1])[2]), each = dim(setTau)[1]))


# Read the Approximated Tau-Leap simulation data from a CSV file
setApp_Tau <- fread(file=paste(pathInput,"/ApproxTauLeap/approxTau_", chrMonApp, ".csv", sep = ""))
# Convert data table to a matrix
matApp_Tau <-  data.frame(rep = rep(seq(1, dim(setApp_Tau)[1]), times = dim(setApp_Tau[, -1])[2]), counts = as.vector(unname(as.matrix(setApp_Tau[,-1]))), nbin = rep(seq(1, dim(setApp_Tau[,-1])[2]), each = dim(setApp_Tau)[1]))
nApp_Tau <- length(setApp_Tau[1,-1])

if(monApp < 10^13){
# Read the Approximated Gillespie simulation data from a CSV file
  setApp_Gill <- fread(file=paste(pathInput,"/ApproxGillespie/ApproxGill_", chrMonApp, ".csv", sep = ""))
  # Convert data table to a matrix
  matApp_Gill <-  data.frame(rep = rep(seq(1, dim(setApp_Gill)[1]), times = dim(setApp_Gill[, -1])[2]), counts = as.vector(unname(as.matrix(setApp_Gill[,-1]))), nbin = rep(seq(1, dim(setApp_Gill[,-1])[2]), each = dim(setApp_Gill)[1]))
  nApp_Gill <- length(setApp_Gill[1,-1])
}



# Perform a sanity check to ensure the total mass at the end is consistent with the initial mass
N <- mon
for (i in unique(matGill$rep)){
  if (sum(setGill[i,-1] * c(1, 1:(N-1))) != N) {
    stop("Error in the mass conservation for Gillespie simulation!")  
  } 
  if(sum(setTau[i,-1]*c(1,1:(N-1))) != mon){
    stop("Error in the mass conservation for Tau-Leap simulation!")  

  }
} 



########################################
#         BUILD BINNED RESULTS
########################################
# Set initial values for repetitions and bins
rep <- dim(setGill)[1]
a <- ceiling(N^(1/4))
bins <- length(setGill[1,-1])%/% a

# Initialize a numeric vector for binned state in Gillespie and Tau-Leap simulations
binned_stateGill <- data.frame(rep = rep(seq(1, rep), times = bins), counts = 0, nbin = rep(seq(1, bins), each = rep))
binned_stateTau <- data.frame(rep = rep(seq(1, rep), times = bins), counts = 0, nbin = rep(seq(1, bins), each = rep))

# Bin together different states in order to visualize better the results of the Gillespie algorithm 
for (r in 1:rep){
    for (b in 1:(bins-1)){
        # Calculate the sum of particles within the b-th bin
        binned_stateGill$counts[which(binned_stateGill$rep == r & binned_stateGill$nbin == b)] <- sum(matGill$counts[matGill$rep == r][1+(1:a) + (b-1)*a])

        binned_stateTau$counts[which(binned_stateTau$rep == r & binned_stateTau$nbin == b)] <- sum(matTau$counts[matTau$rep == r][1+(1:a) + (b-1)*a])

    }
        
    # Calculate the sum of monomers for the last bin 
    binned_stateGill$counts[which(binned_stateGill$rep == r & binned_stateGill$nbin == bins)] <- sum(matGill$counts[matGill$rep == r][1+(1:a-1) + (bins-1)*a])
    binned_stateTau$counts[which(binned_stateTau$rep == r & binned_stateTau$nbin == bins)] <- sum(matTau$counts[matTau$rep == r][1+(1:a-1) + (bins-1)*a])
}






########################################
#                PLOT 
########################################
# set the path where the plots will be saved
pathOutput <- "nanoparticlesProject/graphics"     # path where the plots will be stored



# ----------------------------------------------------------------
# ---------- original Gillespie's and Tau-Leap's results ---------
# ----------------------------------------------------------------

# Initialize variables to find relevant states in the original simulation
minnGill <- N
maxxGill <- 0
minnTau <- N
maxxTau <- 0

# Loop across the different repetitions to find the minimum and the maximum states  with particles'concentration greater than 0 for original Gillespie and Tau-Leap simulations
for (r in 1:rep){ 
    posintGill <- which(matGill$counts[matGill$rep == r] > 0)
    minnGill <- min(c(posintGill, minnGill))
    maxxGill <- max(c(posintGill, maxxGill))

    posintTau <- which(matTau$counts[matTau$rep == r] > 0)
    minnTau <- min(c(posintTau, minnTau))
    maxxTau <- max(c(posintTau, maxxTau))
}

# Create two vectors of relevant states in the original Gillespie and Tau-Leap simulation respectively
relevant_statesGill <- minnGill:maxxGill
relevant_statesTau  <- minnTau:maxxTau


# Scale the states we have found by sqrt(N), such that we obtain a comparable scale when N grows
scaled_relevant_statesGill <- round(relevant_statesGill / sqrt(N), 1)
saveGill <- which(scaled_relevant_statesGill <= 2.05)
relevant_statesGill <- relevant_statesGill[saveGill]
scaled_relevant_statesGill <- scaled_relevant_statesGill[saveGill]
scaled_relevant_statesTau <- round(relevant_statesTau / sqrt(N), 1)
saveTau <- which(scaled_relevant_statesTau <= 2.05)
relevant_statesTau <- relevant_statesTau[saveTau]
scaled_relevant_statesTau <- scaled_relevant_statesTau[saveTau]


# set the graphical elements
limX <- seq(0, max(relevant_statesGill, relevant_statesTau))
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







# plots for the Tau-Leap's algorithm
TauLeapPLOT <- list()
for(r in 1:rep){
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
    
    pdf(file=paste(pathOutput, "TauLeap10_", chrMon,"-",r , ".pdf", sep = ""), width = 9, height = 5)
    TauLeapPLOT[[r]] <- ggplot(subset(matTau, rep == r & nbin %in% limX), aes(x = factor(nbin), y = counts)) + geom_bar(stat="identity",  color="black", fill="black") + ylim(c(0,max(matGill$counts, matTau$counts)))+ labs(x = labX,  y = labY) + scale_x_discrete(breaks = limX[index], labels =  sprintf("%.1f",labels))  + textThemeOriginal + theme(plot.margin = margin(0.5, 0.5, 0.5, 1, unit = "cm")) + annotate(geom="text", x=xText, y=max(matGill$counts, matTau$counts) - 0.5, label=TeX(lab), color="black", size = 8)
    print(TauLeapPLOT[[r]])
    dev.off()
}


# ----------------------------------------------------------------
# ---------- binned Gillespie's and Tau-Leap's results -----------
# ----------------------------------------------------------------

# Initialize variables to find relevant bins in the binned simulation
minnGill_binned <- bins
maxxGill_binned <- 0
minnTau_binned <- bins
maxxTau_binned <- 0

# Loop across the different repetitions to find the minimum and the maximum states  with particles'concentration greater than 0 for binned Gillespie and Tau-Leap simulations
for (r in 1:rep){ 
    posintGill_binned <- which(binned_stateGill$counts[binned_stateGill$rep == r]>0)
    minnGill_binned <- min(c(posintGill_binned, minnGill_binned))
    maxxGill_binned <- max(c(posintGill_binned, maxxGill_binned))

    posintTau_binned <- which(binned_stateTau$counts[binned_stateTau$rep == r]>0)
    minnTau_binned <- min(c(posintTau_binned, minnTau_binned))
    maxxTau_binned <- max(c(posintTau_binned, maxxTau_binned))
}

# Create two vectors of relevant states in the original Gillespie and Tau-Leap simulation respectively
relevant_statesGill_binned <- minnGill_binned:maxxGill_binned
relevant_statesTau_binned  <- minnTau_binned:maxxGill_binned

# Scale the states we have found by sqrt(N)/d, such that we obtain a comparable scale when N grows
scaled_relevant_statesGill_binned <- round(relevant_statesGill_binned*a / sqrt(N), 1)
scaled_relevant_statesTau_binned  <- round(relevant_statesTau_binned* a / sqrt(N), 1)
saveGill <- which(scaled_relevant_statesGill_binned <= 2.05)
relevant_statesGill_binned <- relevant_statesGill_binned[saveGill]
scaled_relevant_statesGill_binned <- scaled_relevant_statesGill_binned[saveGill]
saveTau <- which(scaled_relevant_statesTau_binned <= 2.05)
relevant_statesTau_binned <- relevant_statesTau_binned[saveTau]
scaled_relevant_statesTau_binned <- scaled_relevant_statesTau_binned[saveTau]

# graphical parameters
limX_binned <- seq(min(relevant_statesGill_binned, relevant_statesTau_binned), max(relevant_statesGill_binned, relevant_statesTau_binned))
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
      yText <- max(binned_stateGill$counts, binned_stateTau$counts) - 0.5
    }else if(chrMon == "4"){
      xText <- 7
      yText <- max(binned_stateGill$counts, binned_stateTau$counts) - 0.5
    }else if(chrMon =="8"){
      xText <- 70
      yText <- 115
    }

    lab <- sprintf("Binned process, $N = 10^{%d}, rep = %d$",chrMonInt, r)
    labels <- sprintf("%.1f",labels_binned)
    labels[1] <- ""

    gillespiePLOT_binned[[r]] <- ggplot(subset(binned_stateGill, rep == r & nbin %in% limX_binned), aes(x = factor(nbin), y = counts)) + geom_bar(stat="identity",  color="black", fill="black") + ylim(c(0,max(binned_stateGill$counts, binned_stateTau$counts)))+ labs(x =  labX,  y = labY) + scale_x_discrete(breaks = limX_binned[index_binned], labels =  labels)  + textThemeOriginal + theme(plot.margin = margin(0, 0.5, 0, 0.5, unit = "cm")) +  annotate(geom="text", x=xText, y=yText, label=TeX(lab), color="black", size = 8)
    print(gillespiePLOT_binned[[r]])
    dev.off()
}
# ------------------------------------------
# ------------------------------------------


# plots of the results of the Tau-Leap's algorithm after the binning
tauPLOT_binned <- list()
for(r in 1:rep){
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
      yText <- max(binned_stateGill$counts, binned_stateTau$counts) - 0.5
    }else if(chrMon == "4"){
      xText <- 7
      yText <- max(binned_stateGill$counts, binned_stateTau$counts) - 0.5
    }else if(chrMon =="8"){
      xText <- 70
      yText <- 115
    }

    lab <- sprintf("Binned process, $N = 10^{%d}, rep = %d$",chrMonInt, r)
    labels <- sprintf("%.1f",labels_binned)
    labels[1] <- ""

    pdf(file=paste(pathOutput, "TauLeapBinned10_", chrMon,"-",r , ".pdf", sep = ""), width = 9, height = 5)
    tauPLOT_binned[[r]] <- ggplot(subset(binned_stateTau, rep == r & nbin %in% limX_binned), aes(x = factor(nbin), y = counts)) + geom_bar(stat="identity",  color="black", fill="black") + ylim(c(0,max(binned_stateGill$counts, binned_stateTau$counts)))+ labs(x =  labX,  y = labY) +  scale_x_discrete(breaks = limX_binned[index_binned], labels =  labels)   + textThemeOriginal + theme(plot.margin = margin(0.5, 0.5, 0.5, 1, unit = "cm")) + annotate(geom="text", x=xText, y=yText, label=TeX(lab), color="black", size = 8)
    print(tauPLOT_binned[[r]])
    dev.off()
}




# ----------------------------------------------------------------
# ------ Comparison Binned results obtained with Tau-Leap and ----
# ------ Approximated process obtained with Tau-Leap -------------
# ----------------------------------------------------------------

repApp <- dim(setApp_Tau)[1]
aApp <- ceiling(monApp^(1/4))
binsApp <- length(setApp_Tau[1,-1])

# Initialize variables to find relevant bins in the approximated process
minnApp_Tau <- length(setApp_Tau[1,-1])
maxxApp_Tau <- 0
# Loop across the different repetitions to find the minimum and the maximum bin with particles'concentration greater than 0 for approximated process
for (r in 1:repApp){ 
    posintApp <- which(matApp_Tau$counts[matApp_Tau$rep == r]>0)
    minnApp_Tau <- min(c(posintApp, minnApp_Tau))
    maxxApp_Tau <- max(c(posintApp, maxxApp_Tau))
}
relevant_statesApp_Tau <- minnApp_Tau:maxxApp_Tau
new_relevant_statesApp_Tau <- round(relevant_statesApp_Tau*aApp / sqrt(monApp),3)


# Compute the frequencies of each bin for the binned process
freqBinned_stateGill <- binned_stateGill
for(r in 1:rep){
  freqBinned_stateGill[freqBinned_stateGill$rep == r, "counts"] <- freqBinned_stateGill[freqBinned_stateGill$rep == r, "counts"]/sum(freqBinned_stateGill[freqBinned_stateGill$rep == r, "counts"]) 
}


# ------------- FIGURE 3 (row 1 and 2) OF THE PAPER --------------
# plots with binned results of the exact process and results of the approximated process
if(monApp >= 10^13){
    for(r in 1:rep){
      print(r)
      pdf(file=paste(pathOutput, "cfr_", chrMon,"mon-",chrMonApp ,"monApp-",r,  ".pdf", sep = ""), width = 9.5, height = 5)
      high <- freqBinned_stateGill[freqBinned_stateGill$rep == r & freqBinned_stateGill$nbin %in% relevant_statesGill_binned, "counts"]
      # scale the frequencies of Tau-Leap in order to obtain comparable bin dimension
      highApp_Tau <- matApp_Tau[matApp_Tau$rep == r & matApp_Tau$nbin %in% relevant_statesApp_Tau, "counts"]/sum( matApp_Tau[matApp_Tau$rep == r & matApp_Tau$nbin %in% relevant_statesApp_Tau,"counts"])/bins*binsApp

      df <- c(high, highApp_Tau)
      df <- data.frame(h = df, col = NA)
      df$col <- as.factor(c(rep("1", length(high)), rep("2", length(highApp_Tau))))
      df$states <- c(relevant_statesGill_binned*a / sqrt(N), relevant_statesApp_Tau*aApp / sqrt(monApp)) 

      if(chrMonApp == "4"){
        labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"4"))
      }else if(chrMonApp == "5"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"5"))
      }else if(chrMonApp == "6"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"6"))
      }else if(chrMonApp == "7"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"7"))
      }else if(chrMonApp == "8"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"8"))
      }else if(chrMonApp == "9"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"9"))
      }else if(chrMonApp == "10"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"10"))
      }else if(chrMonApp == "11"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"11"))
      }else if(chrMonApp == "12"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"12"))
      }else if(chrMonApp == "13"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"13"))
      }else if(chrMonApp == "14"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"14"))
      }else if(chrMonApp == "15"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"15"))
      }

      if(chrMon == "4"){
        labExact <- expression("N = "~10^4)
      }else if(chrMon == "5"){
        labExact <- expression("N = "~10^5)
      }else if(chrMon == "6"){
        labExact <- expression("N = "~10^6)
      }else if(chrMon == "7"){
        labExact <- expression("N = "~10^7)
      }else if(chrMon == "8"){
        labExact <- expression("N = "~10^8)
      }

      g1 <- ggplot(data = df, aes(x = states, y = h)) +
      geom_bar(data = filter(df, col == "1"), aes(x = states, y = h, fill = "binned"), stat = "identity") +
      scale_fill_manual(values = c("binned" = "gray55"), labels = labExact) + 
      xlim(0, 2) + 
      labs(
        x = expression(paste("Scaled nanoparticle size (in the unit of" ~sqrt(N)~ ")", sep = "")), 
        y = "Frequency", 
        color = NULL, 
        fill = "Binned original process"
      ) 
      

      g1 <- g1 + ggnewscale::new_scale_fill()  + 
      geom_line(data = filter(df, col == "2"), aes(x = states, y = h, color = "tauApp"), linewidth = 1.5) +  
      scale_color_manual(
        values = c(alpha("red", 0.8)), 
        labels =  expression("Tau-leap")
      ) +    
      labs(color = labApprox) + 
      theme(
        legend.position = c(0.27, 0.7),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, family = "", face = "bold")
      ) + 
      textThemeApprox + 
      theme(
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
        axis.text.y = element_blank(),
        legend.box = "vertical", 
        legend.background = element_blank(),
      ) 

      print(g1)

      dev.off()

  }

}
# ---------------------------------------------------


# -------------------------------------------------------------------------------
# ------ Comparison Binned results obtained with Tau-Leap  ----------------------
# ------ Approximated process obtained with Tau-Leap and with Gillespie ---------
# -------------------------------------------------------------------------------


# ------------- FIGURE 3 (row 3) OF THE PAPER --------------
if(monApp < 10^13){
  # Initialize variables to find relevant bins in the approximated process
  minnApp_Gill <- length(setApp_Gill[1,-1])
  maxxApp_Gill <- 0
  # Loop across the different repetitions to find the minimum and the maximum bin with particles'concentration greater than 0 for approximated process
  for (r in 1:repApp){ 
      posintApp <- which(matApp_Gill$counts[matApp_Gill$rep == r]>0)
      minnApp_Gill <- min(c(posintApp, minnApp_Gill))
      maxxApp_Gill <- max(c(posintApp, maxxApp_Gill))
  }
  relevant_statesApp_Gill <- minnApp_Gill:maxxApp_Gill
  new_relevant_statesApp_Gill <- round(relevant_statesApp_Gill*aApp / sqrt(monApp),3)

  # plots with binned results of the exact process and results of the approximated process (both Gillespie and Tau-Leap)
  for(r in 1:rep){
      print(r)
      pdf(file=paste(pathOutput, "cfr_APPgillTau", chrMonApp,"monApp-r",r,  ".pdf", sep = ""), width = 10, height = 5)

      
      high <- freqBinned_stateGill[freqBinned_stateGill$rep == r & freqBinned_stateGill$nbin %in% relevant_statesGill_binned, "counts"]
      highApp_Gill <- matApp_Gill[matApp_Gill$rep == r & matApp_Gill$nbin %in% relevant_statesApp_Gill, "counts"]/sum( matApp_Gill[matApp_Gill$rep == r & matApp_Gill$nbin %in% relevant_statesApp_Gill,"counts"])/bins*binsApp
      highApp_Tau <- matApp_Tau[matApp_Tau$rep == r & matApp_Tau$nbin %in% relevant_statesApp_Tau, "counts"]/sum( matApp_Tau[matApp_Tau$rep == r & matApp_Tau$nbin %in% relevant_statesApp_Tau,"counts"])/bins*binsApp

      df <- c(high, highApp_Gill, highApp_Tau)
      df <- data.frame(h = df, col = NA)
      df$col <- as.factor(c(rep("1", length(high)), rep("2", length(highApp_Gill)),  rep("3", length(highApp_Tau))))
      df$states <- c(relevant_statesGill_binned*a / sqrt(N), relevant_statesApp_Gill*aApp / sqrt(monApp), relevant_statesApp_Tau*aApp / sqrt(monApp)) 

      if(chrMonApp == "4"){
        labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"4"))
      }else if(chrMonApp == "5"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"5"))
      }else if(chrMonApp == "6"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"6"))
      }else if(chrMonApp == "7"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"7"))
      }else if(chrMonApp == "8"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"8"))
      }else if(chrMonApp == "9"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"9"))
      }else if(chrMonApp == "10"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"10"))
      }else if(chrMonApp == "11"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"11"))
      }else if(chrMonApp == "12"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"12"))
      }else if(chrMonApp == "13"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"13"))
      }else if(chrMonApp == "14"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"14"))
      }else if(chrMonApp == "15"){
          labApprox <- expression(bold(Approximated~process~N~"\u003d"~"10"^"15"))
      }

      if(chrMon == "4"){
        labExact <- expression("N = "~10^4)
      }else if(chrMon == "5"){
        labExact <- expression("N = "~10^5)
      }else if(chrMon == "6"){
        labExact <- expression("N = "~10^6)
      }else if(chrMon == "7"){
        labExact <- expression("N = "~10^7)
      }else if(chrMon == "8"){
        labExact <- expression("N = "~10^8)
      }


      g1 <- ggplot(data = df, aes(x = states, y = h)) + 
      geom_bar(data = filter(df, col == "1"), aes(x = states, y = h, fill = "binned"), stat = "identity") + 
      scale_fill_manual(values = c("binned" = "gray55"), labels = labExact) + 
      xlim(0, 2) +
      labs(
        x = expression(paste("Scaled nanoparticle size (in the unit of " ~ sqrt(N) ~ ")", sep = "")),
        y = "Frequency",
        color = NULL,
        fill = "Binned original process"
      ) 



      g1 <- g1 +  ggnewscale::new_scale_fill()   +  
      geom_line(data = filter(df, col == "2"), aes(x = states, y = h, color = "gillespieApp"), linewidth = 1.5) +
      geom_line(data = filter(df, col == "3"), aes(x = states, y = h, color = "tauApp"), linewidth = 1.5) +   
      scale_color_manual(
        values = c(alpha("blue", 0.8), alpha("red", 0.8)),
        labels = c(expression("Gillespie"), expression("Tau-leap"))
      ) +
      # labs(color = expression(paste("Approximated process, N=" ~ 10^{8}, sep = "")))
      labs(color = labApprox)    +
      theme(
        legend.position = c(0.255, 0.7),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, family = "", face = "bold")
      ) +
      textThemeApprox +
      theme(
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
        axis.text.y = element_blank(),
        legend.box = "vertical", 
        legend.background = element_blank(),
      ) 

    print(g1)
    dev.off()  
      

  }


}
# ---------------------------------------------------
