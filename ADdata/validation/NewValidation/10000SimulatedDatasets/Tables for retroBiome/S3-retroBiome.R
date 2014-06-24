#############################################
#                 retroBiome                #
#                                           #
# An R script for simulating microbiome     #
# samples from historical community data    #
# provided as mean proportions and standard #
# deviation.                                #
#                                           #
# By Sylvia Cheng and Rae A Heitkamp        #
# Version 0.0 - 23JUN2014                   #
# sylvia.cheng07@gmail.com                  #
# rae.heitkamp@gmail.com                    #
#############################################

  #              FUNCTIONS               #

#############################################
# PREMAINDER:CALCULATE SPACE FILLED BY TAXA # 
#   Input:  table containing rows of        #
#           taxa with arithmetic mean       #
#           proportion in descending order  #
#           and standard deviation of       #
#           proportions                     #
#   Output: table containing rows           #
#           of taxa with proportion of      #
#           available space filled and      #
#           standard deviation of original  #
#           proportions                     #    
#############################################

Premainder <- function(taxaMeanAndSD) { 
  numberOfTaxa <- dim(taxaMeanAndSD)[1] 
  taxa <- rownames(taxaMeanAndSD) 
  scaledMean <- matrix(nrow=numberOfTaxa, dimnames=list(taxa,"Mean")) 
  SD <- matrix(taxaMeanAndSD[,2], dimnames=list(taxa,"SD")) 
  for(i in 1:1){
    availableSpace <- apply(X=taxaMeanAndSD,MARGIN=2,sum)[i]
    scaledMean[i,1] <- taxaMeanAndSD[i,1]/availableSpace
    availableSpace <- availableSpace - taxaMeanAndSD[i,1] 
    for(k in 2:numberOfTaxa){ 
      scaledMean[k,i] <- taxaMeanAndSD[k,i]/ availableSpace 
      availableSpace <- availableSpace - taxaMeanAndSD[k,i] 
    }
  } 
  return(cbind(scaledMean,SD))
}

#############################################
# GETBETAPARAMS:DERIVE BETA PARAMETERS      #
#   Input:  an arithmetic mean and standard #
#           deviation                       #
#   Output: a list of type-1 beta function  #
#           parameters with location set to #
#           0 and scale set to 1, and alpha #
#           and beta values as derived from #
#           input mean and standard         #
#           deviation                       # 
#############################################

getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

#############################################
# COLLECTPARAMETERS:COMBINE PARAMETER LISTS #
#   Input:  a matrix with rows of taxa with #
#           proportions of space filled and #
#           standard deviations             #
#   Output: a matrix of the type-1 beta     #
#           function parameters for each    #
#           taxon                           #
#############################################

collectParameters <- function(data) {
  taxa <- rownames(data)
  parameterNames <- list("type","a","b","location","scale")
  BetaParameters <- matrix(data = NA, nrow = length(taxa), ncol = length(parameterNames), dimnames=list(taxa, parameterNames))
  
  for (i in 1:nrow(BetaParameters)) {
    parameters <- getBetaParams(data[i,1],data[i,2]);
    BetaParameters[i,] <- c(parameters$type, parameters$a, parameters$b, parameters$location, parameters$scale)
  }
  
  return(BetaParameters)
}

#############################################
# SPACEFILL:SIMULATE SAMPLES FROM PARAMETERS#
#   Input:  *matrix of proportion space     #
#            filled by each taxon           #
#           *matrix of beta distribution    #
#            parameters by taxon            #
#           *number of samples to simulate  #
#   Output: a matrix of simulated microbiome#
#           samples                         #
#############################################

spaceFill <- function (percentSpaceFillTable, distributionParameters, numberOfSamplesToSimulate){ 
  numberOfTaxa <- dim(percentSpaceFillTable)[1]
  
  outputTable <- matrix(nrow=numberOfSamplesToSimulate,ncol=numberOfTaxa)
  taxaNames <- dimnames(percentSpaceFillTable)[[1]]
  colnames(outputTable) <- taxaNames #list taxa names in columns
  rownames(outputTable) <- rownames(outputTable, do.NULL= FALSE, prefix= "Sample") #call rows samples
  
  for (i in 1:numberOfSamplesToSimulate){
    availableSpace <- 1 
    for (j in 1:(numberOfTaxa-1)){
      spaceTaken <- rpearson(n=1, params=distributionParameters[j,]); 
      outputTable[i,j] <- availableSpace * spaceTaken
      availableSpace <- availableSpace - outputTable[i,j];
    }
    outputTable[i,numberOfTaxa] <- availableSpace;
  }
  return(outputTable)
}

#############################################
 #                 WRAPPER                 #
#############################################
# To execute retroBiome on your data table, #
# uncomment the last lines of this file. Be #
# sure your data is in the proper format as #
# the script does not perform checks.       #
#                                           #
# To uncomment, remove the number sign (#)  #
# from the beginning of the line of code    #
#############################################

retroBiome <- function(summaryTable,outputLabel,numberOfSamplesToSimulate=25) {
  #set.seed(1234) #always use the same seed during testing
  library(PearsonDS) 
  library(HMP) 
  
  percentSpaceFillTable <- Premainder(summaryTable) 
  
  distributionParameters <- collectParameters(percentSpaceFillTable) 
  
  sampleSimulation <- spaceFill(percentSpaceFillTable,distributionParameters,numberOfSamplesToSimulate) 
  
  #Display the simulated data visually
  Barchart.data(sampleSimulation, title=outputLabel)
  
  #Compare the simulated and experimental means
  Simulated.Mean <- colMeans(sampleSimulation)
  Provided.Mean <- array(summaryTable[,1])
  compareSimulatedWithProvidedMean <- cbind(Simulated.Mean, Provided.Mean)  
  
  #Compare the simulated and experimental standard deviations
  Simulated.SD <- apply(sampleSimulation,2,sd) 
  Provided.SD <- array(summaryTable[,2]) 
  compareSimulatedWithProvidedSD <- cbind(Simulated.SD, Provided.SD) 
  
  #Print the comparisons to the console
  listMeanAndSD <- list(compareSimulatedWithProvidedMean, compareSimulatedWithProvidedSD)  
  names(listMeanAndSD) <- c('Mean, simulated vs. provided data', 'Standard Deviation, simulated vs. provided data')
  print(listMeanAndSD)

  return(sampleSimulation)
}

inputData <- read.table("test.txt") #replace "example.txt" with your tab-delineated table of proportion and standard deviation
inputData2 <- read.table("ExampleMeanSD.txt")
outputLabel <- "Simulated Samples" #replace "Simulated Samples" with the title of your graph of simulated samples
retroBiome(inputData, outputLabel, numberOfSamplesToSimulate=25) #if you want to simulate more than 25 samples, change "numberOfSamplesToSimulate" to the number 