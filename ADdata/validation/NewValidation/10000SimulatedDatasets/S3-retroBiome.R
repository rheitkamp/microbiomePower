#CALCULATE PERCENT REMAINDER OF TAXON MEANS ########Fixed this!!!
###########################################

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

#CALCULATE PARAMETERS OF BETA DISTRIBUTION
##########################################

getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

#CALCULATE BETA PARAMETERS FOR EACH TAXON IN DATASET
####################################################

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

#SIMULATE SAMPLES WITH A SPACE-FILLING MODEL FOR A GIVEN DISTRIBUTION
#####################################################################

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


#WRAPPER FUNCTION
#################
retroBiome <- function(summaryTable,outputLabel,numberOfSamplesToSimulate=25) {
  #set.seed(1234) #always use the same seed during testing
  library(PearsonDS) #use pearsons distribution library (beta distribution)
  library(HMP) #use HMP package
  
  #(summaryTable <- read.table(inputFilename)) #read data from file
  
  percentSpaceFillTable <- Premainder(summaryTable) #calculate percent remainder of taxon means
  
  distributionParameters <- collectParameters(percentSpaceFillTable) #calculate parameters of percent remainder beta distribution
  
  sampleSimulation <- spaceFill(percentSpaceFillTable,distributionParameters,numberOfSamplesToSimulate) #run the simulation
  
  Barchart.data(sampleSimulation, title=outputLabel) #display data as barchart
  
  Simulated.Mean <- colMeans(sampleSimulation) # calculate the simulated mean
  Provided.Mean <- array(summaryTable[,1]) # grabs the provided mean
  compareSimulatedWithProvidedMean <- cbind(Simulated.Mean, Provided.Mean) #place the simulated and provided mean in a table 
  
  Simulated.SD <- apply(sampleSimulation,2,sd) #calculates the simulated standard deviation
  Provided.SD <- array(summaryTable[,2]) #grabs the provided standard deviation
  compareSimulatedWithProvidedSD <- cbind(Simulated.SD, Provided.SD) #place the simulated and provided standard deviations in a table
  
  listMeanAndSD <- list(compareSimulatedWithProvidedMean, compareSimulatedWithProvidedSD) #place comparisons in a list 
  names(listMeanAndSD) <- c('Mean, simulated vs. provided data', 'Standard Deviation, simulated vs. provided data')
  print(listMeanAndSD)

  return(sampleSimulation)
}
