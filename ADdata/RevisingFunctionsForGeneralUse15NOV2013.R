#CALCULATE PERCENT REMAINDER OF TAXON MEANS ########Fixed this!!!
###########################################

Premainder <- function(x) {
  rowsize <- dim(x)[1]
  dims <- rownames(x)
  y <- matrix(nrow=rowsize, dimnames=list(dims,"Mean"))
  SD <- matrix(x[,2], dimnames=list(rownames(x),"SD"))
  for(i in 1:1){
    total <- apply(X=x,MARGIN=2,sum)[i]
    y[i,1] <- x[i,1]/total
    total <- total - x[i,1]
    for(k in 2:rowsize){
      y[k,i] <- x[k,i]/ total
      total <- total - x[k,i]
    }
  } 
  z <- cbind(y,SD)
  return(z)
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

#SIMULATE SUBJECTS WITH BROKEN-STICK MODEL FOR A GIVEN DISTRIBUTION
###################################################################

spaceFill <- function (dataMatrix, distParameters, subjects){
  subjects <- subjects #number of subjects for rows
  taxa <- dim(dataMatrix)[1] #number of taxa for columns 
  
  Cdata <- matrix(nrow=subjects,ncol=taxa)
  taxaNames <- dimnames(dataMatrix)[[1]]
  colnames(Cdata) <- taxaNames #list taxa names in columns
  rownames(Cdata) <- rownames(Cdata, do.NULL= FALSE, prefix= "Sample") #call rows samples
  
  for (i in 1:subjects){
    total <- 1
    for (j in 1:(taxa-1)){
      r <- rpearson(n=1, params=distParameters[j,]);
      Cdata[i,j] <- total * r
      total <- total - Cdata[i,j];
    }
    Cdata[i,taxa] <- total;
  }
  return(Cdata)
}


#WRAPPER FUNCTION
#################
simulateBrokenStick <- function(data,outputLabel,numberSubjects=25) {
  #set.seed(1234) #always use the same seed during testing
  library(PearsonDS) #use pearsons distribution library (beta distribution)
  library(HMP) #use HMP package
  
  (data <- read.table(inputFilename)) #read data from file
  
  prdata <- Premainder(data) #calculate percent remainder of taxon means
  
  parameters <- collectParameters(prdata) #calculate parameters of percent remainder beta distribution
  
  brokenStickSim <- spaceFill(prdata,parameters,numberSubjects) #run the simulation
  
  Barchart.data(brokenStickSim, title=outputLabel) #display data as barchart
  
  Simulated.Mean <- t(t(colMeans(brokenStickSim))) # calculate the simulated mean
  Provided.Mean <- array(data[1]) # grabs the provided mean
  SimVSProMean <- cbind(Simulated.Mean, Provided.Mean) #place the simulated and provided mean in a table 
  
  Simulated.SD <- t(t(apply(brokenStickSim,2,sd))) #calculates the simulated standard deviation
  Provided.SD <- array(data[2]) #grabs the provided standard deviation
  SimVSProSD <- cbind(Simulated.SD, Provided.SD) #place the simulated and provided standard deviations in a table
  
  listMSD <- list(SimVSProMean, SimVSProSD) #place comparisons in a list 
  names(listMSD) <- c('Mean, simulated vs. provided data', 'Standard Deviation, simulated vs. provided data')
  print(listMSD)

  return(brokenStickSim)
}