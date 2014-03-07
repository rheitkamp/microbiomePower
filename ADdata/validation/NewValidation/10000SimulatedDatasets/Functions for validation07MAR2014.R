################################
###Simulated Dataset Creation###
################################

###
Premainder <- function(x) {
  rowsize <- dim(x)[1]
  colsize <- dim(x)[2]
  dims <- dimnames(x)[[2]]
  y <- matrix(nrow=rowsize, ncol=colsize)
  for(i in 1:rowsize){
    total <- apply(X=x,MARGIN=1,sum)[i]
    y[i,1] <- x[i,1]/total ### can possible be rewritten as x[i,1]/total
    total <- total - x[i,1]
    for(k in 2:colsize){
      y[i,k] <- x[i,k]/ total
      total <- total - x[i,k]
    }
  } 
  return(matrix(y,nrow=colsize, ncol=rowsize, dimnames=list(dims, "Mean")))
}


###
getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

###
collectParameters <- function(datamatrix) {
  taxa <- rownames(datamatrix)
  parameterNames <- list("type","a","b","location","scale")
  BetaParameters <- matrix(data = NA, nrow = length(taxa), ncol = length(parameterNames), dimnames=list(taxa, parameterNames))
  
  for (i in 1:nrow(BetaParameters)) {
    parameters <- getBetaParams(datamatrix[i,1],datamatrix[i,2]);
    BetaParameters[i,] <- c(parameters$type, parameters$a, parameters$b, parameters$location, parameters$scale)
  }
  
  return(BetaParameters)
}

###
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

###
simulateBrokenStick <- function(rawData, numberSubjects=25) {
  #set.seed(1234) #always use the same seed during testing
  library(PearsonDS) #use pearsons distribution library (beta distribution)
  library(HMP) #use HMP package
  
  #(rawData <- Data) #read data from file
  
  #rowSums(rawData) ##!!validate that column sums equal 1!!
  
  #get standard deviation of each taxon in the data set
  (rawStdDiv <- matrix(apply(rawData,2,sd),dimnames=list(colnames(rawData), "SD"))) 
  
  #get arithmetic mean of each taxon in the dataset
  (rawMean <- matrix(colMeans(rawData), nrow=1, ncol=dim(rawData)[2], dimnames=list("Mean",colnames(rawData))))
  
  (meanPercentRemainder <- Premainder(rawMean)) #get percent remainder for mean proportions of taxa
  
  (pR_MeanStDev <- cbind(meanPercentRemainder, rawStdDiv)) #combine percent remainder means with raw stdev
  
  parameters <- collectParameters(pR_MeanStDev) #calculate parameters of percent remainder beta dist
  
  brokenStickSim <- spaceFill(pR_MeanStDev,parameters,numberSubjects) #run the simulation
  
  #Barchart.data(brokenStickSim, title=outputLabel) #display data as barchart
  
  #print('Standard Deviation, simulated vs. provided data')
  #print(apply(brokenStickSim,2,sd)) #compare simulated standard deviations
  #print(apply(rawData,2,sd)) #to raw reported standard deviations
  
  #print('Mean, simulated vs. provided data')
  #print(colMeans(brokenStickSim)) #compare simulated means
  #print(colMeans(rawData)) #to raw reported means
  #print('')
  #print('')
  return(brokenStickSim)
}

###
validationdatasets <- function(rawData, numberSubjects, numberSimulations){
  listdataset <- replicate(n=length(numberSimulations), expr=list())
  for(i in 1:numberSimulations){
    listdataset[[i]] <- round((simulateBrokenStick(rawData, numberSubjects))*3000)
  }
  return(listdataset)
}

##################################
###Hypothesis Testing Functions###
##################################

meanCheck <- function(simDatalist, observedData, numReads){
  observedData <- round(observedData*numReads)
  fit <- DM.MoM(observedData)
  results <- Xmc.sevsample(group.data=simDatalist, pi0=fit$pi)
  return(results)
}

###
getthetaValues <- function(x){
  numlist <- length(x)
  thetaValues <- numeric(length=numlist)
  for(i in 1:numlist){
    thetaValues[i] <- DM.MoM(x[[i]])$theta
  }
  return(thetaValues)
}

###
thetaTest2tailed <- function(simDatalist, observedData, numReads){
  thetalist <- getthetaValues(simDatalist)
  observedData <- round(observedData*numReads)
  observed_theta <- DM.MoM(observedData)$theta
  pvalue <- (length(which(thetalist > observed_theta))/length(thetalist))/2
  pvalue <- list(pvalue)
  names(pvalue) <- c("p value")
  return(pvalue)
}


