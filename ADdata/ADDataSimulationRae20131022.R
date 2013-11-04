
#############
# FUNCTIONS #
#############

#CALCULATE BETA PARAMETERS FOR EACH TAXON IN DATASET
####################################################

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

#CALCULATE PERCENT REMAINDER OF TAXON MEANS
###########################################

Premainder <- function(x) {
  rowsize <- dim(x)[1]
  colsize <- dim(x)[2]
  dims <- dimnames(x)[[2]]
  y <- matrix(nrow=rowsize, ncol=colsize)
  for(i in 1:rowsize){
    total <- apply(X=x,MARGIN=1,sum)[i]
    y[i,1] <- x[i,1]/100
    total <- total - x[i,1]
    for(k in 2:colsize){
      y[i,k] <- x[i,k]/ total
      total <- total - x[i,k]
    }
  } 
  return(matrix(y,nrow=colsize, ncol=rowsize, dimnames=list(dims, "Mean")))
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
simulateBrokenStick <- function(inputFilename,outputLabel,numberSubjects=25) {
  #set.seed(1234) #always use the same seed during testing
  library(PearsonDS) #use pearsons distribution library (beta distribution)
  library(HMP) #use HMP package
  
  (rawData <- read.table(inputFilename)) #read data from file
  
  colSums(rawData) ##!!validate that column sums equal 1!!
  
  #get standard deviation of each taxon in the data set
  (rawStdDiv <- matrix(apply(rawData,1,sd),dimnames=list(rownames(rawData), "SD"))) 
  
  #get arithmetic mean of each taxon in the dataset
  (rawMean <- matrix(rowMeans(rawData), nrow=1, ncol=dim(rawData)[1], dimnames=list("Mean", rownames(rawData))))

  (meanPercentRemainder <- Premainder(rawMean*100)) #get percent remainder for mean proportions of taxa
  
  (pR_MeanStDev <- cbind(meanPercentRemainder, rawStdDiv)) #combine percent remainder means with raw stdev
  
  parameters <- collectParameters(pR_MeanStDev) #calculate parameters of percent remainder beta dist
  
  brokenStickSim <- spaceFill(pR_MeanStDev,parameters,numberSubjects) #run the simulation
  
  Barchart.data(brokenStickSim, title=outputLabel) #display data as barchart
  
  print('Standard Deviation, simulated vs. provided data')
  print(apply(brokenStickSim,2,sd)) #compare simulated standard deviations
  print(apply(rawData,1,sd)) #to raw reported standard deviations
  
  print('Mean, simulated vs. provided data')
  print(colMeans(brokenStickSim)) #compare simulated means
  print(rowMeans(rawData)) #to raw reported means
  print('')
  print('')
  return(brokenStickSim)
}

#simulateBrokenStick("ADControls.txt","Controls", 5000)
#simulateBrokenStick("ADBaseline.txt","ADBaseline", 300)
#simulateBrokenStick("ADFlareNT.txt","ADFlare", 300)
#simulateBrokenStick("ADFlareT.txt","ADTreatment", 300)
#simulateBrokenStick("ADPostFlare.txt","ADPost")


######GET VALIDATION DATASETS#############
validationdatasets <- function(simdata, numberSimulations, publishNumbersubjects, numReads){ #7.Number of datasets to create/number of times to repeat
  
  BSData <- simdata #1.Simulate 10,000 samples
  
  listdataset <- replicate(n=length(numberSimulations), expr=list()) #2.Creates a place to store everything 
  for(k in 1:numberSimulations){
    publishNumbersubjects <- 22
    col <- dim(BSData)[2]
    table <- matrix(nrow=publishNumbersubjects,ncol=col)
    colnames(table) <- dimnames(BSData)[[2]]
    rownames(table) <- rownames(table, do.NULL= FALSE, prefix= "Sample")
    pick <- sample(x=1:dim(BSData)[1],size=publishNumbersubjects) #3.Pick n sets of samples  
    for(i in 1:publishNumbersubjects){ #4.Put together as one dataset  
      table[i,] <- round(BSData[pick[i],]*numReads)
    }
    listdataset[[k]] <- table #Save the pi stat with dataset
  }
  return(listdataset)
}


#######VARIANCE CHECK, ONE SIM DATASET VS. PUBLISH DATASET AT A TIME#######
variancecheck <- function(inputFilename, datalist, numReads){
  R <- length(datalist)
  pvalues <- numeric(R)
  pubdata <- round(t(read.table(inputFilename))*numReads)
  for(i in 1:R){
    listdata <- list(datalist[[i]],pubdata)
    check <- Xoc.sevsample(listdata)
    pvalues[i] <- check$'p value' 
  }
  return(pvalues)
}

#######MEAN CHECK, ONE SIM DATASET VS. PUBLISH DATASET AT A TIME#######
meancheck <- function(inputFilename, datalist, numReads){
  R <- length(datalist)
  pvalues <- numeric(R)
  pubdata <- round(t(read.table(inputFilename))*numReads)
  for(i in 1:R){
    fit.simdata <- dirmult(datalist[[i]])
    check <- Xsc.onesample(pubdata,fit.simdata$pi)
    pvalues[i] <- check$'p value' 
  }
  return(pvalues)
}

