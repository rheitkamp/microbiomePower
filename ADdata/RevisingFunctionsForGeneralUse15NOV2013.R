#CALCULATE PERCENT REMAINDER OF TAXON MEANS
###########################################

Premainder <- function(x) {
  rowsize <- dim(x)[1]
  colsize <- dim(x)[2]
  dims <- dimnames(x)[[1]]
  total0 <- apply(X=x, MARGIN=2, sum)
  y <- matrix(nrow=rowsize, ncol=colsize)
  for(i in 1:colsize){
    total <- apply(X=x,MARGIN=2,sum)[i]
    y[i,1] <- x[i,1]/total0
    total <- total - x[i,1]
    for(k in 2:rowsize){
      y[k,i] <- x[k,i]/ total
      total <- total - x[k,i]
    }
  } 
  return(matrix(y,nrow=rowsize, ncol=colsize, dimnames=list(dims, "Mean")))
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
  
  #(rawData <- read.table(inputFilename)) #read data from file
  
  rowsums(rawData[,1]) ##!!validate that column sums equal 1!!
  
  #get standard deviation of each taxon in the data set
  (rawStdDiv <- matrix(rawData[,2],dimnames=list(rownames(rawData), "SD"))) 
  
  #get arithmetic mean of each taxon in the dataset
  (rawMean <- matrix(rawData[,1], nrow=1, ncol=dim(rawData)[1], dimnames=list("Mean", rownames(rawData))))
  
  (meanPercentRemainder <- Premainder(rawMean*100)) #get percent remainder for mean proportions of taxa
  
  (pR_MeanStDev <- cbind(meanPercentRemainder, rawStdDiv)) #combine percent remainder means with raw stdev
  
  parameters <- collectParameters(pR_MeanStDev) #calculate parameters of percent remainder beta dist
  
  brokenStickSim <- spaceFill(pR_MeanStDev,parameters,numberSubjects) #run the simulation
  
  Barchart.data(brokenStickSim, title=outputLabel) #display data as barchart
  
  print('Standard Deviation, simulated vs. provided data')
  print(apply(brokenStickSim,2,sd)) #compare simulated standard deviations
  print(apply(rawData[,2])) #to raw reported standard deviations
  
  print('Mean, simulated vs. provided data')
  print(colMeans(brokenStickSim)) #compare simulated means
  print(rawData[,1]) #to raw reported means
  print('')
  print('')
  return(brokenStickSim)
}
