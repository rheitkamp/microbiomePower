library(HMP)
library(PearsonDS)
set.seed(1234)

#####################################
###Function for dataset sumulation###
#####################################

validationdatasets <- function(rawData, numberSubjects, numberSimulations){
  listdataset <- replicate(n=length(numberSimulations), expr=list())
  for(i in 1:numberSimulations){
    listdataset[[i]] <- round((simulateBrokenStick(rawData, numberSubjects))*3000)
  }
  return(listdataset)
}

#########################
###Loading data Phylum###
#########################

###Make sure col=taxon and row=samples###

pcontrol <- t(read.table("ADControls.txt"))
pbaseline <- t(read.table("ADBaseline.txt"))
pflarent <- t(read.table("ADFlareNT.txt"))
pflaret <- t(read.table("ADFlareT.txt"))
ppostflare <- t(read.table("ADPostFlare.txt"))

pControlDatalist <- validationdatasets(rawData=pcontrol,numberSubjects=22, numberSimulations=10000)
pBaselineDatalist <- validationdatasets(rawData=pbaseline, numberSubjects=11, numberSimulations=10000)
pFlarentDatalist <- validationdatasets(rawData=pflarent, numberSubjects=7, numberSimulations=10000)
pFlaretDatalist <- validationdatasets(rawData=pflaret, numberSubjects=5, numberSimulations=10000)
pPostflareDatalist <- validationdatasets(rawData=ppostflare, numberSubjects=11, numberSimulations=10000)

###############################
###Function for Mean Testing###
###############################

meanCheck <- function(simDatalist, observedData, numReads){
  observedData <- round(observedData*numReads)
  fit <- DM.MoM(observedData)
  results <- Xmc.sevsample(group.data=simDatalist, pi0=fit$pi)
  return(results)
}

################################################
###Phylum Probability-Mean Hypothesis Testing###
################################################

#Null Hypothesis: the simulated and provided datasets are the same
#Alternative Hypothesis: the simulated and the provided datasets are different
#Significance Value of 0.05

#Hypothesis testing using Generalized Wald-type Statistics: Several Sample RAD Probability-Mean
##Test comparision with a Known Common Vector (Xmc.sevsample) 
###from La Rosa HMP Package

######################
###Phylum ADControl###
######################

(mpcontrolResults <- meanCheck(simDatalist=pControlDatalist, observedData=pcontrol, numReads=3000))

#######################
###Phylum ADBaseline###
#######################

(mpbaselineResults <- meanCheck(simDatalist=pBaselineDatalist, observedData=pbaseline, numReads=3000))

######################
###Phylum ADFlareNT###
######################

(mpflarentResults <- meanCheck(simDatalist=pFlarentDatalist, observedData=pflarent, numReads=3000))

#####################
###Phylum ADFlareT###
#####################

(mpflaretResults <- meanCheck(simDatalist=pFlaretDatalist, observedData=pflaret, numReads=3000))

########################
###Phylum ADPostflare###
########################

(mppostflareResults <- meanCheck(simDatalist=pPostflareDatalist, observedData=ppostflare, numReads=3000))

#####################################################
##Phylum#Overdispersion-Variance#Hypothesis#Testing##
#####################################################

########################
###Functions for Test###
########################
getthetaValues <- function(x){
  numlist <- length(x)
  thetaValues <- numeric(length=numlist)
  for(i in 1:numlist){
    thetaValues[i] <- DM.MoM(x[[i]])$theta
  }
  return(thetaValues)
}

thetaTest2tailed <- function(simDatalist, observedData, numReads){
  thetalist <- getthetaValues(simDatalist)
  observedData <- round(observedData*numReads)
  observed_theta <- DM.MoM(observedData)$theta
  pvalue <- (length(which(thetalist > observed_theta))/length(thetalist))/2
  pvalue <- list(pvalue)
  names(pvalue) <- c("p value")
  return(pvalue)
}

######################
###Phylum ADControl###
######################
pcontrolsimthetaList <- getthetaValues(pControlDatalist)
pcontroltheta <- DM.MoM(pcontrol*3000)$theta
qqPlot(pcontrolsimthetaList)

(tpcontrolResults <- thetaTest2tailed(simDatalist=pControlDatalist, observedData=pcontrol, numReads=3000))

###Transform
logPcontrol <- log10(pcontrolsimthetaList)
qqPlot(logPcontrol)

sqrtPcontrol <- sqrt(pcontrolsimthetaList)
qqPlot(sqrtPcontrol)

asinPcontrol <- asin(pcontrolsimthetaList)
qqPlot(asinPcontrol)
###
(length(which(sqrtPcontrol > sqrt(pcontroltheta)))/length(sqrtPcontrol))/2

#######################
###Phylum ADBaseline###
#######################
pbaselinesimthetaList <- getthetaValues(pBaselineDatalist)
pbaselinetheta <- DM.MoM(pbaseline*3000)$theta
qqPlot(pbaselinesimthetaList)

(tpbaselineResults <- thetaTest2tailed(simDatalist=pBaselineDatalist, observedData=pbaseline, numReads=3000))

###Transform
logPbaseline <- log10(pbaselinesimthetaList)
qqPlot(logPbaseline)

sqrtPbaseline <- sqrt(pbaselinesimthetaList)
qqPlot(sqrtPbaseline)

asinPbaseline <- asin(pbaselinesimthetaList)
qqPlot(asinPbaseline)
###
(length(which(sqrtPbaseline > sqrt(pbaselinetheta)))/length(sqrtPbaseline))/2

######################
###Phlyum ADFlareNT###
######################
pflarentsimthetaList <- getthetaValues(pFlarentDatalist)
pflarenttheta <- DM.MoM(pflarent*3000)$theta
qqPlot(pflarentsimthetaList)

(tpflarentResults <- thetaTest2tailed(simDatalist=pFlarentDatalist, observedData=pflarent, numReads=3000))

###Transform
logPflarent <- log10(pflarentsimthetaList)
qqPlot(logPflarent)

sqrtPflarent <- sqrt(pflarentsimthetaList)
qqPlot(sqrtPflarent)

hist(ppostflaresimthetaList, 25)
hist(sqrtPflarent, 25)

qqPlot(exp(pflarentsimthetaList))

asinPflarent <- asin(pflarentsimthetaList)
qqPlot(asinPflarent)
###
(length(which(sqrtPflarent > sqrt(pflarenttheta)))/length(sqrtPflarent))/2

#####################
###Phlyum ADFlareT###
#####################
pflaretsimthetaList <- getthetaValues(pFlaretDatalist)
pflarettheta <- DM.MoM(pflaret*3000)$theta
qqPlot(pflaretsimthetaList)

(tpflaretResults <- thetaTest2tailed(simDatalist=pFlaretDatalist, observedData=pflaret, numReads=3000))

###Transform
logPflaret <- log10(pflaretsimthetaList)
qqPlot(logPflaret)

sqrtPflaret <- sqrt(pflaretsimthetaList)
qqPlot(sqrtPflaret)

asinPflaret <- asin(pflaretsimthetaList)
qqPlot(asinPflaret)
###
(length(which(sqrtPflaret > sqrt(pflarettheta)))/length(sqrtPflaret))/2

########################
###Phylum ADPostflare###
########################
ppostflaresimthetaList <- getthetaValues(pPostflareDatalist)
ppostflaretheta <- DM.MoM(ppostflare*3000)$theta
qqPlot(ppostflaresimthetaList)

(tppostflareResults <- thetaTest2tailed(simDatalist=pPostflareDatalist, observedData=ppostflare, numReads=3000))

###Transform
logPpostflare <- log10(ppostflaresimthetaList)
qqPlot(logPpostflare)

sqrtPpostflare <- sqrt(ppostflaresimthetaList)
qqPlot(sqrtPpostflare)

asinPpostflare <- asin(ppostflaresimthetaList)
qqPlot(asinPpostflare)
###
(length(which(sqrtPpostflare > sqrt(ppostflaretheta)))/length(sqrtPpostflare))/2
