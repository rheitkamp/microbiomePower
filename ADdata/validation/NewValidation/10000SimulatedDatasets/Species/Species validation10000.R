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
###Loading data Species###
#########################

###Make sure col=taxon and row=samples
###Data should be ranked

scontrol <- read.table("ADControl.txt")
sbaseline <- read.table("ADBaseline.txt")
sflarent <- read.table("ADFlareNT.txt")
sflaret <- read.table("ADFlareT.txt")
spostflare <- read.table("ADPostFlare.txt")

sControlDatalist <- validationdatasets(rawData=scontrol,numberSubjects=22, numberSimulations=10000)
sBaselineDatalist <- validationdatasets(rawData=sbaseline, numberSubjects=12, numberSimulations=10000)
sFlarentDatalist <- validationdatasets(rawData=sflarent, numberSubjects=7, numberSimulations=10000)
sFlaretDatalist <- validationdatasets(rawData=sflaret, numberSubjects=5, numberSimulations=10000)
sPostflareDatalist <- validationdatasets(rawData=spostflare, numberSubjects=11, numberSimulations=10000)

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
###Species Probability-Mean Hypothesis Testing###
################################################

#Null Hypothesis: the simulated and provided datasets are the same
#Alternative Hypothesis: the simulated and the provided datasets are different
#Significance Value of 0.05

#Hypothesis testing using Generalized Wald-type Statistics: Several Sample RAD Probability-Mean
##Test comparision with a Known Common Vector (Xmc.sevsample) 
###from La Rosa HMP Package

######################
###Species ADControl###
######################

(mscontrolResults <- meanCheck(simDatalist=sControlDatalist, observedData=scontrol, numReads=3000))

#######################
###Species ADBaseline###
#######################

(msbaselineResults <- meanCheck(simDatalist=sBaselineDatalist, observedData=sbaseline, numReads=3000))

######################
###Species ADFlareNT###
######################

(msflarentResults <- meanCheck(simDatalist=sFlarentDatalist, observedData=sflarent, numReads=3000))

#####################
###Species ADFlareT###
#####################

(msflaretResults <- meanCheck(simDatalist=sFlaretDatalist, observedData=sflaret, numReads=3000))

########################
###Species ADPostflare###
########################

(mspostflareResults <- meanCheck(simDatalist=sPostflareDatalist, observedData=spostflare, numReads=3000))

#####################################################
##Species#Overdispersion-Variance#Hypothesis#Testing##
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
###Species ADControl###
######################

scontrolsimthetaList <- getthetaValues(sControlDatalist)
hist(scontrolsimthetaList, 25)
scontroltheta <- DM.MoM(scontrol*3000)$theta
abline(v=scontroltheta)

qqPlot(scontrolsimthetaList)

(tscontrolResults <- thetaTest2tailed(simDatalist=sControlDatalist, observedData=scontrol, numReads=3000))

###Transformation ADControl###

logScontrol <- log10(scontrolsimthetaList)
qqPlot(logScontrol)

sqrtScontrol <- sqrt(scontrolsimthetaList)
qqPlot(sqrtScontrol)

asinScontrol <- asin(scontrolsimthetaList)
qqPlot(asinScontrol)
###
(length(which(sqrtScontrol > sqrt(scontroltheta)))/length(sqrtScontrol))/2

#######################
###Species ADBaseline###
#######################

sbaselinesimthetaList <- getthetaValues(sBaselineDatalist)
hist(sbaselinesimthetaList, 25)
sbaselinetheta <- DM.MoM(sbaseline*3000)$theta
abline(v=sbaselinetheta)

qqPlot(sbaselinesimthetaList)

(tsbaselineResults <- thetaTest2tailed(simDatalist=sBaselineDatalist, observedData=sbaseline, numReads=3000))

###ADBaseline Transformation###

logSbaseline <- log10(sbaselinesimthetaList)
qqPlot(logSbaseline)

sqrtSbaseline <- sqrt(sbaselinesimthetaList)
qqPlot(sqrtSbaseline)

asinSbaseline <- asin(sbaselinesimthetaList)
qqPlot(asinSbaseline)
###
(length(which(sqrtSbaseline > sqrt(sbaselinetheta)))/length(sqrtSbaseline))/2

######################
###Phlyum ADFlareNT###
######################

sflarentsimthetaList <- getthetaValues(sFlarentDatalist)
hist(sflarentsimthetaList, 25)
sflarenttheta <- DM.MoM(sflarent*3000)$theta
abline(v=sflarenttheta)

qqPlot(sflarentsimthetaList)

(tsflarentResults <- thetaTest2tailed(simDatalist=sFlarentDatalist, observedData=sflarent, numReads=3000))

###ADFlareNT Transform###
logSflarent <- log10(sflarentsimthetaList)
qqPlot(logSflarent)

sqrtSflarent <- sqrt(sflarentsimthetaList)
qqPlot(sqrtSflarent)

asinSflarent <- asin(sflarentsimthetaList)
qqPlot(asinSflarent)
###
(length(which(sqrtSflarent > sqrt(sflarenttheta)))/length(sqrtSflarent))/2

#####################
###Phlyum ADFlareT###
#####################
sflaretsimthetaList <- getthetaValues(sFlaretDatalist)
sflarettheta <- DM.MoM(sflaret*3000)$theta

hist(sflaretsimthetaList, 25)
hist(sqrt(sflaretsimthetaList), 25)
qqPlot(sflaretsimthetaList)
###
(tsflaretResults <- thetaTest2tailed(simDatalist=sFlaretDatalist, observedData=sflaret, numReads=3000))

###ADFlaret Transform###
logSflaret <- log10(sflaretsimthetaList)
qqPlot(logSflaret)

sqrtSflaret <- sqrt(sflaretsimthetaList)
qqPlot(sqrtSflaret)

asinSflaret <- asin(sflaretsimthetaList)
qqPlot(asinSflaret)
###
(length(which(sqrtSflaret > sqrt(sflarettheta)))/length(sqrtSflaret))/2

########################
###Species ADPostflare###
########################
spostflaresimthetaList <- getthetaValues(sPostflareDatalist)
spostflaretheta <- DM.MoM(spostflare*3000)$theta
qqPlot(spostflaresimthetaList)

(tspostflareResults <- thetaTest2tailed(simDatalist=sPostflareDatalist, observedData=spostflare, numReads=3000))

###ADPostflare Transform###
logSpostflare <- log10(spostflaresimthetaList)
qqPlot(logSpostflare)

sqrtSpostflare <- sqrt(spostflaresimthetaList)
qqPlot(sqrtSpostflare)

asinSpostflare <- asin(spostflaresimthetaList)
qqPlot(asinSpostflare)
####
(length(which(sqrtSpostflare > sqrt(spostflaretheta)))/length(sqrtSpostflare))/2