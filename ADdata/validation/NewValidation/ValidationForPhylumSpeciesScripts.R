###For my computer
setwd("~/data sim/ADdata/validation/NewValidation")
###

###Load data
setwd("~/data sim/ADdata/validation/Phylum/RData files")
###

save(vcontrol,file="ADControlPhylumVdataset1000.RData")
save(vbaseline, file="ADBaselinePhylumVdataset1000.RData")
save(vflarent, file="ADFlareNTPhylumVdataset1000.RData")
save(vflaret, file="ADFlareTPhylumVdataset1000.RData")
save(vpostflare, file="ADPostFlarePhylumVdataset1000.RData")

save(vdscontrols,file="ADControlSpeciesVdataset1000.RData")
save(vdsbaseline, file="ADBaselineSpeciesVdataset1000.RData")
save(vdsflarent, file="ADFlareNTSpeciesVdataset1000.RData")
save(vdsflaret, file="ADFlareTSpeciesVdataset1000.RData")
save(vdspostflare, file="ADPostFlareSpeciesVdataset1000.RData")

############
pcontrol <- read.table("ADControls.txt")
(pcontrol <- round(pcontrol*3000))
(pcontrol <- t(pcontrol))
pbaseline <- read.table("ADBaseline.txt")
(pbaseline <- round(pbaseline*3000))
(pbaseline <- t(pbaseline))
pflarent <- read.table("ADFlareNT.txt")
(pflarent <- round(pflarent*3000))
(pflarent <- t(pflarent))
pflaret <- read.table("ADFlareT.txt")
(pflaret <- round(pflaret*3000))
(pflaret <- t(pflaret))
ppostflare <- read.table("ADPostFlare.txt")
(ppostflare <- round(ppostflare*3000))
(ppostflare <- t(ppostflare))

scontrol <- read.table("ADControl.txt")
(scontrol <- round(scontrol*3000))
sbaseline <- read.table("ADBaseline.txt")
(sbaseline <- round(sbaseline*3000))
sflarent <- read.table("ADFlareNT.txt")
sflarent <- round(sflarent*3000)
sflaret <- read.table("ADFlareT.txt")
(sflaret <- round(sflaret*3000))
spostflare <- read.table("ADPostflare.txt")
(spostflare <- round(spostflare*3000))

##############################################
##Phylum#Probability-Mean#Hypothesis#Testing##
##############################################

#Null Hypothesis: the simulated and provided datasets are the same
#Alternative Hypothesis: the simulated and the provided datasets are different

#Hypothesis testing using Generalized Wald-type Statistics: Several Sample RAD Probability-Mean
##Test comparision with a Known Common Vector (Xmc.sevsample)

####################
##Phylum#ADControl##
####################

fit.pcontrol <- DM.MoM(pcontrol)
pcontrolresults <- Xmc.sevsample(vcontrol, fit.pcontrol$pi)

#####################
##Phylum#ADBaseline##
#####################

fit.pbaseline <- DM.MoM(pbaseline)
pbaselineresults <- Xmc.sevsample(vbaseline, fit.pbaseline$pi)

####################
##Phylum#ADFlareNT##
####################

fit.pflarent <- DM.MoM(pflarent)
pflarentresults <- Xmc.sevsample(vflarent, fit.pflarent$pi)

###################
##Phylum#ADFlareT##
###################

fit.pflaret <- DM.MoM(pflaret)
pflaretresult <- Xmc.sevsample(vflaret, fit.pflaret$pi)

######################
##Phylum#ADPostflare##
######################

fit.ppostflare <- DM.MoM(ppostflare)
ppostflareresults <- Xmc.sevsample(vpostflare, fit.ppostflare$pi)

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

thetaTest2tailed <- function(simDatalist, observedData){
  thetalist <- getthetaValues(simDatalist)
  observed_theta <- DM.MoM(observedData)$theta
  pvalue <- (length(which(thetalist > observed_theta))/length(thetalist))/2
  return(pvalue)
}

####################
##Phylum#ADControl##
####################

(pcontroltheta <- DM.MoM(pcontrol)$theta)
pcontrolthetalist <- getthetaValues(vcontrol)
hist(pcontrolthetalist)
abline(v=0.1383081)

thetaTest2tailed(vcontrol, pcontrol)

#####################
##Phylum#ADBaseline##
#####################

(pbaselinetheta <- DM.MoM(pbaseline)$theta)
pbaselinethetalist <- getthetaValues(vbaseline)
hist(pbaselinethetalist)

thetaTest2tailed(vbaseline, pbaseline)

####################
##Phylum#ADFlareNT##
####################

(pflarenttheta <- DM.MoM(pflarent)$theta)
pflarentthetalist <- getthetaValues(vflarent)
hist(pflarentthetalist)

thetaTest2tailed(vflarent, pflarent)

###################
##Phylum#ADFlareT##
###################

(pflarettheta <- DM.MoM(pflaret)$theta)
pflaretthetalist <- getthetaValues(vflaret)
hist(pflaretthetalist)

thetaTest2tailed(vflaret, pflaret)

######################
##Phylum#ADPostflare##
######################

(ppostflaretheta <- DM.MoM(ppostflare)$theta)
ppostflarethetalist <- getthetaValues(vpostflare)
hist(ppostflarethetalist)

thetaTest2tailed(vpostflare, ppostflare)

##############################################
##Species#Probability-Mean#Hypothesis#Testing##
##############################################

#Null Hypothesis: the simulated and provided datasets are the same
#Alternative Hypothesis: the simulated and the provided datasets are different

#Hypothesis testing using Generalized Wald-type Statistics: Several Sample RAD Probability-Mean
##Test comparision with a Known Common Vector (Xmc.sevsample)

#####################
##Species#ADControl##
#####################

fit.scontrol <- DM.MoM(scontrol)
(scontrolresults <- Xmc.sevsample(vdscontrols, fit.scontrol$pi))

######################
##Species#ADBaseline##
######################

fit.sbaseline <- DM.MoM(sbaseline)
(sbaselineresults <- Xmc.sevsample(vdsbaseline, fit.sbaseline$pi))

#####################
##Species#ADFlareNT##
#####################

fit.sflarent <- DM.MoM(sflarent)
(sflarentresults <- Xmc.sevsample(vdsflarent, fit.sflarent$pi))

####################
##Species#ADFlareT##
####################

fit.sflaret <- DM.MoM(sflaret)
(sflaretresults <- Xmc.sevsample(vdsflaret, fit.sflaret$pi))

#######################
##Species#ADPostflare##
#######################

fit.spostflare <- DM.MoM(spostflare)
(spostflareresults <- Xmc.sevsample(vdspostflare, fit.spostflare$pi))

######################################################
##Species#Overdispersion-Variance#Hypothesis#Testing##
######################################################

#####################
##Species#ADControl##
#####################

(scontroltheta <- DM.MoM(scontrol)$theta)
scontrolthetalist <- getthetaValues(vdscontrols)
hist(scontrolthetalist)

thetaTest2tailed(vdscontrols, scontrol)

######################
##Species#ADBaseline##
######################

(sbaselinetheta <- DM.MoM(sbaseline)$theta)
sbaselinethetalist <- getthetaValues(vdsbaseline)
hist(sbaselinethetalist)
abline(v=sbaselinetheta)

thetaTest2tailed(vdsbaseline, sbaseline)

#####################
##Species#ADFlareNT##
#####################

(sflarenttheta <- DM.MoM(sflarent)$theta)
sflarentthetalist <- getthetaValues(vdsflarent)
hist(sbaselinethetalist)

thetaTest2tailed(vdsflarent, sflarent)

####################
##Species#ADFlareT##
####################

(sflarettheta <- DM.MoM(sflaret)$theta)
sflaretthetalist <- getthetaValues(vdsflaret)
hist(sflaretthetalist)

thetaTest2tailed(vdsflaret, sflaret)

#######################
##Species#ADPostflare##
#######################

(spostflaretheta <- DM.MoM(spostflare)$theta)
spostflarethetalist <- getthetaValues(vdspostflare)
hist(spostflarethetalist)
plot(density(spostflarethetalist))
abline(v=0.1388843)

thetaTest2tailed(vdspostflare, spostflare)
