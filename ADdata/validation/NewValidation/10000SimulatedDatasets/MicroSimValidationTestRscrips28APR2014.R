library(HMP)
library(survcomp)


###Function to compare mean proportions and overdispersion using Xdc.sevsample test
###Prints out vector of pvalues 
meanVarCheckPvalues <- function(simlist, observedData, reads){
  observedData <- round(observedData*reads)
  listLength <- length(simlist)
  pValues <- numeric(length=listLength)
  for(i in 1:listLength){
    pValues[i] <- Xdc.sevsample(list(simlist[[i]], observedData))$'p value'    
  }
  return(pValues)
}

###AD datasets vs Simulated AD Datasets (10,000 sets) Phylum level
pSimObsControlPvalues <- meanVarCheckPvalues(simlist=pControlDatalist, observedData=pcontrol, reads=3000)
mean(pSimObsControlPvalues)
length(which(pSimObsControlPvalues > 0.05))/length(pSimObsControlPvalues)


pSimObsBaselinePvalues <- meanVarCheckPvalues(simlist=pBaselineDatalist, observedData=pbaseline, reads=3000)
mean(pSimObsBaselinePvalues)
length(which(pSimObsBaselinePvalues > 0.05))/length(pSimObsBaselinePvalues)


pSimObsFlarentPvalues <- meanVarCheckPvalues(simlist=pFlarentDatalist, observedData=pflarent, reads=3000)
mean(pSimObsFlarentPvalues)
length(which(pSimObsFlarentPvalues > 0.05))/length(pSimObsFlarentPvalues)


pSimObsFlaretPvalues <- meanVarCheckPvalues(simlist=pFlaretDatalist, observedData=pflaret, reads=3000)
mean(pSimObsFlaretPvalues)
length(which(pSimObsFlaretPvalues > 0.05))/length(pSimObsFlaretPvalues)


pSimObsPostflarePvalues <- meanVarCheckPvalues(simlist=pPostflareDatalist, observedData=ppostflare, reads=3000)
mean(pSimObsPostflarePvalues)
length(which(pSimObsPostflarePvalues > 0.05))/length(pSimObsPostflarePvalues)


###Combined p-values using Stouffer's z-score method (survcomp package)
combine.test(p=pSimObsControlPvalues, method="z.transform")
combine.test(p=pSimObsBaselinePvalues, method="z.transform")

pSimObsFlarentPvaluesEdited <- pSimObsFlarentPvalues
pSimObsFlarentPvaluesEdited[pSimObsFlarentPvaluesEdited == 0] <- 0.000000000000001
combine.test(p=pSimObsFlarentPvaluesEdited, method="z.transform") ###Change to small float value

combine.test(p=pSimObsFlaretPvalues, method="z.transform")
combine.test(p=pSimObsPostflarePvalues, method="z.transform")

hist(pSimObsControlPvalues, 100)
hist(pSimObsBaselinePvalues, 100)
hist(pSimObsFlarentPvalues, 100)
hist(pSimObsFlaretPvalues, 100)
hist(pSimObsPostflarePvalues, 100)


###Saving data objects in .RData format
save(pControlDatalist,file="pControlDatalist.RData")
save(pBaselineDatalist, file="pBaselineDatalist.RData")
save(pFlarentDatalist, file="pFlarentDatalist.RData")
save(pFlaretDatalist, file="pFlaretDatalist.RData")
save(pPostflareDatalist, file="pPostflareDatalist.RData")

save(pcontrol, pbaseline, pflarent, pflaret, ppostflare, file="PhylumObservedDatasetsAD.RData")

save(pSimObsControlPvalues, file="pSimObsControlPvalues.RData")
save(pSimObsBaselinePvalues, file="pSimObsBaselinePvalue.RData")
save(pSimObsFlarentPvalues, file="pSimObsFlarentPvalues.RData")
save(pSimObsFlaretPvalues, file="pSimObsFlaretPvalues.RData")
save(pSimObsPostflarePvalues, file="pSimObsPostflarePvalue.RData")

######################################################################
###AD datasets vs Simulated AD Datasets (10,000 sets) Species level###

sSimObsControlPvalues <- meanVarCheckPvalues(simlist=sControlDatalist, observedData=scontrol, reads=3000)
mean(sSimObsControlPvalues)
length(which(sSimObsControlPvalues > 0.05))/length(sSimObsControlPvalues)

sSimObsBaselinePvalues <- meanVarCheckPvalues(simlist=sBaselineDatalist, observedData=sbaseline, reads=3000)
mean(sSimObsBaselinePvalues)
length(which(sSimObsBaselinePvalues > 0.05))/length(sSimObsBaselinePvalues)

sSimObsFlarentPvalues <- meanVarCheckPvalues(simlist=sFlarentDatalist, observedData=sflarent, reads=3000)
mean(sSimObsFlarentPvalues)
length(which(sSimObsFlarentPvalues > 0.05))/length(sSimObsFlarentPvalues)

sSimObsFlaretPvalues <- meanVarCheckPvalues(simlist=sFlaretDatalist, observedData=sflaret, reads=3000)
mean(sSimObsFlaretPvalues)
length(which(sSimObsFlaretPvalues > 0.05))/length(sSimObsFlaretPvalues)

sSimObsPostflarePvalues <- meanVarCheckPvalues(simlist=sPostflareDatalist, observedData=spostflare, reads=3000)
mean(sSimObsPostflarePvalues)
length(which(sSimObsPostflarePvalues > 0.05))/length(sSimObsPostflarePvalues)

########
###Combined p-values using Stouffer's z-score method (survcomp package) ### change all 0 to small float
###Change zeros to very small float (0.000000000000001) for the ones that produce NaN
sSimObsControlPvaluesEdited <- sSimObsControlPvalues
sSimObsControlPvaluesEdited[sSimObsControlPvaluesEdited == 0] <- 0.000000000000001
combine.test(p=sSimObsControlPvaluesEdited, method="z.transform")

sSimObsBaselinePvaluesEdited <- sSimObsBaselinePvalues
sSimObsBaselinePvaluesEdited[sSimObsBaselinePvaluesEdited == 0] <- 0.000000000000001
combine.test(p=sSimObsBaselinePvaluesEdited, method="z.transform")

sSimObsFlarentPvaluesEdited <- sSimObsFlarentPvalues
sSimObsFlarentPvaluesEdited[sSimObsFlarentPvaluesEdited == 0] <- 0.000000000000001
combine.test(p=sSimObsFlarentPvaluesEdited, method="z.transform") 

sSimObsFlaretPvaluesEdited <- sSimObsFlaretPvalues
sSimObsFlaretPvaluesEdited[sSimObsFlaretPvaluesEdited == 0] <- 0.000000000000001
combine.test(p=sSimObsFlaretPvaluesEdited, method="z.transform")

sSimObsPostflarePvaluesEdited <- sSimObsPostflarePvalues
sSimObsPostflarePvaluesEdited[sSimObsPostflarePvaluesEdited == 0] <- 0.000000000000001
combine.test(p=sSimObsPostflarePvaluesEdited, method="z.transform")

hist(sSimObsControlPvalues, 100)
hist(sSimObsBaselinePvalues, 100)
hist(sSimObsFlarentPvalues, 100)
hist(sSimObsFlaretPvalues, 100)
hist(sSimObsPostflarePvalues, 100)


save(sControlDatalist, file="sControlDatalist.RData")
save(sBaselineDatalist, file="sBaselineDatalist.RData")
save(sFlarentDatalist, file="sFlarentDatalist.RData")
save(sFlaretDatalist, file="sFlaretDatalist.RData")
save(sPostflareDatalist, file="sPostflareDatalist.RData")

save(scontrol, sbaseline, sflarent, sflaret, spostflare, file="SpeciesObservedDatasetsAD.RData")

save(sSimObsControlPvalues, file="sSimObsControlPvalues.RData")
save(sSimObsBaselinePvalues, file="sSimObsBaselinePvalue.RData")
save(sSimObsFlarentPvalues, file="sSimObsFlarentPvalues.RData")
save(sSimObsFlaretPvalues, file="sSimObsFlaretPvalues.RData")
save(sSimObsPostflarePvalues, file="sSimObsPostflarePvalue.RData")

########################################################
###RAT DATA#####
################
##Note to self: use ...LSCsorted instead of pLS... data because the ...LSCsorted is ranked and pLS is not.
##To run our simulations, the data must be ranked first, ...LSCsorted is the ranked dataset
###of the published data

###Runs the hypothesis testing
d2LSCpValues <- meanVarCheckPvalues(simlist=vLSCCED2, observedData=d2LSCsorted, reads=3000)
mean(d2LSCpValues)
length(which(d2LSCpValues > 0.05))/length(d2LSCpValues)

d4LSCCpValues <- meanVarCheckPvalues(simlist=vLSCCD4, observedData=d4LSCCsorted, reads=3000)
mean(d4LSCCpValues)
length(which(d4LSCCpValues > 0.05))/length(d4LSCCpValues)

d4LSCEpValues <- meanVarCheckPvalues(simlist=vLSCED4, observedData=d4LSCEsorted, reads=3000)
mean(d4LSCEpValues)
length(which(d4LSCEpValues > 0.05))/length(d4LSCCpValues)

###Changes the zeros to small floats to run the z-score method of combining p values
d2LSCpValuesEdited <- d2LSCpValues
d2LSCpValuesEdited[d2LSCpValuesEdited == 0] <- 0.000000000000001
combine.test(p=d2LSCpValuesEdited, method="z.transform")
hist(d2LSCpValuesEdited, 100)

d4LSCCpValuesEdited <- d4LSCCpValues
d4LSCCpValuesEdited[d4LSCCpValuesEdited == 0] <- 0.000000000000001
combine.test(p=d4LSCCpValuesEdited, method="z.transform")
hist(d4LSCCpValuesEdited, 100)

combine.test(p=d4LSCEpValues, method="z.transform")
hist(d4LSCEpValues, 100)

save(LSCCD4, file="LSCCD4Pub.RData")
save(LSCED4, file="LSCED4Pub.RData")
LSCCED2 <- LsericataD2
save(LSCCED2, file="LSCCED2Pub.RData")

save(pLSCCD4, file="pLSCCD4.RData")
save(pLSCED4, file="pLSCED4.RData")
pLSCCED2 <- pLsericataD2
save(pLSCCED2, file="pLSCCED2.RData")

save(d2LSCsorted, file="d2LSCsorted.RData")
save(d4LSCCsorted, file="d4LSCCsorted.RData")
save(d4LSCEsorted, file="d4LSCEsorted.RData")

save(vLSCCD4, file="vLSCCD4.RData")
save(vLSCCED2, file="vLSCCED2.RData")
save(vLSCED4, file="vLSCED4.RData")


