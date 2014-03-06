library(HMP)


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
saveRDS(object=pControlDatalist, file="pControlDatalist.RData")
saveRDS(object=pBaselineDatalist, file="pBaselineDatalist.RData")
saveRDS(object=pFlarentDatalist, file="pFlarentDatalist.RData")
saveRDS(object=pFlaretDatalist, file="pFlaretDatalist.RData")
saveRDS(object=pPostflareDatalist, file="pPostflareDatalist.RData")

save(pcontrol, pbaseline, pflarent, pflaret, ppostflare, file="PhylumObservedDatasetsAD.RData")

saveRDS(object=pSimObsControlPvalues, file="pSimObsControlPvalues.RData")
saveRDS(object=pSimObsBaselinePvalues, file="pSimObsBaselinePvalue.RData")
saveRDS(object=pSimObsFlarentPvalues, file="pSimObsFlarentPvalues.RData")
saveRDS(object=pSimObsFlaretPvalues, file="pSimObsFlaretPvalues.RData")
saveRDS(object=pSimObsPostflarePvalues, file="pSimObsPostflarePvalue.RData")

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


saveRDS(object=sControlDatalist, file="sControlDatalist.RData")
saveRDS(object=sBaselineDatalist, file="sBaselineDatalist.RData")
saveRDS(object=sFlarentDatalist, file="sFlarentDatalist.RData")
saveRDS(object=sFlaretDatalist, file="sFlaretDatalist.RData")
saveRDS(object=sPostflareDatalist, file="sPostflareDatalist.RData")

save(scontrol, sbaseline, sflarent, sflaret, spostflare, file="SpeciesObservedDatasetsAD.RData")

saveRDS(object=sSimObsControlPvalues, file="sSimObsControlPvalues.RData")
saveRDS(object=sSimObsBaselinePvalues, file="sSimObsBaselinePvalue.RData")
saveRDS(object=sSimObsFlarentPvalues, file="sSimObsFlarentPvalues.RData")
saveRDS(object=sSimObsFlaretPvalues, file="sSimObsFlaretPvalues.RData")
saveRDS(object=sSimObsPostflarePvalues, file="sSimObsPostflarePvalue.RData")

########################################################
