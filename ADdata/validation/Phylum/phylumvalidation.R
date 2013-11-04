####Simulate 10,000 data points
simcontrol <- simulateBrokenStick("ADControls.txt", "Control", 10000)
simbaseline <- simulateBrokenStick("ADBaseline.txt", "ADBaseline", 10000)
simflarent <- simulateBrokenStick("ADFlareNT.txt", "ADFlareNT", 10000)
simflaret <- simulateBrokenStick("ADFlareT.txt", "ADFlareT", 10000)
simpostflare <- simulateBrokenStick("ADPostFlare.txt", "Postflare", 10000)

####Save
write.table(x=simcontrol, file="simcontrol10000samples.txt")
write.table(x=simbaseline, file="simbaseline10000samples.txt")
write.table(x=simflarent, file="simflarent10000samples.txt")
write.table(x=simflaret, file="simflaret10000samples.txt")
write.table(x=simpostflare, file="simpostflare10000samples.txt")
####
picheckcontrol <- picheck("ADControls.txt", simcontrol, 3000)
picheckbaseline <- picheck("ADBaseline.txt", simbaseline, 3000)
picheckflarent <- picheck("ADFlareNT.txt", simflarent, 3000)
picheckflaret <- picheck("ADFlareT.txt", simflaret, 3000)
picheckpostflare <- picheck("ADPostFlare.txt", simpostflare, 3000)

####Save
write.table(x=picheckcontrol, file="picheckcontrol.txt")
write.table(x=picheckbaseline, file="picheckbaseline.txt")
write.table(x=picheckflarent, file="picheckflarent.txt")
write.table(x=picheckflaret, file="picheckflaret.txt")
write.table(x=picheckpostflare, file="picheckpostflare.txt")

####
ocheckcontrol <- Ocheck("ADControls.txt", simcontrol, 3000)
ocheckbaseline <- Ocheck("ADBaseline.txt", simbaseline, 3000)
ocheckflarent <- Ocheck("ADFlareNT.txt", simflarent, 3000)
ocheckflaret <- Ocheck("ADFlareT.txt", simflaret, 3000)
ocheckpostflare <- Ocheck("ADPostFlare.txt", simpostflare, 3000)

####Save
write.table(x=ocheckcontrol, file="varscheckcontrol.txt")
write.table(x=ocheckbaseline, file="varscheckbaseline.txt")
write.table(x=ocheckflarent, file="varscheckflarent.txt")
write.table(x=ocheckflaret, file="varscheckflaret.txt")
write.table(x=ocheckpostflare, file="varscheckpostflare.txt")
####


####Data check
picheck <- function(inputFilename, simdata, numReads){
  simdatareads <- round(simdata*numReads)
  pubdata <- round(t(read.table(inputFilename))*numReads)
  fit.simdata <- dirmult(simdatareads)
  Check <- Xsc.onesample(pubdata,fit.simdata$pi)
  print(head(simdatareads))
  print(head(pubdata))
  return(Check)
}
####

####
Ocheck <- function(inputFilename, simdata, numReads){
  simdatareads <- round(simdata*numReads)
  pubdata <- round(t(read.table(inputFilename))*numReads)
  listdata <- list(simdatareads,pubdata)
  oCheck <- Xoc.sevsample(listdata)
  print(head(simdatareads))  
  print(head(pubdata))
  return(oCheck)
}
####
####1. Simulate dataset for all
vcontrol <- validationdatasets(simdata=simcontrol,numberSimulations=1000,publishNumbersubjects=22,numReads=3000)
write.table(x=vcontrol,file="vcontroldatasets1000.txt")

vbaseline <- validationdatasets(simdata=simbaseline,numberSimulations=1000,publishNumbersubjects=11,numReads=3000)
write.table(x=vbaseline, file="vbaselinedatasets1000.txt")

vflarent <- validationdatasets(simdata=simflarent, numberSimulations=1000, publishNumbersubjects= 7, numReads=3000) 
write.table(x=vflarent, file="vflarentdatasets1000.txt")

vflaret <- validationdatasets(simdata=simflaret, numberSimulations=1000, publishNumbersubjects=5, numReads=3000) 
write.table(x=vflaret, file="vflaretdatasets1000.txt")

vpostflare <- validationdatasets(simdata=simpostflare, numberSimulations=1000, publishNumbersubjects=11, numReads=3000)
write.table(x=vpostflare, file="vpostflaredatasets1000.txt")

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


meancheckcontrol <- meancheck("ADControls.txt", vcontrol, 3000)
mean(meancheckcontrol < 0.05)
write.table(meancheckcontrol, "meancheckcontrolpvalues1000.txt")

meancheckbaseline <- meancheck("ADBaseline.txt", vbaseline, 3000)
mean(meancheckbaseline < 0.05)
write.table(meancheckbaseline, "meancheckbaselinepvalues1000.txt")

meancheckflarent <- meancheck("ADFlareNT.txt", vflarent, 3000)
mean(meancheckflarent < 0.05)
write.table(meancheckflarent, "meancheckflarentpvalues1000.txt")

meancheckflaret <- meancheck("ADFlareT.txt", vflaret, 3000)
mean(meancheckflaret < 0.05)
write.table(meancheckflaret, "meancheckflaretpvalues1000.txt")

meancheckpostflare <- meancheck("ADPostFlare.txt", vpostflare, 3000)
mean(meancheckpostflare < 0.05)
write.table(meancheckpostflare, "meancheckpostflarepvalues1000.txt")

###Varcheck  
varcheckcontrol <- variancecheck("ADControls.txt",vcontrol,3000)
mean(varcheckcontrol < 0.05)
mean(varcheckcontrol)
median(varcheckcontrol)
write.table(varcheckcontrol, "varcheckcontrolpvalues1000.txt")

varcheckbaseline <- variancecheck("ADBaseline.txt",vbaseline,3000)
mean(varcheckbaseline2 < 0.05)
mean(varcheckbaseline)
median(varcheckbaseline)
write.table(varcheckbaseline, "varcheckbaselinepvalues1000.txt")

varcheckflarent <- variancecheck("ADFlareNT.txt",vflarent,3000)
mean(varcheckflarent < 0.05)
mean(varcheckflarent)
median(varcheckflarent)
write.table(varcheckflarent, "varcheckflarentpvalues1000.txt")

varcheckflaret <- variancecheck("ADFlareT.txt", vflaret, 3000)
mean(varcheckflaret < 0.05)
mean(varcheckflaret)
median(varcheckflaret)
write.table(varcheckflaret, "varcheckflaretpvalues1000.txt")

varcheckpostflare <- variancecheck("ADPostFlare.txt", vpostflare, 3000)
mean(varcheckpostflare < 0.05)
mean(varcheckpostflare)
median(varcheckpostflare)
write.table(varcheckpostflare, "varcheckpostflare.txt")

####
#24Oct2013
#####Validation datasets without drawing from a pre created dataset
validationdatasets2 <- function(inputFilename,outputLabel,numberSubjects,numberSimulations){
  listdataset <- replicate(n=length(numberSimulations), expr=list())
  for(i in 1:numberSimulations){
    listdataset[[i]] <- round((simulateBrokenStick(inputFilename,outputLabel,numberSubjects))*3000)
  }
  return(listdataset)
}

control <- read.table("ADControls.txt")
baseline <- read.table("ADBaseline.txt")
flarent <- read.table("ADFlareNT.txt")
flaret <- read.table("ADFlareT.txt")
postflare <- read.table("ADBaseline.txt")
dim(read.table("ADControls.txt"))[2]

v2control <- validationdatasets2(inputFilename="ADControls.txt",outputLabel="Control",numberSubjects=22,numberSimulations=1000)
v2baseline <- validationdatasets2(inputFilename="ADBaseline.txt",outputLabel="Baseline",numberSubjects=11,numberSimulations=1000)
v2flarent <- validationdatasets2(inputFilename="ADFlareNT.txt",outputLabel="ADFlareNT",numberSubjects=7,numberSimulations=1000)
v2flaret <- validationdatasets2(inputFilename="ADFlareT.txt",outputLabel="ADFlareT",numberSubjects=5,numberSimulations=1000)
v2postflare <- validationdatasets2(inputFilename="ADPostFlare.txt",outputLabel="ADPostflare",numberSubjects=11,numberSimulations=1000)

write.table(x=v2control,file="v2controldatasets1000.csv", sep=",")
write.table(x=v2baseline,file="v2baselinedatasets1000.txt")
write.table(x=v2flarent, file="v2flarentdatasets1000.csv", sep=",")
write.table(x=v2flaret, file="v2flaretdatasets1000.csv", sep=",")
write.table(x=v2postflare, file="v2postflaredatasets1000.csv", sep=",")

#####
meancheckcontrol2 <- meancheck("ADControls.txt", v2control, 3000)
mean(meancheckcontrol2 < 0.05)
write.table(meancheckcontrol2, "meancheckcontrol2pvalues1000.csv", sep=",")

meancheckbaseline2 <- meancheck("ADBaseline.txt", v2baseline, 3000)
mean(meancheckbaseline2 < 0.05)
write.table(meancheckbaseline2, "meancheckbaseline2pvalues1000.csv", sep=",")

meancheckflarent2 <- meancheck("ADFlareNT.txt", v2flarent, 3000)
mean(meancheckflarent2 < 0.05)
write.table(meancheckflarent2, "meancheckflarent2pvalues1000.csv", sep=",")

meancheckflaret2 <- meancheck("ADFlareT.txt", v2flaret, 3000)
mean(meancheckflaret2 < 0.05)
write.table(meancheckflaret2, "meancheckflaret2pvalues1000.csv", sep=",")

meancheckpostflare2 <- meancheck("ADPostFlare.txt", v2postflare, 3000)
mean(meancheckpostflare2 < 0.05)
write.table(meancheckpostflare2, "meancheckpostflare2pvalues1000.csv", sep=",")

###Varcheck  
varcheckcontrol2 <- variancecheck("ADControls.txt",v2control,3000)
mean(varcheckcontrol2 < 0.05)
mean(varcheckcontrol2)
median(varcheckcontrol2)
write.table(varcheckcontrol2, "varcheckcontrol2pvalues1000.csv", sep=",")

varcheckbaseline2 <- variancecheck("ADBaseline.txt",v2baseline,3000)
mean(varcheckbaseline2 < 0.05)
mean(varcheckbaseline2)
median(varcheckbaseline2)
write.table(varcheckbaseline2, "varcheckbaseline2pvalues1000.csv", sep=",")

varcheckflarent2 <- variancecheck("ADFlareNT.txt",v2flarent,3000)
mean(varcheckflarent2 < 0.05)
mean(varcheckflarent2)
median(varcheckflarent2)
write.table(varcheckflarent2, "varcheckflarent2pvalues1000.csv", sep=",")

varcheckflaret2 <- variancecheck("ADFlareT.txt", v2flaret, 3000)
mean(varcheckflaret2 < 0.05)
mean(varcheckflaret2)
median(varcheckflaret2)
write.table(varcheckflaret2, "varcheckflaret2pvalues1000.csv", sep=",")

varcheckpostflare2 <- variancecheck("ADPostFlare.txt", v2postflare, 3000)
mean(varcheckpostflare2 < 0.05)
mean(varcheckpostflare2)
median(varcheckpostflare2)
write.table(varcheckpostflare2, "varcheckpostflare2pvalues1000.csv", sep=",")
