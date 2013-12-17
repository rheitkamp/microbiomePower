set.seed(1234)
library(HMP)
library(PearsonDS)

####DATA
controls <- read.table("tADControl.txt")
simcontrols <- simulateBrokenStick(inputFilename="tADControl.txt",  outputLabel="Controls", numberSubjects=10000)
vdscontrols <- validationdatasets(simdata=simcontrols,numberSimulations=1000,publishNumbersubjects=22,numReads=3000)
write.table(x=simcontrols, file="simcontroldata10000.txt")
write.table(x=vdscontrols, file="validationcontroldatasets1000.txt")

baseline <- read.table("tADBaseline.txt")
simbaseline <- simulateBrokenStick(inputFilename="tADBaseline.txt", outputLabel="Baseline", numberSubjects=10000)
vdsbaseline2nd <- validationdatasets(simdata=simbaseline, numberSimulations=1000,publishNumbersubjects=12,numReads=3000)
write.table(x=simbaseline, file="simbaselinedata10000.txt")
write.table(x=vdsbaseline, file="validationbaselinedataset1000.txt")

flarent <- read.table("tADFlareNT.txt")
simflarent <- simulateBrokenStick(inputFilename="tADFlareNT.txt", outputLabel="FlareNT", numberSubjects=10000)
vdsflarent <- validationdatasets(simdata=simflarent,numberSimulations=1000,publishNumbersubjects=7,numReads=3000)
write.table(x=simflarent, file="simflarentdata10000.txt")
write.table(x=vdsflarent, file="validationflarentdataset1000.txt")

flaret <- read.table("tADFlareT.txt")
simflaret <- simulateBrokenStick(inputFilename="tADFlareT.txt", outputLabel="FlareT", numberSubjects=10000)
vdsflaret <- validationdatasets(simdata=simflaret,numberSimulations=1000,publishNumbersubjects=5,numReads=3000)
write.table(x=simflaret, file="simflaretdata10000.txt")
write.table(x=vdsflaret, file="validationflaretdataset1000.txt")

postflare <- read.table("tADPostflare.txt")
simpostflare <- simulateBrokenStick(inputFilename="tADPostflare.txt", outputLabel="Postflare", numberSubjects=10000)
vdspostflare <- validationdatasets(simdata=simpostflare,numberSimulations=1000,publishNumbersubjects=11,numReads=3000)
write.table(x=simpostflare, file="simpostflaredata10000.txt")
write.table(x=vdspostflare, file="validationpostflaredatasets.txt")

#######################Hypothesis Testing
tcontrol <- round(t(controls)*3000)
simcontrol3000 <- round(simcontrols*3000)
listcontrol <- list(tcontrol,simcontrol3000)
simvspubcontrolvars <- Xoc.sevsample(listcontrol)

tbaseline <- round(t(baseline)*3000)
simbaseline3000 <- round(simbaseline*3000)
listbaseline <- list(tbaseline,simbaseline3000)
simvspubbaselinevars <- Xoc.sevsample(listbaseline)

tflarent <- round(t(flarent)*3000)
simflarent3000 <- round(simflarent*3000)
listflarent <- list(tflarent, simflarent3000)
simvspubflarentvars <- Xoc.sevsample(listflarent)

tflaret <- round(t(flaret)*3000)
simflaret3000 <- round(simflaret*3000)
listflaret <- list(tflaret, simflarent3000)
simvspubflaretvars <- Xoc.sevsample(listflaret)

tpostflare <- round(t(postflare)*3000)
simpostflare3000 <- round(simpostflare*3000)
listpostflare <- list(tpostflare,simpostflare3000)
simvspubpostflarevars <- Xoc.sevsample(listpostflare)

####Mean check
meancheckcontrol <- meancheck("tADControl.txt",datalist=vdscontrols,numReads=3000)
mean(meancheckcontrol > 0.05)
mean(meancheckcontrol)
median(meancheckcontrol)
write.table(x=meancheckcontrol, "meanpvaluescontrols1000set.txt")

meancheckbaseline <- meancheck("tADBaseline.txt",datalist=vdsbaseline,numReads=3000)
mean(meancheckbaseline > 0.05)
mean(meancheckbaseline)
median(meancheckbaseline)
write.table(x=meancheckbaseline, "meanpvaluesbaseline1000set.txt")

meancheckflarent <- meancheck("tADFlareNT.txt",datalist=vdsflarent,numReads=3000)
mean(meancheckflarent > 0.05)
mean(meancheckflarent)
median(meancheckflarent)
write.table(x=meancheckflarent, "meanpvaluesflarent1000set.txt")

meancheckflaret <- meancheck("tADFlareT.txt",datalist=vdsflaret,numReads=3000)
mean(meancheckflaret > 0.05)
mean(meancheckflaret)
median(meancheckflaret)
write.table(x=meancheckflaret, "meanpvaluesflaret1000set.txt")

meancheckpostflare <- meancheck("tADPostflare.txt",datalist=vdspostflare,numReads=3000)
mean(meancheckpostflare > 0.05)
mean(meancheckpostflare)
median(meancheckpostflare)
write.table(x=meancheckpostflare, "meanpvaluespostflare1000set.txt")

#####Variance check
variancecheckcontrol <- variancecheck("tADControl.txt",datalist=vdscontrols,numReads=3000)
hist(variancecheckcontrol,100)
mean(variancecheckcontrol > 0.05)
mean(variancecheckcontrol)
median(variancecheckcontrol)
write.table(x=variancecheckcontrol, "variancepvaluescontrol1000set.txt")

variancecheckbaseline <- variancecheck("tADBaseline.txt",datalist=vdsbaseline,numReads=3000)
hist(variancecheckbaseline,100)
mean(variancecheckbaseline > 0.05)
mean(variancecheckbaseline)
median(variancecheckbaseline)
write.table(x=variancecheckbaseline, "variancepvaluesbaseline1000set.txt")

variancecheckflarent <- variancecheck("tADFlareNT.txt",datalist=vdsflarent, numReads=3000)
hist(variancecheckflarent,100)
mean(variancecheckflarent > 0.05)
mean(variancecheckflarent)
median(variancecheckflarent)
write.table(x=variancecheckflarent, "variancepvaluesflarent1000set.csv", sep=",")

variancecheckflaret <- variancecheck("tADFlareT.txt", datalist=vdsflaret, numReads=3000)
hist(variancecheckflaret,100)
mean(variancecheckflaret > 0.05)
mean(variancecheckflaret)
median(variancecheckflaret)
write.table(x=variancecheckflaret, "variancepvaluesflaret1000set.txt")

variancecheckpostflare <- variancecheck("tADPostflare.txt", datalist=vdsflaret, numReads=3000)
hist(variancecheckpostflare, 100)
mean(variancecheckpostflare > 0.05)
mean(variancecheckpostflare)
median(variancecheckpostflare)
write.table(x=variancecheckpostflare, "variancepvalues1000set.txt")

#####attempt two with this dataset
variancecheckbaseline <- variancecheck("tADBaseline.txt",datalist=vdsbaseline2nd,numReads=3000)
mean(variancecheckbaseline < 0.05)
mean(variancecheckbaseline)
median(variancecheckbaseline)
write.table(x=variancecheckbaseline, "variancepvaluesbaseline1000set.txt")

#########25OCT2013
v2scontrol <- validationdatasets2(inputFilename="tADControl.txt",outputLabel="Control",numberSubjects=22,numberSimulations=1000)
v2sbaseline <- validationdatasets2(inputFilename="tADBaseline.txt", outputLabel="Baseline",numberSubjects=12,numberSimulations=1000)
v2sflarent <- validationdatasets2(inputFilename="tADFlareNT.txt", outputLabel="FlareNT", numberSubjects=7, numberSimulations=1000)
v2sflaret <- validationdatasets2(inputFilename="tADFlareT.txt", outputLabel="FlareT", numberSubjects=5,numberSimulations=1000)
v2spostflare <- validationdatasets2(inputFilename="tADPostflare.txt",outputLabel="Postflare", numberSubjects=11,numberSimulations=1000)

meancheck2scontrol <- meancheck(inputFilename="tADControl.txt",datalist=v2scontrol,numReads=3000)
mean(meancheck2scontrol < 0.05)
mean(meancheck2scontrol)
median(meancheck2scontrol)

meancheck2sbaseline <- meancheck(inputFilename="tADBaseline.txt",datalist=v2sbaseline,numReads=3000)
mean(meancheck2sbaseline < 0.05)
mean(meancheck2sbaseline)
median(meancheck2sbaseline)

meancheck2flarent <- meancheck(inputFilename="tADFlareNT.txt",datalist=v2sflarent,numReads=3000)
mean(meancheck2sflarent < 0.05)
mean(meancheck2sflarent)
median(meancheck2sflarent)

meancheck2flaret <- meancheck(inputFilename="tADFlareT.txt",datalist=v2sflaret,numReads=3000)
meancheck2postflare <- meancheck(inputFilename="tADPostflare.txt",datalist=v2spostflare,numReads=3000)

varcheck2scontrol <- variancecheck("tADControl.txt",datalist=v2scontrol,numReads=3000)
mean(varcheck2scontrol < 0.05)
mean(varcheck2scontrol)
median(varcheck2scontrol)
varcheck2sbaseline <- variancecheck("tADBaseline.txt",datalist=v2sbaseline,numReads=3000)
mean(varcheck2sbaseline < 0.05)
mean(varcheck2sbaseline)
median(varcheck2sbaseline)
