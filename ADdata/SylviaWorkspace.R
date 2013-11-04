set.seed(1234)
library(PearsonDS)
library(HMP)

(Control <- read.table("ADControls.txt"))
(Baseline <- read.table("ADBaseline.txt"))
(FlareNT <- read.table("ADFlareNT.txt"))
(FlareT <- read.table("ADFlareT.txt"))
(PostFlare <- read.table("ADPostFlare.txt"))

colSums(Control)
colSums(Baseline)
colSums(FlareNT)
colSums(FlareT)
colSums(PostFlare)

#######
(ControlSD <- matrix(apply(Control,1,sd),dimnames=list(rownames(Control), "SD")))
(BaselineSD <- matrix(apply(Baseline,1,sd), dimnames=list(rownames(Baseline), "SD")))
(FlareNTSD <- matrix(apply(FlareNT,1,sd), dimnames=list(rownames(FlareNT), "SD")))
(FlareTSD <- matrix(apply(FlareT,1,sd), dimnames=list(rownames(FlareT), "SD")))
(PostFlareSD <- matrix(apply(PostFlare,1,sd), dimnames=list(rownames(PostFlare), "SD")))

(ControlMean <- matrix(rowMeans(Control), dimnames=list(rownames(Control), "Mean")))
(BaselineMean <- matrix(rowMeans(Baseline), dimnames=list(rownames(Baseline), "Mean")))
(FlareNTMean <- matrix(rowMeans(FlareNT), dimnames=list(rownames(FlareNT), "Mean")))
(FlareTMean <- matrix(rowMeans(FlareT), dimname=list(rownames(FlareT), "Mean")))
(PostFlareMean <- matrix(rowMeans(PostFlare), dimnames=list(rownames(PostFlare), "Mean")))

(ControlMSD <- cbind(ControlMean, ControlSD))
(ControlMSD <- ControlMSD*100)
(BaselineMSD <- cbind(BaselineMean, BaselineSD))
(BaselineMSD <- BaselineMSD*100)
(FlareNTMSD <- cbind(FlareNTMean, FlareNTSD))
(FlareNTMSD <- FlareNTMSD*100)
(FlareTMSD <- cbind(FlareTMean, FlareTSD))
(FlareTMSD <- FlareTMSD*100)
(PostFlareMSD <- cbind(PostFlareMean, PostFlareSD))
(PostFlareMSD <- PostFlareMSD*100)

####Beta Dist before

#####
pMSD <- function(x){
  (nrow <- dim(x)[1])
  (ncol <- dim(x)[2])
  (y <- matrix(nrow=nrow,ncol=ncol, dimnames=list(rownames(x), colnames(x))))
  (total <- 100)
  (y[1,2] <- (x[1,2])/100)
  (y[1,1] <- (x[1,1])/100)
  for(i in 2:nrow){
    (y[i,2] <- x[i,2]/(total - x[i-1,1]))
    (total <- total - x[i-1,1])
    (y[i,1] <- x[i,1]/total)
  }
  return(y)
}

#####

#####
(nControl <- pMSD(ControlMSD))
(nBaseline <- pMSD(BaselineMSD))
(nFlareNT <- pMSD(FlareNTMSD))
(nFlareT <- pMSD(FlareTMSD))
(nPostFlare <- pMSD(PostFlareMSD))

####
pcontroltn <- urnorm(n=1000, mean=nControl[2,1], sd=nControl[2,2],lb=0,ub=1)
hist(pcontroltn,100)
sd(pcontroltn)
mean(pcontroltn)
####
numrow <- 25
numcol <- dim(ControlMSD)[1]
Cdata <- matrix(nrow=numrow,ncol=numcol)
colnames(Cdata) <- dimnames(ControlMSD)[[1]]
rownames(Cdata) <- rownames(Cdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(Cdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  Cdata[i,1] <- urnorm(n=1,mean=nControl[1,1],sd=nControl[1,2],lb=0,ub=1)
  total <- total - Cdata[i,1]
  
  ###Actinobacteria
  r <- urnorm(n=1, mean=nControl[2,1], sd=nControl[2,2], lb=0, ub=1)
  Cdata[i,2] <- total*r
  total <- total - Cdata[i,2]
  
  ###Proteobacteria
  r <- urnorm(n=1, mean=nControl[3,1], sd=nControl[3,2], lb=0, ub=1)
  Cdata[i,3] <- total*r
  total <- total - Cdata[i,3]
  
  ###Bacteriodetes
  r <- urnorm(n=1, mean=nControl[4,1], sd=nControl[4,2], lb=0, ub=1)
  Cdata[i,4] <- total*r
  total <- total - Cdata[i,4]
  
  ###Other
  Cdata[i,5] <- total
  
}


Cdata
colMeans(Cdata)
apply(Cdata,2,sd)
ControlMSD

#####
getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

########Parameters

(CparFirm <- getBetaParams(nControl[1,1], nControl[1,2]))
rpearson(n=1, params=CparFirm)
(CparActino <- getBetaParams(nControl[2,1], nControl[2,2]))
rpearson(n=1, params=CparActino)
(CparProteo <- getBetaParams(nControl[3,1], nControl[3,2]))
rpearson(n=1, params=CparProteo)
(CparBact <- getBetaParams(nControl[4,1], nControl[4,2]))
rpearson(n=1, params=CparBact)

(BparFirm <- getBetaParams(nBaseline[1,1], nBaseline[1,2]))
rpearson(n=1, params=BparFirm)
(BparActino <- getBetaParams(nBaseline[2,1], nBaseline[2,2]))
rpearson(n=1, params=BparActino)
(BparProteo <- getBetaParams(nBaseline[3,1], nBaseline[3,2]))
rpearson(n=1, params=BparProteo)
(BparBact <- getBetaParams(nBaseline[4,1], nBaseline[4,2]))
rpearson(n=1, params=BparBact)

(FNTparFirm <- getBetaParams(nFlareNT[1,1], nFlareNT[1,2]))
rpearson(n=1, params=FNTparFirm)
(FNTparActino <- getBetaParams(nFlareNT[2,1], nFlareNT[2,2]))
rpearson(n=1, params=FNTparActino)
(FNTparProteo <- getBetaParams(nFlareNT[3,1], nFlareNT[3,2]))
rpearson(n=1, params=FNTparProteo)
(FNTparBact <- getBetaParams(nFlareNT[4,1], nFlareNT[4,2]))
rpearson(n=1, params=FNTparBact)

(FTparFirm <- getBetaParams(nFlareT[1,1], nFlareT[1,2]))
rpearson(n=1, params=FTparFirm)
(FTparActino <- getBetaParams(nFlareT[2,1], nFlareT[2,2]))
rpearson(n=1, params=FTparActino)
(FTparProteo <- getBetaParams(nFlareT[3,1], nFlareT[3,2]))
rpearson(n=1, params=FTparProteo)
(FTparBact <- getBetaParams(nFlareT[4,1], nFlareT[4,2]))
rpearson(n=1, params=FTparBact)

(PFparFirm <- getBetaParams(nPostFlare[1,1], nPostFlare[1,2]))
rpearson(n=1, params=PFparFirm)
(PFparActino <- getBetaParams(nPostFlare[2,1], nPostFlare[2,2]))
rpearson(n=1, params=PFparActino)
(PFparProteo <- getBetaParams(nPostFlare[3,1], nPostFlare[3,2]))
rpearson(n=1, params=PFparProteo)
(PFparBact <- getBetaParams(nPostFlare[4,1], nPostFlare[4,2]))
rpearson(n=1, params=PFparBact)

####PostFlare
numrow <- 25
numcol <- dim(PostFlareMSD)[1]
PFdata <- matrix(nrow=numrow,ncol=numcol)
colnames(PFdata) <- dimnames(PostFlareMSD)[[1]]
rownames(PFdata) <- rownames(PFdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(PFdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  PFdata[i,1] <- rpearson(n=1, params=PFparFirm)
  total <- total - PFdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=PFparActino)
  PFdata[i,2] <- total*r
  total <- total - PFdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=PFparProteo)
  PFdata[i,3] <- total*r
  total <- total - PFdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=PFparBact)
  PFdata[i,4] <- total*r
  total <- total - PFdata[i,4]
  
  ###Other
  PFdata[i,5] <- total
  
}
colMeans(PFdata)
apply(PFdata,2,sd)
PostFlareMSD







####Control
numrow <- 25
numcol <- dim(ControlMSD)[1]
Cdata <- matrix(nrow=numrow,ncol=numcol)
colnames(Cdata) <- dimnames(ControlMSD)[[1]]
rownames(Cdata) <- rownames(Cdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(Cdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  Cdata[i,1] <- rpearson(n=1, params=CparFirm)
  total <- total - Cdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=CparActino)
  Cdata[i,2] <- total*r
  total <- total - Cdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=CparProteo)
  Cdata[i,3] <- total*r
  total <- total - Cdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=CparBact)
  Cdata[i,4] <- total*r
  total <- total - Cdata[i,4]
  
  ###Other
  Cdata[i,5] <- total
  
}
Cdata




































####15OCT2013
################################################################################
###Create list of Simulation
listdataset <- function(inputFilename,outputLabel,numdataset, numberSubjects){
  sim <- replicate(n=numdataset,expr=list())
  for(i in 1:length(sim)){
    sim[[i]] <- simulateBrokenStick(inputFilename,outputLabel, numberSubjects)
  }
  return(sim)
}

####
sim <- listdataset(inputFilename="ADControls.txt",outputLabel="Controls",numdataset=10,numberSubjects=25)
getMeanSD(sim)
####

###Mean and SD of all 
getMeanSD <- function(simdatalist){
  MSD <- replicate(n=length(simdatalist), expr=list())
  for(i in 1:10){
    Mean <- apply(simdatalist[[i]],2,mean)
    SD <- apply(simdatalist[[i]],2,sd)
    MSD[[i]] <- rbind(Mean, SD)
  }
  return(MSD)
}

###

###validation datasets
validationdatasets <- function(inputFilename, outputLabel, numberSimsubjects, numberSimulations, publishNumbersubjects, numberReads){ #7.Number of datasets to create/number of times to repeat
  
  BSData <- simulateBrokenStick(inputFilename,outputLabel, numberSimsubjects) #1.Simulate 10,000 samples
  listdataset <- replicate(n=length(numberSimulations*2), expr=list()) #2.Creates a place to store everything 
  datalist <- seq(from=1, to=numberSimulations*2, by=2)
  pilist <- seq(from=2, to=numberSimulations*2, by=2)
  
  for(k in 1:numberSimulations){
    publishNumbersubjects <- 22
    col <- dim(BSData)[2]
    table <- matrix(nrow=publishNumbersubjects,ncol=col)
    colnames(table) <- dimnames(BSData)[[2]]
    rownames(table) <- rownames(table, do.NULL= FALSE, prefix= "Sample")
    pick <- sample(x=1:10000,size=publishNumbersubjects) #3.Pick n sets of samples  
    for(i in 1:publishNumbersubjects){ #4.Put together as one dataset  
      table[i,] <- round(BSData[pick[i],]*numberReads)
      pio <- DM.MoM(table)$theta #5.Get La Rosa pi statistics for set
    }
    listdataset[datalist[k]:pilist[k]] <- list(table, pio) #Save the pi stat with dataset
  }
  return(listdataset)
}

###
###Update
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

###
control <- validationdatasets(inputFilename="ADControls.txt", outputLabel="Controls", numberSimsubjects=10000, numberSimulations=1000, publishNumbersubjects=22,numberReads=3000)
b <- validationdatasets(inputFilename="ADControls.txt", outputLabel="Controls", numberSimsubjects=10000, numberSimulations=10, publishNumbersubjects=22)
###

###
d <- seq(1,20,2)
t <- seq(2,2000,2)

###Gathers all theta values into a table
thetatable <- function(x){
  llength <- length(x)
  hllength <- llength/2
  tloc <- seq(from=2, to=llength, by=2)
  ttable <- matrix(nrow=hllength,ncol=1)
  colnames(ttable) <- "theta"
  for(i in 1:hllength){
    ttable[i] <- x[[tloc[i]]]
  }
  hist(ttable)
  return(ttable)
}

test <- thetatable(a)
####

####Gathers all simulated datasets
getDatatables <- function(x){
  dlength <- length(x)
  hdlength <- length(x)/2
  dloc <- seq(from=1, to=dlength, by=2 )
  dlist <- replicate(n=hdlength,expr=list())
  for(i in 1:hdlength){
    dlist[[i]] <- x[[dloc[i]]]
  }
  return(dlist)
}
####
pubControl <- dirmult(round(t(read.table("ADControls.txt"))*3000))

listcontrols <- getDatatables(control) 
listcontrols[2]
data <- listcontrols
(controlcheck <- Xmc.sevsample(listcontrols, pubControl$pi))
(check <- Xsc.onesample(listcontrols[1],pubControl$pi))
####

####
simcontrol <- simulateBrokenStick("ADControls.txt", "Control", 10000)
simcontrolreads <- round(simcontrol*3000)
pubControl <- round(t(read.table("ADControls.txt"))*3000)
fit.simcontrol <- dirmult(simcontrolreads)

listcontrol <- list(simcontrolreads,pubControl)
OVControlCheck <- Xoc.sevsample(listcontrol)

ControlCheck <- Xsc.onesample(pubControl,fit.simcontrol$pi)
BaselineCheck <- check(inputFilename="ADBaseline.txt", outputLabel="ADbaseline", numberSubjects=10000, numReads=3000)
FlareNTCheck <- check("ADFlareNT.txt", "ADFlareNT", 10000, 3000)
FlareTCheck <- check("ADFlareT.txt", "ADFlareT", 10000, 3000)
PostflareCheck <- check("ADPostFlare.txt","ADPostflare", 10000, 3000)
####
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
mean(varcheckbaseline < 0.05)
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











BaselineOCheck <- Ocheck(inputFilename="ADBaseline.txt", outputLabel="ADbaseline", numberSubjects=10000, numReads=3000)
FlareNToCheck <- Ocheck("ADFlareNT.txt", "ADFlareNT", 10000, 3000)
FlareToCheck <- Ocheck("ADFlareT.txt", "ADFlareT", 10000, 3000)
PostflareoCheck <- Ocheck("ADPostFlare.txt","ADPostflare", 10000, 3000)

####
a <- mean(test)
s <- sd(test)
n <- 22
xbar <- DM.MoM(t(read.table("ADControls.txt"))*3000)$theta
(z <- (xbar-a)/(s/sqrt(n)))
p <- 2*pnorm(-abs(z))
read.table("ADBaseline.txt")
simulateBrokenStick("ADBaseline.txt","ADBaseline", 300)
simulateBrokenStick("ADFlareNT.txt","ADFlare", 300)
simulateBrokenStick("ADFlareT.txt","ADTreatment", 300)
simulateBrokenStick("ADPostFlare.txt","ADPost")