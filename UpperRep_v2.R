####Loading Package####

library(PearsonDS)
library(HMP)

#########################################################################
#########################################################################
#######################Fall/Winter#Data#Simulation#######################
#########################################################################
#########################################################################

Fallmean <- read.table("UpperRFall.txt")
FallSD <- read.table("FallSD.txt")
SDDimRow <- dimnames(FallSD)[[2]]

FallSD <- matrix(FallSD,nrow=13,ncol=1,dimnames=list(SDDimRow,"SD"))
mode(FallSD) <- "numeric"
FallSD <- FallSD/100

###Function for calculating percent remainder
###Format for original data: 
###Rows: Samples/Subjects, Columns: Taxa

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

PRemainder <- Premainder(Fallmean)
FallPSD <- cbind(PRemainder,FallSD)
mode(FallPSD) <- "numeric"

##Mean = alpha/alpha+beta
##Variance = alpha*beta/(alpha+beta)^2 * (alpha + beta + 1)
##SD^2 = alpha*beta/(alpha+beta)^2 * (alpha + beta + 1)
## a= alpha b=beta
##using WolframAlpha to fine alpha and beta

getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

FallparProteo <- getBetaParams(FallPSD[1,1], FallPSD[1,2])
FallparFirm <- getBetaParams(FallPSD[2,1], FallPSD[2,2])
FallparBact <- getBetaParams(FallPSD[3,1], FallPSD[3,2])
FallparActino <- getBetaParams(FallPSD[4,1], FallPSD[4,2])
FallparFuso <- getBetaParams(FallPSD[5,1], FallPSD[5,2])
FallparCyan <- getBetaParams(FallPSD[6,1], FallPSD[6,2])
FallparOD1 <- getBetaParams(FallPSD[7,1], FallPSD[7,2])
FallparTM7 <- getBetaParams(FallPSD[8,1], FallPSD[8,2])
FallparDT <- getBetaParams(FallPSD[9,1], FallPSD[9,2])
FallparNitro <- getBetaParams(FallPSD[10,1], FallPSD[10,2])
FallparPlanct <- getBetaParams(FallPSD[11,1], FallPSD[11,2])

###Loop
numrow <- 50
numcol <- dim(PRemainder)[1]
UpperRep <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRep) <- dimnames(PRemainder)[[1]]
rownames(UpperRep) <- rownames(UpperRep, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRep)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRep[i,1] <- rpearson(n=1, params=FallparProteo) 
  total <- total - UpperRep[i,1]
  
  ###Firmicutes
  r <- rpearson(n=1, params=FallparFirm)
  UpperRep[i,2] <- total*r
  total <- total - UpperRep[i,2]
  
  ###Bacteroidetes
  r <- rpearson(n=1, params=FallparBact)
  UpperRep[i,3] <- total*r
  total <- total - UpperRep[i,3]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=FallparActino)
  UpperRep[i,4] <- total*r
  total <- total - UpperRep[i,4]
  
  ###Fusobacteria
  r <- rpearson(n=1, params=FallparFuso)
  UpperRep[i,5] <- total*r
  total <- total - UpperRep[i,5]
  
  ###Cyanobacteria
  r <- rpearson(n=1, params=FallparCyan)
  UpperRep[i,6] <- total*r
  total <- total - UpperRep[i,6]
  
  ###OD1
  r <- rpearson(n=1, params=FallparOD1)
  UpperRep[i,7] <- total*r
  total <- total - UpperRep[i,7]
  
  ###TM7
  r <- rpearson(n=1, params=FallparTM7)
  UpperRep[i,8] <- total*r
  total <- total - UpperRep[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearson(n=1, params=FallparDT)
  UpperRep[i,9] <- total*r
  total <- total - UpperRep[i,9]
  
  ###Nitrospira
  r <- rpearson(n=1, params=FallparNitro)
  UpperRep[i,10] <- total*r
  total <- total - UpperRep[i,10]
  
  ###Planctomycetes
  r <- rpearson(n=1, params=FallparPlanct)
  UpperRep[i,11] <- total*r
  
  ###Chloroflexi
  UpperRep[i,12] <- total - UpperRep[i,11]
  
  ###BRC1
  UpperRep[1:numrow,13] <- 0
}

write.table(x=UpperRep,file="UpperRepFallData.csv",sep=",")

#barchart(x=UpperRep,horizontal=FALSE, col=rainbow(13))

##########################################################################
##########################################################################
########################Effect#size#data#Simulation#######################
##########################################################################
##########################################################################


set.seed(123)

Fall <- read.csv("pUpperRepFall.csv")
Spring <- read.csv("pUpperRepSpring.csv")

Fallavg <- matrix(data=colMeans(Fall), nrow=1, ncol=13)
colnames(Fallavg) <- colnames(Fall)

S1 <- Fallavg[1] - 0.7

sum3taxa <- sum(Fallavg[,2:4])
Firmspace <- Fallavg[,2]/sum3taxa
Bactspace <- Fallavg[,3]/sum3taxa
Actinospace <- Fallavg[,4]/sum3taxa

ProteoS1 <- 0.7
FirmS1 <- Firmspace * S1
BactS1 <- Bactspace * S1
ActinoS1 <- Actinospace * S1

S1mean <- Fallavg
S1mean[,1] <- ProteoS1
S1mean[,2] <- FirmS1
S1mean[,3] <- BactS1
S1mean[,4] <- ActinoS1

S2 <- Fallavg[,1] - 0.65

ProteoS2 <- 0.65
FirmS2 <- Firmspace * S2
BactS2 <- Bactspace * S2
ActinoS2 <- Actinospace * S2

S2mean <- Fallavg
S2mean[,1] <- ProteoS2
S2mean[,2] <- FirmS2
S2mean[,3] <- BactS2
S2mean[,4] <- ActinoS2

S3 <- Fallavg[,1] - 0.6

ProteoS3 <- 0.60
FirmS3 <- Firmspace * S3
BactS3 <- Bactspace * S3
ActinoS3 <- Actinospace * S3

S3mean <- Fallavg
S3mean[,1] <- ProteoS3
S3mean[,2] <- FirmS3
S3mean[,3] <- BactS3
S3mean[,4] <- ActinoS3

S4 <- Fallavg[,1] - 0.55

ProteoS4 <- 0.55
FirmS4 <- Firmspace * S4
BactS4 <- Bactspace * S4
ActinoS4 <- Actinospace * S4

S4mean <- Fallavg
S4mean[,1] <- ProteoS4
S4mean[,2] <- FirmS4
S4mean[,3] <- BactS4
S4mean[,4] <- ActinoS4

S5 <- Fallavg[,1] - 0.50

ProteoS5 <- 0.50
FirmS5 <- Firmspace * S5
BactS5 <- Bactspace * S5
ActinoS5 <- Actinospace * S5

S5mean <- Fallavg
S5mean[,1] <- ProteoS5
S5mean[,2] <- FirmS5
S5mean[,3] <- BactS5
S5mean[,4] <- ActinoS5

S1mean
S2mean
S3mean
S4mean
S5mean

FallSD <- read.table("FallSD.txt")
SpringSD <- read.table("SpringSD.txt")

FirmSD <- (SpringSD[2] - FallSD[2])/5
BactSD <- (SpringSD[3] - FallSD[3])/5
ActinoSD <- (SpringSD[4] - FallSD[4])/5

###Proteobacteria SD stays same

FirmSDS1 <- FallSD[2] + FirmSD
FirmSDS2 <- FirmSDS1 + FirmSD
FirmSDS3 <- FirmSDS2 + FirmSD
FirmSDS4 <- FirmSDS3 + FirmSD
FirmSDS5 <- FirmSDS4 + FirmSD

BactSDS1 <- FallSD[3] + BactSD
BactSDS2 <- BactSDS1 + BactSD
BactSDS3 <- BactSDS2 + BactSD
BactSDS4 <- BactSDS3 + BactSD
BactSDS5 <- BactSDS4 + BactSD

ActinoSDS1 <- FallSD[4] + ActinoSD
ActinoSDS2 <- ActinoSDS1 + ActinoSD
ActinoSDS3 <- ActinoSDS2 + ActinoSD
ActinoSDS4 <- ActinoSDS3 + ActinoSD
ActinoSDS5 <- ActinoSDS4 + ActinoSD

SDS1 <- matrix(data=c(FallSD[1], FirmSDS1, BactSDS1, ActinoSDS1, FallSD[5], FallSD[6], FallSD[7], 
                   FallSD[8], FallSD[9], FallSD[10], FallSD[11], FallSD[12], FallSD[13]), dimnames=list(colnames(Fall), "SD"))
SDS2 <- matrix(data=c(FallSD[1], FirmSDS2, BactSDS2, ActinoSDS2, FallSD[5], FallSD[6], FallSD[7],
                   FallSD[8], FallSD[9], FallSD[10], FallSD[11], FallSD[12], FallSD[13]), dimnames=list(colnames(Fall), "SD"))
SDS3 <- matrix(data=c(FallSD[1], FirmSDS3, BactSDS3, ActinoSDS3, FallSD[5], FallSD[6], FallSD[7],
                   FallSD[8], FallSD[9], FallSD[10], FallSD[11], FallSD[12], FallSD[13]),  dimnames=list(colnames(Fall), "SD"))
SDS4 <- matrix(data=c(FallSD[1], FirmSDS4, BactSDS4, ActinoSDS4, FallSD[5], FallSD[6], FallSD[7],
                   FallSD[8], FallSD[9], FallSD[10], FallSD[11], FallSD[12], FallSD[13]), dimnames=list(colnames(Fall), "SD"))
SDS5 <- matrix(data=c(FallSD[1], FirmSDS5, BactSDS5, ActinoSDS5, FallSD[5], FallSD[6], FallSD[7],
                   FallSD[8], FallSD[9], FallSD[10], FallSD[11], FallSD[12], FallSD[13]), dimnames=list(colnames(Fall), "SD"))

mode(SDS1) <- "numeric"
mode(SDS2) <- "numeric"
mode(SDS3) <- "numeric"
mode(SDS4) <- "numeric"
mode(SDS5) <- "numeric"

SDS1 <- SDS1/100
SDS2 <- SDS2/100
SDS3 <- SDS3/100
SDS4 <- SDS4/100
SDS5 <- SDS5/100

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

PUpperFallS1 <- Premainder((S1mean*100))
PUpperFallS2 <- Premainder((S2mean*100))
PUpperFallS3 <- Premainder((S3mean*100))
PUpperFallS4 <- Premainder((S4mean*100))
PUpperFallS5 <- Premainder((S5mean*100))

S1PSD <- cbind(PUpperFallS1,SDS1)
S2PSD <- cbind(PUpperFallS2, SDS2)
S3PSD <- cbind(PUpperFallS3, SDS3)
S4PSD <- cbind(PUpperFallS4, SDS4)
S5PSD <- cbind(PUpperFallS5, SDS5)

######
getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}
######

S1parProteo <- getBetaParams(S1PSD[1,1], S1PSD[1,2])
S1parFirm <- getBetaParams(S1PSD[2,1], S1PSD[2,2])
S1parBact <- getBetaParams(S1PSD[3,1], S1PSD[3,2])
S1parActino <- getBetaParams(S1PSD[4,1], S1PSD[4,2])

S2parProteo <- getBetaParams(S2PSD[1,1], S2PSD[1,2])
S2parFirm <- getBetaParams(S2PSD[2,1], S2PSD[2,2])
S2parBact <- getBetaParams(S2PSD[3,1], S2PSD[3,2])
S2parActino <- getBetaParams(S2PSD[4,1], S2PSD[4,2])

S3parProteo <- getBetaParams(S3PSD[1,1], S3PSD[1,2])
S3parFirm <- getBetaParams(S3PSD[2,1], S3PSD[2,2])
S3parBact <- getBetaParams(S3PSD[3,1], S3PSD[3,2])
S3parActino <- getBetaParams(S3PSD[4,1], S3PSD[4,2])

S4parProteo <- getBetaParams(S4PSD[1,1], S4PSD[1,2])
S4parFirm <- getBetaParams(S4PSD[2,1], S4PSD[2,2])
S4parBact <- getBetaParams(S4PSD[3,1], S4PSD[3,2])
S4parActino <- getBetaParams(S4PSD[4,1], S4PSD[4,2])

S5parProteo <- getBetaParams(S5PSD[1,1], S5PSD[1,2])
S5parFirm <- getBetaParams(S5PSD[2,1], S5PSD[2,2])
S5parBact <- getBetaParams(S5PSD[3,1], S5PSD[3,2])
S5parActino <- getBetaParams(S5PSD[4,1], S5PSD[4,2])


###### ######
##S1#####70##
###### ######
numrow <- 50
numcol <- dim(PUpperFallS1)[1]
UpperRepS1 <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRepS1) <- dimnames(PUpperFallS1)[[1]]
rownames(UpperRepS1) <- rownames(UpperRepS1, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRepS1)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRepS1[i,1] <- rpearson(n=1, params=S1parProteo)
  total <- total - UpperRepS1[i,1]
  
  ###Firmicutes
  r <- rpearson(n=1, params=S1parFirm)
  UpperRepS1[i,2] <- total*r
  total <- total - UpperRepS1[i,2]
  
  
  ###Bacteroidetes
  r <- rpearson(n=1, params=S1parBact)
  UpperRepS1[i,3] <- total*r
  total <- total - UpperRepS1[i,3]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=S1parActino)
  UpperRepS1[i,4] <- total*r
  total <- total - UpperRepS1[i,4]
  
  ###Fusobacteria
  r <- rpearson(n=1, params=FallparFuso)
  UpperRepS1[i,5] <- total*r
  total <- total - UpperRepS1[i,5]
  
  ###Cyanobacteria
  r <- rpearson(n=1, params=FallparCyan)
  UpperRepS1[i,6] <- total*r
  total <- total - UpperRepS1[i,6]
  
  ###OD1
  r <- rpearson(n=1, params=FallparOD1)
  UpperRepS1[i,7] <- total*r
  total <- total - UpperRepS1[i,7]
  
  ###TM7
  r <- rpearson(n=1, params=FallparTM7)
  UpperRepS1[i,8] <- total*r
  total <- total - UpperRepS1[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearson(n=1, params=FallparDT)
  UpperRepS1[i,9] <- total*r
  total <- total - UpperRepS1[i,9]
  
  ###Nitrospira
  r <- rpearson(n=1, params=FallparNitro)
  UpperRepS1[i,10] <- total*r
  total <- total - UpperRepS1[i,10]
  
  ###Planctomycetes
  r <- rpearson(n=1, params=FallparPlanct)
  UpperRepS1[i,11] <- total*r
  
  ###Chloroflexi
  UpperRepS1[i,12] <- total - UpperRepS1[i,11]
  
  ###BRC1
  UpperRepS1[1:numrow,13] <- 0
}

colMeans(UpperRepS1)
barchart(x=UpperRepS1,horizontal=FALSE, col=rainbow(13))
write.table(x=UpperRepS1, file="UpperRepS1.csv", sep=",")

###### ######
##S2#####65##
###### ######
numrow <- 50
numcol <- dim(PUpperFallS2)[1]
UpperRepS2 <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRepS2) <- dimnames(PUpperFallS2)[[1]]
rownames(UpperRepS2) <- rownames(UpperRepS2, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRepS2)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRepS2[i,1] <- rpearson(n=1, params=S2parProteo) 
  total <- total - UpperRepS2[i,1]
  
  ###Firmicutes
  r <- rpearson(n=1, params=S2parFirm)
  UpperRepS2[i,2] <- total*r
  total <- total - UpperRepS2[i,2]
  
  
  ###Bacteroidetes
  r <- rpearson(n=1, params=S2parBact)
  UpperRepS2[i,3] <- total*r
  total <- total - UpperRepS2[i,3]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=S2parActino)
  UpperRepS2[i,4] <- total*r
  total <- total - UpperRepS2[i,4]
  
  ###Fusobacteria
  r <- rpearson(n=1, params=FallparFuso)
  UpperRepS2[i,5] <- total*r
  total <- total - UpperRepS2[i,5]
  
  ###Cyanobacteria
  r <- rpearson(n=1, params=FallparCyan)
  UpperRepS2[i,6] <- total*r
  total <- total - UpperRepS2[i,6]
  
  ###OD1
  r <- rpearson(n=1, params=FallparOD1)
  UpperRepS2[i,7] <- total*r
  total <- total - UpperRepS2[i,7]
  
  ###TM7
  r <- rpearson(n=1, params=FallparTM7)
  UpperRepS2[i,8] <- total*r
  total <- total - UpperRepS2[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearson(n=1, params=FallparDT)
  UpperRepS2[i,9] <- total*r
  total <- total - UpperRepS2[i,9]
  
  ###Nitrospira
  r <- rpearson(n=1, params=FallparNitro)
  UpperRepS2[i,10] <- total*r
  total <- total - UpperRepS2[i,10]
  
  ###Planctomycetes
  r <- rpearson(n=1, params=FallparPlanct)
  UpperRepS2[i,11] <- total*r
  
  ###Chloroflexi
  UpperRepS2[i,12] <- total - UpperRepS2[i,11]
  
  ###BRC1
  UpperRepS2[1:numrow,13] <- 0
}

colMeans(UpperRepS2)
barchart(x=UpperRepS2,horizontal=FALSE, col=rainbow(13))
write.table(x=UpperRepS2, file="UpperRepS2.csv", sep=",")

###### ######
##S3#####60##
###### ######
numrow <- 50
numcol <- dim(PUpperFallS3)[1]
UpperRepS3 <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRepS3) <- dimnames(PUpperFallS3)[[1]]
rownames(UpperRepS3) <- rownames(UpperRepS3, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRepS3)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRepS3[i,1] <- rpearson(n=1, params=S3parProteo) 
  total <- total - UpperRepS3[i,1]
  
  ###Firmicutes
  r <- rpearson(n=1, params=S3parFirm)
  UpperRepS3[i,2] <- total*r
  total <- total - UpperRepS3[i,2]
  
  
  ###Bacteroidetes
  r <- rpearson(n=1, params=S3parBact)
  UpperRepS3[i,3] <- total*r
  total <- total - UpperRepS3[i,3]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=S3parActino)
  UpperRepS3[i,4] <- total*r
  total <- total - UpperRepS3[i,4]
  
  ###Fusobacteria
  r <- rpearson(n=1, params=FallparFuso)
  UpperRepS3[i,5] <- total*r
  total <- total - UpperRepS3[i,5]
  
  ###Cyanobacteria
  r <- rpearson(n=1, params=FallparCyan)
  UpperRepS3[i,6] <- total*r
  total <- total - UpperRepS3[i,6]
  
  ###OD1
  r <- rpearson(n=1, params=FallparOD1)
  UpperRepS3[i,7] <- total*r
  total <- total - UpperRepS3[i,7]
  
  ###TM7
  r <- rpearson(n=1, params=FallparTM7)
  UpperRepS3[i,8] <- total*r
  total <- total - UpperRepS3[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearson(n=1, params=FallparDT)
  UpperRepS3[i,9] <- total*r
  total <- total - UpperRepS3[i,9]
  
  ###Nitrospira
  r <- rpearson(n=1, params=FallparNitro)
  UpperRepS3[i,10] <- total*r
  total <- total - UpperRepS3[i,10]
  
  ###Planctomycetes
  r <- rpearson(n=1, params=FallparPlanct)
  UpperRepS3[i,11] <- total*r
  
  ###Chloroflexi
  UpperRepS3[i,12] <- total - UpperRepS3[i,11]
  
  ###BRC1
  UpperRepS3[1:numrow,13] <- 0
}

colMeans(UpperRepS3)
barchart(x=UpperRepS3,horizontal=FALSE, col=rainbow(13))
write.table(x=UpperRepS3, file="UpperRepS3.csv", sep=",")

###### ######
##S4#####55##
###### ######
numrow <- 50
numcol <- dim(PUpperFallS4)[1]
UpperRepS4 <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRepS4) <- dimnames(PUpperFallS4)[[1]]
rownames(UpperRepS4) <- rownames(UpperRepS4, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRepS4)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRepS4[i,1] <- rpearson(n=1, params=S4parProteo) 
  total <- total - UpperRepS4[i,1]
  
  ###Firmicutes
  r <- rpearson(n=1, params=S4parFirm)
  UpperRepS4[i,2] <- total*r
  total <- total - UpperRepS4[i,2]
  
  
  ###Bacteroidetes
  r <- rpearson(n=1, params=S4parBact)
  UpperRepS4[i,3] <- total*r
  total <- total - UpperRepS4[i,3]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=S4parActino)
  UpperRepS4[i,4] <- total*r
  total <- total - UpperRepS4[i,4]
  
  ###Fusobacteria
  r <- rpearson(n=1, params=FallparFuso)
  UpperRepS4[i,5] <- total*r
  total <- total - UpperRepS4[i,5]
  
  ###Cyanobacteria
  r <- rpearson(n=1, params=FallparCyan)
  UpperRepS4[i,6] <- total*r
  total <- total - UpperRepS4[i,6]
  
  ###OD1
  r <- rpearson(n=1, params=FallparOD1)
  UpperRepS4[i,7] <- total*r
  total <- total - UpperRepS4[i,7]
  
  ###TM7
  r <- rpearson(n=1, params=FallparTM7)
  UpperRepS4[i,8] <- total*r
  total <- total - UpperRepS4[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearson(n=1, params=FallparDT)
  UpperRepS4[i,9] <- total*r
  total <- total - UpperRepS4[i,9]
  
  ###Nitrospira
  r <- rpearson(n=1, params=FallparNitro)
  UpperRepS4[i,10] <- total*r
  total <- total - UpperRepS4[i,10]
  
  ###Planctomycetes
  r <- rpearson(n=1, params=FallparPlanct)
  UpperRepS4[i,11] <- total*r
  
  ###Chloroflexi
  UpperRepS4[i,12] <- total - UpperRepS4[i,11]
  
  ###BRC1
  UpperRepS4[1:numrow,13] <- 0
}

colMeans(UpperRepS4)
barchart(x=UpperRepS4,horizontal=FALSE, col=rainbow(13))
write.table(x=UpperRepS4, file="UpperRepS4.csv", sep=",")

###### ######
##S5#####50##
###### ######
numrow <- 50
numcol <- dim(PUpperFallS5)[1]
UpperRepS5 <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRepS5) <- dimnames(PUpperFallS5)[[1]]
rownames(UpperRepS5) <- rownames(UpperRepS5, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRepS5)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRepS5[i,1] <- rpearson(n=1, params=S5parProteo) 
  total <- total - UpperRepS5[i,1]
  
  ###Firmicutes
  r <- rpearson(n=1, params=S5parFirm)
  UpperRepS5[i,2] <- total*r
  total <- total - UpperRepS5[i,2]
  
  
  ###Bacteroidetes
  r <- rpearson(n=1, params=S5parBact)
  UpperRepS5[i,3] <- total*r
  total <- total - UpperRepS5[i,3]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=S5parActino)
  UpperRepS5[i,4] <- total*r
  total <- total - UpperRepS5[i,4]
  
  ###Fusobacteria
  r <- rpearson(n=1, params=FallparFuso)
  UpperRepS5[i,5] <- total*r
  total <- total - UpperRepS5[i,5]
  
  ###Cyanobacteria
  r <- rpearson(n=1, params=FallparCyan)
  UpperRepS5[i,6] <- total*r
  total <- total - UpperRepS5[i,6]
  
  ###OD1
  r <- rpearson(n=1, params=FallparOD1)
  UpperRepS5[i,7] <- total*r
  total <- total - UpperRepS5[i,7]
  
  ###TM7
  r <- rpearson(n=1, params=FallparTM7)
  UpperRepS5[i,8] <- total*r
  total <- total - UpperRepS5[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearson(n=1, params=FallparDT)
  UpperRepS5[i,9] <- total*r
  total <- total - UpperRepS5[i,9]
  
  ###Nitrospira
  r <- rpearson(n=1, params=FallparNitro)
  UpperRepS5[i,10] <- total*r
  total <- total - UpperRepS5[i,10]
  
  ###Planctomycetes
  r <- rpearson(n=1, params=FallparPlanct)
  UpperRepS5[i,11] <- total*r
  
  ###Chloroflexi
  UpperRepS5[i,12] <- total - UpperRepS5[i,11]
  
  ###BRC1
  UpperRepS5[1:numrow,13] <- 0
}

colMeans(UpperRepS5)
barchart(x=UpperRepS5,horizontal=FALSE, col=rainbow(13))
write.table(x=UpperRepS5, file="UpperRepS5.csv", sep=",")

