library("PearsonDS")

###AD Controls %remaining dataset
control <- read.table("ADcontrolsredist.txt")

###0-1 AD Control extracted each individual taxa dataset (proportions)
Firm <- as.numeric(control[1,])
Actino <- as.numeric(control[2,])
Proteo <- as.numeric(control[3,])
Bacter <- as.numeric(control[4,])

###Distribution parameters 0-1 (proportions)
FirmADCparam <- pearsonFitML(Firm)
ActinoADCparam <- pearsonFitML(Actino)
ProteoADCparam <- pearsonFitML(Proteo)
BacterADCparam <- pearsonFitML(Bacter)

###Histogram of taxa distributions
hist(rpearson(n=1000, params=FirmADCparam), breaks=100)
hist(rpearson(n=1000, params=ActinoADCparam), breaks=100)
hist(rpearson(n=1000, params=ProteoADCparam), breaks=100)
hist(rpearson(n=1000, params=BacterADCparam), breaks=100)

###Generates AD Control dataset (proportions)###
numrow <- 100 ###Number of subjects/samples you want to generate
numcol <- 5 ###Number of taxa
ADControl <- matrix(nrow=numrow,ncol=numcol)
ADControl <- data.frame(ADControl)
colnames(ADControl) <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Others")
size <- dim(ADControl)[1]

for (i in 1:size){
  total <- 1
  ###Generate number for Firmicutes
  ADControl$Firmicutes[i] <- rpearson(n=1, params=FirmADCparam)
  total <- total - ADControl$Firmicutes[i]
  
  ###Generate number for Actinobacteria
  rActino <- rpearson(n=1, params=ActinoADCparam)
  ADControl$Actinobacteria[i] <- (total*rActino)
  total <- total - ADControl$Actinobacteria[i]
  
  ###Generate number for Proteobacteria
  rProteo <- rpearson(n=1, params=ProteoADCparam)
  ADControl$Proteobacteria[i] <- (total*rProteo)
  total <- total - ADControl$Proteobacteria[i]
  
  ###Generate number for Bacteridetes
  rBacter <- rpearson(n=1, params=BacterADCparam)
  ADControl$Bacteroidetes[i] <- (total*rBacter)
  
  ###Remainder is placed in Others
  ADControl$Others[i] <- total - ADControl$Bacteroidetes[i]
  
}

###0-100 AD Control extracted individual taxa dataset (percentage)
Firm1 <- as.numeric(control[1,])*100
Actino1 <- as.numeric(control[2,])*100
Proteo1 <- as.numeric(control[3,])*100
Bacter1 <- as.numeric(control[4,])*100

###Distribution parameters 0-100 (percentage)
FirmADCparam <- pearsonFitML(Firm1)
ActinoADCparam <- pearsonFitML(Actino1)
ProteoADCparam <- pearsonFitML(Proteo1)
BacterADCparam <- pearsonFitML(Bacter1)

###Histograme of each taxa distributions 0-100 (percentage) 
hist(rpearson(1000, params=FirmADCparam), breaks=100)
hist(rpearson(1000, params=ActinoADCparam), breaks=100)
hist(rpearson(1000, params=ProteoADCparam), breaks=100)
hist(rpearson(1000, params=BacterADCparam), breaks=100)

###Generates AD control dataset (percentage)###
numrow <- 100 ###Number of subjects/samples you want to generate
numcol <- 5 ###Number of taxa
ADControl <- matrix(nrow=numrow,ncol=numcol)
ADControl <- data.frame(ADControl)
colnames(ADControl) <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Others")
size <- dim(ADControl)[1]

for (i in 1:size){
  total <- 100
  ###Generate number for Firmicutes
  ADControl$Firmicutes[i] <- rpearson(n=1, params=FirmADCparam)
  total <- total - ADControl$Firmicutes[i]
  
  ###Generate number for Actinobacteria
  rActino <- rpearson(n=1, params=ActinoADCparam)
  ADControl$Actinobacteria[i] <- (total*rActino)/100
  total <- total - ADControl$Actinobacteria[i]
  
  ###Generate number for Proteobacteria
  rProteo <- rpearson(n=1, params=ProteoADCparam)
  ADControl$Proteobacteria[i] <- (total*rProteo)/100
  total <- total - ADControl$Proteobacteria[i]
  
  ###Generate number for Bacteridetes
  rBacter <- rpearson(n=1, params=BacterADCparam)
  ADControl$Bacteroidetes[i] <- (total*rBacter)/100
  
  ###Remainder is placed in Others
  ADControl$Others[i] <- total - ADControl$Bacteroidetes[i]
  
}


##########################

###AD Flares with No Treatment %remaining dataset
ADFlaresNT <- read.table("ADFlareNoTreatment.txt")

###0-1 AD Flare No Treatment extracted each individual taxa dataset
FirmFNT <- as.numeric(ADFlaresNT[1,])
ActinoFNT <- as.numeric(ADFlaresNT[2,])
ProteoFNT <- as.numeric(ADFlaresNT[3,])
BacterFNT <- as.numeric(ADFlaresNT[4,])

###Fitting distribution 0-1
unlist(pearsonFitML(FirmFNT))
unlist(pearsonFitML(ActinoFNT))
unlist(pearsonFitML(ProteoFNT))
unlist(pearsonFitML(BacterFNT))

###Fitted distribution taxa generation

FirmicutesFNT <- rpearsonI(n=1000, a=0.1166093, b=0.1248777, location=0.5898664, scale=0.3871436)
ActinobacteriaFNT <- rpearsonI(n=1000, a=0.81383910,  b=0.34911813, location=-0.02031676, scale=0.73061272)
ProteobacteriaFNT <- rpearsonII(n=1000, a=0.4902809, location=0.1733658, scale=0.6121342)
BacteroidetesFNT <- rpearsonI(n=1000, a=0.4347489, b=0.2513799, location=0.1557466, scale=0.8442535)

###Histogram of generated taxa 0-1
hist(FirmicutesFNT, 100)
hist(ActinobacteriaFNT, 100)
hist(ProteobacteriaFNT, 100)
hist(BacteroidetesFNT, 100)

###0-100 AD Flare No Treatment extracted individual taxa dataset
FirmFNT1 <- as.numeric(ADFlaresNT[1,])*100
ActinoFNT1 <- as.numeric(ADFlaresNT[2,])*100
ProteoFNT1 <- as.numeric(ADFlaresNT[3,])*100
BacterFNT1 <- as.numeric(ADFlaresNT[4,])*100

###Fitting distribution 0-100
unlist(pearsonFitML(FirmFNT1))
unlist(pearsonFitML(ActinoFNT1))
unlist(pearsonFitML(ProteoFNT1))
unlist(pearsonFitML(BacterFNT1))

###Fitted distribution dataset generation
FirmicutesFNT1 <- rpearsonII(n=1000, a=0.1835496, location=73.3078333, scale=24.3931667)
ActinobacteriaFNT1 <- rpearsonI(n=1000, a=0.813595, b=0.346429, location=12.491984, scale=58.537612)
ProteobacteriaFNT1 <- rpearsonI(n=1000, a=0.4793018,  b=0.6781994, location=22.6042533, scale=58.4019944)
BacteroidetesFNT1 <- rpearsonI(n=1000, a=0.4361804, b=0.2796059, location=29.6034631, scale=70.3965525)

###Histogram of generated taxa 0-100
hist(FirmicutesFNT1, 100)
hist(ActinobacteriaFNT1, 100)
hist(ProteobacteriaFNT1, 100)
hist(BacteroidetesFNT1, 100)
