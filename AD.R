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
colnames(ADControl) <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Others") ###Names of columns
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

###Distribution parameters 0-1 (proportions)
FirmFNTparam <- pearsonFitML(FirmFNT)
ActinoFNTparam <- pearsonFitML(ActinoFNT)
ProteoFNTparam <- pearsonFitML(ProteoFNT)
BacterFNTparam <- pearsonFitML(BacterFNT)


###Histogram of taxa distributions
hist(rpearson(n=1000, params=FirmFNTparam), breaks=100)
hist(rpearson(n=1000, params=ActinoFNTparam), breaks=100)
hist(rpearson(n=1000, params=ProteoFNTparam), breaks=100)
hist(rpearson(n=1000, params=BacterFNTparam), breaks=100)

###Generates AD Control dataset (proportions)###
numrow <- 100 ###Number of subjects/samples you want to generate
numcol <- 5 ###Number of taxa
ADFNT <- matrix(nrow=numrow,ncol=numcol)
ADFNT <- data.frame(ADFNT)
colnames(ADFNT) <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Others") ###Names of columns
size <- dim(ADFNT)[1]

for (i in 1:size){
  total <- 1
  ###Generate number for Firmicutes
  ADFNT$Firmicutes[i] <- rpearson(n=1, params=FirmFNTparam)
  total <- total - ADFNT$Firmicutes[i]
  
  ###Generate number for Actinobacteria
  rActino <- rpearson(n=1, params=ActinoFNTparam)
  ADFNT$Actinobacteria[i] <- (total*rActino)
  total <- total - ADFNT$Actinobacteria[i]
  
  ###Generate number for Proteobacteria
  rProteo <- rpearson(n=1, params=ProteoFNTparam)
  ADFNT$Proteobacteria[i] <- (total*rProteo)
  total <- total - ADFNT$Proteobacteria[i]
  
  ###Generate number for Bacteridetes
  rBacter <- rpearson(n=1, params=BacterFNTparam)
  ADFNT$Bacteroidetes[i] <- (total*rBacter)
  
  ###Remainder is placed in Others
  ADFNT$Others[i] <- total - ADFNT$Bacteroidetes[i]
  
}

###AD Flare No Treatment Loop
library("PearsonDS")
ADFlaresNT <- read.table("ADFlareNoTreatment.txt")

###0-100 AD Flare No Treatment extracted individual taxa dataset (percentage)
FirmFNT1 <- as.numeric(ADFlaresNT[1,])*100
ActinoFNT1 <- as.numeric(ADFlaresNT[2,])*100
ProteoFNT1 <- as.numeric(ADFlaresNT[3,])*100
BacterFNT1 <- as.numeric(ADFlaresNT[4,])*100

###Fitting distribution 0-100 (percentage)
FirmFNTparam <- pearsonFitML(FirmFNT1)
ActinoFNTparam <-pearsonFitML(ActinoFNT1)
ProteoFNTparam <- pearsonFitML(ProteoFNT1)
BacterFNTparam <- pearsonFitML(BacterFNT1)

###Histogram of taxa distributions
hist(rpearson(n=1000, params=FirmFNTparam), breaks=100)
hist(rpearson(n=1000, params=ActinoFNTparam), breaks=100)
hist(rpearson(n=1000, params=ProteoFNTparam), breaks=100)
hist(rpearson(n=1000, params=BacterFNTparam), breaks=100)

numrow=100
numcol=5
ADFlareNT <- matrix(nrow=numrow,ncol=numcol)
ADFlareNT <- data.frame(ADFlareNT)
colnames(ADFlareNT) <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Others")
size <- dim(ADFlareNT)[1]

for (i in 1:size){
  total <- 100
  ###Generate number for Firmicute and place into table 
  ADFlareNT$Firmicutes[i] <- rpearson(n=1, params=FirmFNTparam)
  total <- total - ADFlareNT$Firmicutes[i]
  
  ###Generate number for Actinobacteria and place into table
  rActinoFNT <- rpearson(n=1, params=ActinoFNTparam)
  ADFlareNT$Actinobacteria[i] <- (total*rActinoFNT)/100
  total <- total - ADFlareNT$Actinobacteria[i] 
  
  ###Generate number for Proteobacteria and place into table 
  rProteoFNT <- rpearson(n=1, params=ProteoFNTparam)
  ADFlareNT$Proteobacteria[i] <- (total*rProteoFNT)/100
  total <- total - ADFlareNT$Proteobacteria[i]
  
  ###Generate number for Bacteroidetes and place into table
  rBacterFNT <- rpearson(n=1, params=BacterFNTparam)
  ADFlareNT$Bacteroidetes[i] <- (total*rBacterFNT)/100
  
  ###Remainder is placed into Others
  ADFlareNT$Others[i] <- total - ADFlareNT$Bacteroidetes[i]
  
  
}

