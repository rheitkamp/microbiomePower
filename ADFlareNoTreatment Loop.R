###AD Flare No Treatment Loop
library("PearsonDS")
ADFlaresNT <- read.table("ADFlareNoTreatment.txt")

###0-100 AD Flare No Treatment extracted individual taxa dataset
FirmFNT1 <- as.numeric(ADFlaresNT[1,])*100
ActinoFNT1 <- as.numeric(ADFlaresNT[2,])*100
ProteoFNT1 <- as.numeric(ADFlaresNT[3,])*100
BacterFNT1 <- as.numeric(ADFlaresNT[4,])*100

###Fitting distribution 0-100
FirmFNTparam <- pearsonFitML(FirmFNT1)
ActinoFNTparam <-pearsonFitML(ActinoFNT1)
ProteoFNTparam <- pearsonFitML(ProteoFNT1)
BacterFNTparam <- pearsonFitML(BacterFNT1)

###Fitted distribution dataset generation
FirmicutesFNT <- rpearson(n=1000, params=FirmFNTparam)
ActinobacteriaFNT <- rpearson(n=1000, params=ActinoFNTparam)
ProteobacteriaFNT <- rpearson(n=1000, params=ProteoFNTparam)
BacteroidetesFNT <- rpearson(n=1000, params=BacterFNTparam)

numrow=1000
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


