###AD Control Loop

library("PearsonDS")

###AD Controls %remaining dataset
control <- read.table("ADcontrolsredist.txt")

###0-100 AD Control extracted individual taxa dataset
Firm1 <- as.numeric(control[1,])*100
Actino1 <- as.numeric(control[2,])*100
Proteo1 <- as.numeric(control[3,])*100
Bacter1 <- as.numeric(control[4,])*100

###Fitting distribution 0-100
FirmADCparam <- pearsonFitML(Firm1)
ActinoADCparam <- pearsonFitML(Actino1)
ProteoADCparam <- pearsonFitML(Proteo1)
BacterADCparam <- pearsonFitML(Bacter1)

###Loop
numrow <- 100
numcol <- 5
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

