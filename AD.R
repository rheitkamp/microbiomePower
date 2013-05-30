library("PearsonDS")

###AD Controls
control <- read.table("ADcontrolsredist.txt")

###0-1 dataset
Firm <- as.numeric(control[1,])
Actino <- as.numeric(control[2,])
Proteo <- as.numeric(control[3,])
Bacter <- as.numeric(control[4,])

###Fitting distribution 0-1
unlist(pearsonFitML(Firm))
unlist(pearsonFitML(Actino))
unlist(pearsonFitML(Proteo))
unlist(pearsonFitML(Bacter))


###fitted distribution dataset generation 0-1
Firmicutes <- rpearsonI(n=1000, a=0.9977008, b=2.1005555, location=0.1063188, scale=0.9611892)
Actinobacteria<- rpearsonI(n=1000, a=0.99874103, b=2.47718315, location=0.09764036, scale=1.03623409)
Proteobacteria <- rpearsonI(n=1000, a=1.3711653, b=1.4778528, location=0.4542903, scale=0.4134097)
Bacteroidetes <- rpearsonI(n=1000, a=2.2113205, b=0.6628346, location=-0.2730789, scale=1.2522455)

###Histograms 0-1
hist(Firmicutes,100)
hist(Actinobacteria,100)
hist(Proteobacteria,100)
hist(Bacteroidetes, 100)

###0-100 dataset
Firm1 <- as.numeric(control[1,])*100
Actino1 <- as.numeric(control[2,])*100
Proteo1 <- as.numeric(control[3,])*100
Bacter1 <- as.numeric(control[4,])*100

###Fitting distribution 0-100
unlist(pearsonFitML(Firm1))
unlist(pearsonFitML(Actino1))
unlist(pearsonFitML(Proteo1))
unlist(pearsonFitML(Bacter1))

###fitted distribution dataset generation 0-100
Firmicutes1 <- rpearsonI(n=1000, a=0.9970149, b=2.0062602, location=10.6318818, scale=89.3067879)
Actinobacteria1 <- rpearsonII(n=1000, a=0.840323, location=9.764036, scale=81.726210)
Proteobacteria1 <- rpearsonI(n=1000, a=1.371166,  b=1.477854, location=45.429029, scale=41.340968)
Bacteroidetes1 <- rpearsonI(n=1000, a=2.2098166, b=0.6798847, location=-12.4973479, scale=110.4140148)

###Histograms 0-100
hist(Firmicutes1, 100)
hist(Actinobacteria1, 100)
hist(Proteobacteria1, 100)
hist(Bacteroidetes1, 100)



##########################

###AD Flares with No Treatment

ADFlaresNT <- read.table("ADFlareNoTreatment.txt")

###0-1 AD Flare No Treatment
FirmFNT <- as.numeric(ADFlaresNT[1,])
ActinoFNT <- as.numeric(ADFlaresNT[2,])
ProteoFNT <- as.numeric(ADFlaresNT[3,])
BacterFNT <- as.numeric(ADFlaresNT[4,])

###Fitting distribution 0-1
unlist(pearsonFitML(FirmFNT))
unlist(pearsonFitML(ActinoFNT))
unlist(pearsonFitML(ProteoFNT))
unlist(pearsonFitML(BacterFNT))

###Fitted distribution dataset generation

FirmicutesFNT <- rpearsonI(n=1000, a=0.1166093, b=0.1248777, location=0.5898664, scale=0.3871436)
ActinobacteriaFNT <- rpearsonI(n=1000, a=0.81383910,  b=0.34911813, location=-0.02031676, scale=0.73061272)
ProteobacteriaFNT <- rpearsonII(n=1000, a=0.4902809, location=0.1733658, scale=0.6121342)
BacteroidetesFNT <- rpearsonI(n=1000, a=0.4347489, b=0.2513799, location=0.1557466, scale=0.8442535)

###Histogram of generated taxa 0-1
hist(FirmicutesFNT, 100)
hist(ActinobacteriaFNT, 100)
hist(ProteobacteriaFNT, 100)
hist(BacteroidetesFNT, 100)

###0-100 AD Flare No Treatment
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
