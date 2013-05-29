control <- read.table("ADcontrolsredist.txt")
control <- data.frame(controldataset)
sapply(controldata,mode)



##0-1 dataset
controlFirmicutes <- array(c(0.359244469, 0.671044383, 0.276331639, 0.269166758, 0.150069156, 0.333149184, 0.443099196, 0.516356896, 0.357908236, 0.436807996, 0.382578589, 0.285195937, 0.125090059, 0.106318818, 0.875177396, 0.691970877, 0.645803699, 0.554053544, 0.26962963, 0.474536405, 0.579606773, 0.554370957))
controlActinobacteria <- array(c(0.138580221, 0.525894166, 0.250027094, 0.432774087, 0.097640358, 0.148109752, 0.429304816,	0.527600906, 0.285880301, 0.419387503, 0.4735928, 0.452791878, 0.88929961, 0.887595165, 0.357144544, 0.448015498, 0.144578314, 0.421637846,	0.172413792,	0.428535216, 0.376990496,	0.434151517))



##0-100 dataset
controlFirmicutes <- array(c(0.359244469, 0.671044383, 0.276331639, 0.269166758, 0.150069156, 0.333149184, 0.443099196, 0.516356896, 0.357908236, 0.436807996, 0.382578589, 0.285195937, 0.125090059, 0.106318818, 0.875177396, 0.691970877, 0.645803699, 0.554053544, 0.26962963, 0.474536405, 0.579606773, 0.554370957))*100
controlActinobacteria <- array(c(0.138580221, 0.525894166, 0.250027094, 0.432774087, 0.097640358, 0.148109752, 0.429304816,  0.527600906, 0.285880301, 0.419387503, 0.4735928, 0.452791878, 0.88929961, 0.887595165, 0.357144544, 0.448015498, 0.144578314, 0.421637846,	0.172413792,	0.428535216, 0.376990496,	0.434151517))*100


##fitted distribution 0-1
Firmicutes <- rpearsonI(n=100, 0.9977008, 2.1005555, 0.1063188, 0.9611892)


##fitted distribution 0-100
Firmicutes <- rpearsonI(n=1000, 0.9922127,  2.0054740, 10.6318818, 89.3081162)



##########################

library("PearsonDS")
controldataset <- read.table("ADcontrolsredist.txt")
controldata <- data.frame(controldataset)
sapply(controldata,mode)


controldata <- array(controldataset)
controldata

Firm <- as.numeric(controldata[1,])
Firm
Actino <- as.numeric(controldata[2,])
Proteo <- as.numeric(controldata[3,])
Bacter <- as.numeric(controldata[4,])
hist(Bacter,10)
Other <- as.numeric(controldata[5,])


library("PearsonDS")

unlist(pearsonFitML(Firm))
unlist(pearsonFitML(Actino))
unlist(pearsonFitML(Proteo))
unlist(pearsonFitML(Bacter))
unlist(pearsonFitML(Other))

Firmicutes <- rpearsonI(n=1000, a=0.9977008, b=2.1005555, location=0.1063188, scale=0.9611892)
hist(Firmicutes,100)
Actinobacteria<- rpearsonI(n=1000, a=0.99874103, b=2.47718315, location=0.09764036, scale=1.03623409)
hist(Actinobacteria,100)
Proteobacteria <- rpearsonI(n=1000, a=1.3711653, b=1.4778528, location=0.4542903, scale=0.4134097)
hist(Proteobacteria,100)
Bacteroidetes <- rpearsonI(n=1000, a=2.2113205, b=0.6628346, location=-0.2730789, scale=1.2522455)
hist(Bacteroidetes, 100)
Others <- rpearsonIV(n=1000, m=3.207085e+00, nu=1.255142e+00, location=1.000000e+00, scale=1.673377e-07)
hist(Others,100)



library("HMP")
data(saliva)