a <- read.table("UpperRmean.txt")

rowsize <- dim(a)[1]
colsize <- dim(a)[2]
dims <- dimnames(a)

PRemainder <- matrix(nrow=rowsize, ncol=colsize, dimnames=dims)


for(i in 1:rowsize){
  total <- apply(X=a,MARGIN=1,sum)
  PRemainder[i,1] <- a[i,1]/100
  total <- total - a[i,1]
  for(k in 2:colsize){
    PRemainder[1,k] <- a[1,k]/ total
    total <- total - a[1,k]
  }
} 


##Mean = alpha/alpha+beta
##Variance = alpha*beta/(alpha+beta)^2 * (alpha + beta + 1)
##SD^2 = alpha*beta/(alpha+beta)^2 * (alpha + beta + 1)
## a= alpha b=beta
##using WolframAlpha to fine alpha and beta

Proteobacteria <- rpearsonI(n=1000, a=1.71339, b=0.559312, location=0, scale=1)
Firmicutes <- rpearsonI(n=1000, a=4.65184, b=2.00553, location=0, scale=1)
Bacteroidetes <- rpearsonI(n=1000, a=13.5619, b=16.1382, location=0, scale=1)
Actinobacteria <- rpearsonI(n=1000, a=102.259, b=117.162, location=0, scale=1)
Fusobacteria <- rpearsonI(n=1000, a=17.8716, b=4.12579, location=0, scale=1)
Cyanobacteria <- rpearsonI(n=1000, a=7651.47, b=3797.4, location=0, scale=1)
OD1 <- rpearsonI(n=1000, a=169010, b=169010, location=0, scale=1)
TM7 <- rpearsonI(n=1000, a=15922.3, b=7784.23, location=0, scale=1)
DT <- rpearsonI(n=1000, a=1.6435e6, b=766967, location=0, scale=1)
Nitrospira <- rpearsonI(n=1000, a=1.01231e7, b=4.04924e6, location=0, scale=1)
Planctomycetes <- rpearsonI(n=1000, a=1.95312e7, b=1.95312e7, location=0, scale=1)

###Loop
numrow <- 100
numcol <- dim(PRemainder)[2]
names <- dimnames(PRemainder)[2]
UpperRep <- matrix(nrow=numrow,ncol=numcol, dimnames=list(c(1:numrow),
                                                          c(names)))

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