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
numrow <- 50
numcol <- dim(PRemainder)[2]
names <- dimnames(PRemainder)[2]
UpperRep <- matrix(nrow=numrow,ncol=numcol)

size <- dim(UpperRep)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRep[i,1] <- rpearsonI(n=1, a=1.71339, b=0.559312, location=0, scale=1)  
  total <- total - UpperRep[i,1]
  
  ###Firmicutes
  r <- rpearsonI(n=1, a=4.65184, b=2.00553, location=0, scale=1)
  UpperRep[i,2] <- total*r
  total <- total - UpperRep[i,2]
  
  
  ###Bacteroidetes
  r <- rpearsonI(n=1, a=13.5619, b=16.1382, location=0, scale=1)
  UpperRep[i,3] <- total*r
  total <- total - UpperRep[i,3]
  
  ###Actinobacteria
  r <- rpearsonI(n=1, a=102.259, b=117.162, location=0, scale=1)
  UpperRep[i,4] <- total*r
  total <- total - UpperRep[i,4]
  
  ###Fusobacteria
  r <- rpearsonI(n=1, a=17.8716, b=4.12579, location=0, scale=1)
  UpperRep[i,5] <- total*r
  total <- total - UpperRep[i,5]
  
  ###Cyanobacteria
  r <- rpearsonI(n=1, a=7651.47, b=3797.4, location=0, scale=1)
  UpperRep[i,6] <- total*r
  total <- total - UpperRep[i,6]
  
  ###OD1
  r <- rpearsonI(n=1, a=169010, b=169010, location=0, scale=1)
  UpperRep[i,7] <- total*r
  total <- total - UpperRep[i,7]
  
  ###TM7
  r <- rpearsonI(n=1, a=15922.3, b=7784.23, location=0, scale=1)
  UpperRep[i,8] <- total*r
  total <- total - UpperRep[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearsonI(n=1, a=1.6435e6, b=766967, location=0, scale=1)
  UpperRep[i,9] <- total*r
  total <- total - UpperRep[i,9]
 
  ###Nitrospira
  r <- rpearsonI(n=1, a=1.01231e7, b=4.04924e6, location=0, scale=1)
  UpperRep[i,10] <- total*r
  total <- total - UpperRep[i,10]
  
  ###Planctomycetes
  r <- rpearsonI(n=1, a=1.95312e7, b=1.95312e7, location=0, scale=1)
  UpperRep[i,11] <- total*r
  
  ###Chloroflexi
  UpperRep[i,12] <- total - UpperRep[i,11]
  
  ###BRC1
  UpperRep[1:numrow,13] <- 0
}

barchart(x=UpperRep,horizontal=FALSE)