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
UpperRep <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRep) <- dimnames(PRemainder)[2]


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


######################

Proteobacteria <- rpearsonI(n=1000, a=1.7133906, b=0.5593122, location=0, scale=1)
hist(Proteobacteria)
Firmicutes <- rpearsonI(n=1000, a=0.7207408, b=3.4647458, location=0, scale=1)
hist(Firmicutes)
Bacteroidetes <- rpearsonI(n=1000, a=0.1034732, b=2.9488343, location=0, scale=1)
hist(Bacteroidetes)
Actinobacteria <- rpearsonI(n=1000, a=0.2883812, b=15.0510462, location=0, scale=1)
hist(Actinobacteria)
Fusobacteria <- rpearsonI(n=1000, a=0.02791087, b=1.56699613, location=0, scale=1)
hist(Fusobacteria)
Cyanobacteria <- rpearsonI(n=1000, a=0.3728329, b=137.7134274, location=0, scale=1)
hist(Cyanobacteria)
OD1 <- rpearsonI(n=1000, a=0.605873, b=903.682256, location=0, scale=1)
hist(OD1)
TM7 <- rpearsonI(n=1000, a=0.02130855, b=47.33101610, location=0, scale=1)
hist(TM7)
DT <- rpearsonI(n=1000, a=0.2498125, b=1665.1668542, location=0, scale=1)
hist(DT)
Nitrospira <- rpearsonI(n=1000, a=0.1735524, b=3470.8750587, location=0, scale=1)
hist(Nitrospira)
Planctomycetes <- rpearsonI(n=1000, a=1.561484e-02, b=1.561469e+03, location=0, scale=1)
hist(Planctomycetes)
Chloroflexi <- rpearsonI(n=1000, a=3.99896e-02, b=3.99892e+03, location=0, scale=1)
hist(Chloroflexi)


numrow <- 100
numcol <- dim(PRemainder)[2]
data <- matrix(nrow=numrow, ncol=numcol)
names <- colnames(PRemainder)
