###Fall/Winter
a <- read.table("UpperRmean.txt")
c <- read.table("ADControlP.txt")
###Function need to fix

Premainder <- function(x) {
  rowsize <- dim(x)[1]
  colsize <- dim(x)[2]
  dims <- dimnames(x)
  y <- matrix(nrow=rowsize, ncol=colsize, dimnames=dims)
  for(i in 1:rowsize){
    total <- apply(X=x,MARGIN=1,sum)[i]
    y[i,1] <- x[i,1]/100
    total <- total - x[i,1]
    for(k in 2:colsize){
      y[1,k] <- x[1,k]/ total
      total <- total - x[1,k]
    }
  } 
  return(matrix(y,nrow=rowsize, ncol=colsize, dimnames=dims))
}
  
rowsize <- dim(c)[1]
colsize <- dim(c)[2]
dims <- dimnames(c)
y <- matrix(nrow=rowsize, ncol=colsize, dimnames=dims)
for(i in 1:rowsize){
  total <- apply(X=c,MARGIN=1,sum)[2]
  y[2,1] <- c[2,1]/100
  total <- total - c[2,1]
  for(k in 2:colsize){
    y[1,2] <- c[1,2]/ total
    total <- total - c[1,2]

##Mean = alpha/alpha+beta
##Variance = alpha*beta/(alpha+beta)^2 * (alpha + beta + 1)
##SD^2 = alpha*beta/(alpha+beta)^2 * (alpha + beta + 1)
## a= alpha b=beta
##using WolframAlpha to fine alpha and beta

Proteobacteria <- rpearsonI(n=1000, a=1.7133906, b=0.5593122, location=0, scale=1)
Firmicutes <- rpearsonI(n=1000, a=4.651839, b=2.005531, location=0, scale=1)
Bacteroidetes <- rpearsonI(n=1000, a=13.56189, b=16.13825, location=0, scale=1)
Actinobacteria <- rpearsonI(n=1000, a=102.2587, b=117.1624, location=0, scale=1)
Fusobacteria <- rpearsonI(n=1000, a=381.28608,  b=88.02261, location=0, scale=1)
Cyanobacteria <- rpearsonI(n=1000, a=7651.471, b=3797.397, location=0, scale=1)
OD1 <- rpearsonI(n=1000, a=169009.8, b=169009.8, location=0, scale=1)
TM7 <- rpearsonI(n=1000, a=15922.282, b=7784.227, location=0, scale=1)
DT <- rpearsonI(n=1000, a=1643500.4, b=766966.9, location=0, scale=1)
Nitrospira <- rpearsonI(n=1000, a=10123096, b=4049238, location=0, scale=1)
Planctomycetes <- rpearsonI(n=1000, a=19531249, b=19531249, location=0, scale=1)

###Loop
numrow <- 50
numcol <- dim(PRemainder)[2]
UpperRep <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRep) <- dimnames(PRemainder)[[2]]
rownames(UpperRep) <- rownames(UpperRep, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRep)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRep[i,1] <- rpearsonI(n=1, a=1.7133906, b=0.5593122, location=0, scale=1) 
  total <- total - UpperRep[i,1]
  
  ###Firmicutes
  r <- rpearsonI(n=1, a=4.651839, b=2.005531, location=0, scale=1)
  UpperRep[i,2] <- total*r
  total <- total - UpperRep[i,2]
  
  
  ###Bacteroidetes
  r <- rpearsonI(n=1, a=13.56189, b=16.13825, location=0, scale=1)
  UpperRep[i,3] <- total*r
  total <- total - UpperRep[i,3]
  
  ###Actinobacteria
  r <- rpearsonI(n=1, a=102.2587, b=117.1624, location=0, scale=1)
  UpperRep[i,4] <- total*r
  total <- total - UpperRep[i,4]
  
  ###Fusobacteria
  r <- rpearsonI(n=1, a=381.28608,  b=88.02261, location=0, scale=1)
  UpperRep[i,5] <- total*r
  total <- total - UpperRep[i,5]
  
  ###Cyanobacteria
  r <- rpearsonI(n=1, a=7651.471, b=3797.397, location=0, scale=1)
  UpperRep[i,6] <- total*r
  total <- total - UpperRep[i,6]
  
  ###OD1
  r <- rpearsonI(n=1, a=169009.8, b=169009.8, location=0, scale=1)
  UpperRep[i,7] <- total*r
  total <- total - UpperRep[i,7]
  
  ###TM7
  r <- rpearsonI(n=1, a=15922.282, b=7784.227, location=0, scale=1)
  UpperRep[i,8] <- total*r
  total <- total - UpperRep[i,8]
  
  ###Deinococcus-Thermus
  r <- rpearsonI(n=1, a=1643500.4, b=766966.9, location=0, scale=1)
  UpperRep[i,9] <- total*r
  total <- total - UpperRep[i,9]
 
  ###Nitrospira
  r <- rpearsonI(n=1, a=10123096, b=4049238, location=0, scale=1)
  UpperRep[i,10] <- total*r
  total <- total - UpperRep[i,10]
  
  ###Planctomycetes
  r <- rpearsonI(n=1, a=19531249, b=19531249, location=0, scale=1)
  UpperRep[i,11] <- total*r
  
  ###Chloroflexi
  UpperRep[i,12] <- total - UpperRep[i,11]
  
  ###BRC1
  UpperRep[1:numrow,13] <- 0
}


barchart(x=UpperRep,horizontal=FALSE, col=rainbow(13))


######Working#Backward######

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


numrow <- 50
numcol <- dim(PRemainder)[2]
data <- matrix(nrow=numrow, ncol=numcol)
colnames(data) <- dimnames(PRemainder)[[2]]
rownames(data) <- rownames(data, do.NULL= FALSE, prefix= "Sample")


size <- dim(data)[1]

for (i in 1:size){
  total <- 1
  
  ###BRC1
  data[i,13] <- 0
  total <- total - data[i,13]
  
  ###Chloroflexi
  data[i,12] <- rpearsonI(n=1, a=3.99896e-02, b=3.99892e+03, location=0, scale=1)
  total <- total - data[i,12]
  
  ###Planctomycetes
  data[i,11] <- rpearsonI(n=1, a=1.561484e-02, b=1.561469e+03, location=0, scale=1)
  total <- total - data[i,11]
  
  ###Nitrospira
  data[i,10] <- rpearsonI(n=1, a=0.1735524, b=3470.8750587, location=0, scale=1)
  total <- total - data[i,10]
  
  ###Deinococcus-Thermus
  data[i,9] <- rpearsonI(n=1, a=0.2498125, b=1665.1668542, location=0, scale=1)
  total <- total - data[i,9]
  
  ###TM7
  data[i,8] <- rpearsonI(n=1, a=0.02130855, b=47.33101610, location=0, scale=1)
  total <- total - data[i,8]
  
  ###OD1
  data[i,7] <- rpearsonI(n=1, a=0.605873, b=903.682256, location=0, scale=1)
  total <- total - data[i,7]
  
  ###Cyanobacteria
  data[i,6] <- rpearsonI(n=1, a=0.3728329, b=137.7134274, location=0, scale=1)
  total <- total - data[i,6]
  
  ###Fusobacteria
  data[i,5] <- rpearsonI(n=1, a=0.02791087, b=1.56699613, location=0, scale=1)
  total <- total - data[i, 5]
  
  ###Bacteroidetes
  data[i,3] <- rpearsonI(n=1, a=0.1034732, b=2.9488343, location=0, scale=1)
  total <- total - data[i,3]
  
  ###Actinobacteria
  data[i,4] <- rpearsonI(n=1, a=0.2883812, b=15.0510462, location=0, scale=1)
  total <- total - data[i,4]
  
  ###Firmicutes
  data[i,2] <- rpearsonI(n=1, a=0.7207408, b=3.4647458, location=0, scale=1)  
  total <- total - data[i,2]
  
  ###Proteobacteria
  data[i,1] <- total
  
}

barchart(data,horizontal=FALSE, col=rainbow(13))



###Spring
a <- read.table("UpperRepSpring.txt")

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

Proteobacteria <- rpearsonI(n=1000, a=1.749882, b=1.672536, location=0, scale=1)
Firmicutes <- rpearsonI(n=1000, a=1.483752, b=1.446118, location=0, scale=1)
Bacteroidetes <- rpearsonI(n=1000, a=3.831617, b=0.956813, location=0, scale=1)
Actinobacteria <- rpearsonI(n=1000, a=38.86014, b=12.19806, location=0, scale=1)
Fusobacteria <- rpearsonI(n=1000, a=23.80443, b=4.17792, location=0, scale=1)
Cyanobacteria <- rpearsonI(n=1000, a=77236.73, b=70371.24, location=0, scale=1)
OD1 <- rpearsonI(n=1000, a=473792.9, b=370794.5, location=0, scale=1)
TM7 <- rpearsonI(n=1000, a=64644.76, b=36538.34, location=0, scale=1)
DT <- rpearsonI(n=1000, a=7431283, b=4644552, location=0, scale=1)
Nitrospira <- rpearsonI(n=1000, a=12800000, b=51199999, location=0, scale=1)
Planctomycetes <- rpearsonI(n=1000, a=7396449, b=7396449, location=0, scale=1)
Chloroflexi <- rpearsonI(n=1000, a=5e+07, b=5e+07, location=0, scale=1)


###Loop
numrow <- 50
numcol <- dim(PRemainder)[2]
UpperRepSpringSpring <- matrix(nrow=numrow,ncol=numcol)
colnames(UpperRepSpringSpring) <- dimnames(PRemainder)[[2]]
rownames(UpperRepSpringSpring) <- rownames(UpperRepSpringSpring, do.NULL= FALSE, prefix= "Sample")


size <- dim(UpperRepSpring)[1]

for (i in 1:size){
  total <- 1
  ###Proteobacteria
  UpperRepSpring[i,1] <- rpearsonI(n=1, a=1.749882, b=1.672536, location=0, scale=1)  
  total <- total - UpperRepSpring[i,1]
  
  ###Firmicutes
  r <- rpearsonI(n=1, a=1.483752, b=1.446118, location=0, scale=1)
  UpperRepSpring[i,2] <- total*r
  total <- total - UpperRepSpring[i,2]
  
  
  ###Bacteroidetes
  r <- rpearsonI(n=1, a=3.831617, b=0.956813, location=0, scale=1)
  UpperRepSpring[i,3] <- total*r
  total <- total - UpperRepSpring[i,3]
  
  ###Actinobacteria
  r <- rpearsonI(n=1, a=38.86014, b=12.19806, location=0, scale=1)
  UpperRepSpring[i,4] <- total*r
  total <- total - UpperRepSpring[i,4]
  
  ###Fusobacteria
  r <- rpearsonI(n=1, a=23.80443, b=4.17792, location=0, scale=1)
  UpperRepSpring[i,5] <- total*r
  total <- total - UpperRepSpring[i,5]
  
  ###Cyanobacteria
  r <- rpearsonI(n=1, a=77236.73, b=70371.24, location=0, scale=1)
  UpperRepSpring[i,6] <- total*r
  total <- total - UpperRepSpring[i,6]
  
  ###OD1
  r <- rpearsonI(n=1, a=473792.9, b=370794.5, location=0, scale=1)
  UpperRepSpring[i,7] <- total*r
  total <- total - UpperRepSpring[i,7]
  
  ###TM7
  r <- rpearsonI(n=1, a=64644.76, b=36538.34, location=0, scale=1)
  UpperRepSpring[i,8] <- total*r
  total <- total - UpperRepSpring[i,8]
  
  ###Deinococcus-Thermus
  r <- DT <- rpearsonI(n=1, a=7431283, b=4644552, location=0, scale=1)
  UpperRepSpring[i,9] <- total*r
  total <- total - UpperRepSpring[i,9]
  
  ###Nitrospira
  r <- rpearsonI(n=1, a=12800000, b=51199999, location=0, scale=1)
  UpperRepSpring[i,10] <- total*r
  total <- total - UpperRepSpring[i,10]
  
  ###Planctomycetes
  r <- rpearsonI(n=1, a=7396449, b=7396449, location=0, scale=1)
  UpperRepSpring[i,11] <- total*r
  
  ###Chloroflexi
  r <- rpearsonI(n=1, a=5e+07, b=5e+07, location=0, scale=1)  
  UpperRepSpring[i,12] <- total*r
  
  ###BRC1
  UpperRepSpring[i,13] <- total - UpperRepSpring[i,12]
}


barchart(x=UpperRepSpring,horizontal=FALSE, col=rainbow(13))


######Working#Backward######

Proteobacteria <- rpearsonI(n=1000, a=1.749882, b=1.672536, location=0, scale=1)
hist(Proteobacteria)
Firmicutes <- rpearsonI(n=1000, a=0.4776035, b=1.4513283, location=0, scale=1)
hist(Firmicutes)
Bacteroidetes <- rpearsonI(n=1000, a=0.8961366, b=3.7446536, location=0, scale=1)
hist(Bacteroidetes)
Actinobacteria <- rpearsonI(n=1000, a=0.336087, b=8.804163, location=0, scale=1)
hist(Actinobacteria)
Fusobacteria <- rpearsonI(n=1000, a=0.01189997, b=1.20238249, location=0, scale=1)
hist(Fusobacteria)
Cyanobacteria <- rpearsonI(n=1000, a=0.4779586, b=530.5871302, location=0, scale=1)
hist(Cyanobacteria)
OD1 <- rpearsonI(n=1000, a=0.7248578, b=1575.0529200, location=0, scale=1)
hist(OD1)
TM7 <- rpearsonI(n=1000, a=0.0229654, b=99.8266026, location=0, scale=1)
hist(TM7)
DT <- rpearsonI(n=1000, a=0.3264245, b=4079.9796980, location=0, scale=1)
hist(DT)
Nitrospira <- rpearsonI(n=1000, a=3.99896e-02, b=3.99892e+0, location=0, scale=1)
hist(Nitrospira)
Planctomycetes <- rpearsonI(n=1000, a=2.364817e-02, b=1.182385e+03, location=0, scale=1)
hist(Planctomycetes)
Chloroflexi <- rpearsonI(n=1000, a=3.99896e-02, b=3.99892e+03, location=0, scale=1)
hist(Chloroflexi)
BRC1 <- rpearsonI(n=1000, a=1.233556e-02, b=1.233543e+03, location=0, scale=1)


numrow <- 50
numcol <- dim(PRemainder)[2]
data <- matrix(nrow=numrow, ncol=numcol)
colnames(data) <- dimnames(PRemainder)[[2]]
rownames(data) <- rownames(data, do.NULL= FALSE, prefix= "Sample")


size <- dim(data)[1]

for (i in 1:size){
  total <- 1
  
  ###BRC1
  data[i,13] <- rpearsonI(n=1, a=1.233556e-02, b=1.233543e+03, location=0, scale=1)
  total <- total - data[i,13]
  
  ###Chloroflexi
  data[i,12] <- rpearsonI(n=1, a=3.99896e-02, b=3.99892e+03, location=0, scale=1)
  total <- total - data[i,12]
  
  ###Planctomycetes
  data[i,11] <- rpearsonI(n=1, a=2.364817e-02, b=1.182385e+03, location=0, scale=1)
  total <- total - data[i,11]
  
  ###Nitrospira
  data[i,10] <- rpearsonI(n=1, a=3.99896e-02, b=3.99892e+0, location=0, scale=1)
  total <- total - data[i,10]
  
  ###Deinococcus-Thermus
  data[i,9] <- rpearsonI(n=1, a=0.3264245, b=4079.9796980, location=0, scale=1)
  total <- total - data[i,9]
  
  ###TM7
  data[i,8] <- rpearsonI(n=1, a=0.0229654, b=99.8266026, location=0, scale=1)
  total <- total - data[i,8]
  
  ###OD1
  data[i,7] <- rpearsonI(n=1, a=0.7248578, b=1575.0529200, location=0, scale=1)
  total <- total - data[i,7]
  
  ###Cyanobacteria
  data[i,6] <- rpearsonI(n=1, a=0.4779586, b=530.5871302, location=0, scale=1)
  total <- total - data[i,6]
  
  ###Fusobacteria
  data[i,5] <- rpearsonI(n=1, a=0.01189997, b=1.20238249, location=0, scale=1)
  total <- total - data[i,5]
  
  ###Actinobacteria
  data[i,4] <- rpearsonI(n=1, a=0.336087, b=8.804163, location=0, scale=1)
  total <- total - data[i,4]
  
  ###Bacteroidetes
  data[i,3] <- rpearsonI(n=1, a=0.8961366, b=3.7446536, location=0, scale=1)
  total <- total - data[i,3]
  
  
  ###Firmicutes
  data[i,2] <- rpearsonI(n=1, a=0.4776035, b=1.4513283, location=0, scale=1) 
  total <- total - data[i,2]
  
  ###Proteobacteria
  data[i,1] <- total
  
}

barchart(data,horizontal=FALSE, col=rainbow(13))