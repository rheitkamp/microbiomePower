###Creates UpperRep %remainder dataset
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

###Creates ADControl %remainder dataset
b <- read.table("ADControlsData.txt")

rowsize <- dim(b)[1]
colsize <- dim(b)[2]
dims <- dimnames(b)

PR <- matrix(nrow=rowsize, ncol=colsize, dimnames=dims)

for(i in 1:colsize){
  total <- apply(X=b, MARGIN=2, sum)[i]
  PR[1,i] <- b[1,i]
  total <- total - b[1,i]
  for(k in 2:rowsize){
    PR[k,i] <- b[k,i]/total
    total <- total - b[k,i]
  }
}

###Creates ADFNT %remainder dataset
b <- read.table("ADFNTdata.txt")

rowsize <- dim(b)[1]
colsize <- dim(b)[2]
dims <- dimnames(b)

PRFNT <- matrix(nrow=rowsize, ncol=colsize, dimnames=dims)

for(i in 1:colsize){
  total <- apply(X=b, MARGIN=2, sum)[i]
  PRFNT[1,i] <- b[1,i]
  total <- total - b[1,i]
  for(k in 2:rowsize){
    PRFNT[k,i] <- b[k,i]/total
    total <- total - b[k,i]
  }
}