a <- read.table("UpperRmean.txt")


rowsize <- dim(a)[1]
colsize <- dim(a)[2]
dims <- dimnames(a)

PRemainder <- matrix(nrow=rowsize, ncol=colsize, dimnames=dims)

PRemainder <- data.frame(pmean)

total <- 100
k <- 1

PRemainder[k,1] <- a[k,1]/100

total <- total - a[k,1]

k <- k + 1

