library(PearsonDS)
library(HMP)

set.seed(123)

###AD Control creates the %remaining dataset
b <- read.table("ADControldata.csv")

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

###Params
F <- pearsonFitML(as.numeric(PR[1,]))
A <- pearsonFitML(as.numeric(PR[2,]))
P <- pearsonFitML(as.numeric(PR[3,]))
B <- pearsonFitML(as.numeric(PR[4,]))

###AD Control dataset
colsize <- dim(PR)[1]
rowsize <- 20
ADControl <- matrix(nrow=rowsize, ncol=colsize)
colnames(ADControl) <- dimnames(PR)[[1]]
rownames(ADControl) <- rownames(ADControl, do.NULL= FALSE, prefix= "Sample")

for (i in 1:rowsize){
  total <- 1
  ADControl[i,1] <- rpearson(n=1, params=F)
  total <- total - ADControl[i,1]
  R <- rpearson(n=1, params=A)
  ADControl[i,2] <- total * R
  total <- total - ADControl[i,2]
  R <- rpearson(n=1, params=P)
  ADControl[i,3] <- total * R
  total <- total - ADControl[i,3]
  R <- rpearson(n=1, params=B)
  ADControl[i,4] <- total * R
  ADControl[i,5] <- total - ADControl[i,4]
}

###ADFNT creates %remaining dataset
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

###Params
FFNT <- pearsonFitML(as.numeric(PRFNT[1,]))
AFNT <- pearsonFitML(as.numeric(PRFNT[2,]))
PFNT <- pearsonFitML(as.numeric(PRFNT[3,]))
BFNT <- pearsonFitML(as.numeric(PRFNT[4,]))

###ADFNT dataset
colsize <- dim(PRFNT)[1]
rowsize <- 20
ADFNT <- matrix(nrow=rowsize, ncol=colsize)
colnames(ADFNT) <- dimnames(PRFNT)[[1]]
rownames(ADFNT) <- rownames(ADFNT, do.NULL= FALSE, prefix= "Sample")

for (i in 1:rowsize){
  total <- 1
  ADFNT[i,1] <- rpearson(n=1, params=FFNT)
  total <- total - ADFNT[i,1]
  R <- rpearson(n=1, params=AFNT)
  ADFNT[i,2] <- total * R
  total <- total - ADFNT[i,2]
  R <- rpearson(n=1, params=PFNT)
  ADFNT[i,3] <- total * R
  total <- total - ADFNT[i,3]
  R <- rpearson(n=1, params=BFNT)
  ADFNT[i,4] <- total * R
  ADFNT[i,5] <- total - ADFNT[i,4]
}

#####
ADControl500 <- ADControl*500
ADFNT500 <- ADFNT*500

fit.ADControl500 <- DM.MoM(ADControl500)
fit.ADFNT500 <- DM.MoM(ADFNT500)

MC <- 1000

Nrs1 <- rep(500,5)
Nrs2 <- rep(500,5)

group.Nrs <- list(Nrs1, Nrs2)

group.alphap <- rbind(fit.ADControl500$gamma, fit.ADFNT500$gamma)

P <- MC.Xdc.statistics(group.Nrs, MC, group.alphap, 2, "ha", 0.05)

group.pi <- rbind(fit.ADControl500$pi, fit.ADFNT500$pi)
group.theta <- c(fit.ADControl500$theta, fit.ADFNT500$theta)
P <- MC.Xmc.statistics(group.Nrs, MC, fit.ADControl500$pi, group.pi, group.theta, "ha")

P