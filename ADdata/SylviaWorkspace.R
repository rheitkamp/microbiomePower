set.seed(1234)
library(PearsonDS)
library(HMP)

(Control <- read.table("ADControls.txt"))
(Baseline <- read.table("ADBaseline.txt"))
(FlareNT <- read.table("ADFlareNT.txt"))
(FlareT <- read.table("ADFlareT.txt"))
(PostFlare <- read.table("ADPostFlare.txt"))

colSums(Control)
colSums(Baseline)
colSums(FlareNT)
colSums(FlareT)
colSums(PostFlare)

#######
(ControlSD <- matrix(apply(Control,1,sd),dimnames=list(rownames(Control), "SD")))
(BaselineSD <- matrix(apply(Baseline,1,sd), dimnames=list(rownames(Baseline), "SD")))
(FlareNTSD <- matrix(apply(FlareNT,1,sd), dimnames=list(rownames(FlareNT), "SD")))
(FlareTSD <- matrix(apply(FlareT,1,sd), dimnames=list(rownames(FlareT), "SD")))
(PostFlareSD <- matrix(apply(PostFlare,1,sd), dimnames=list(rownames(PostFlare), "SD")))

(ControlMean <- matrix(rowMeans(Control), dimnames=list(rownames(Control), "Mean")))
(BaselineMean <- matrix(rowMeans(Baseline), dimnames=list(rownames(Baseline), "Mean")))
(FlareNTMean <- matrix(rowMeans(FlareNT), dimnames=list(rownames(FlareNT), "Mean")))
(FlareTMean <- matrix(rowMeans(FlareT), dimname=list(rownames(FlareT), "Mean")))
(PostFlareMean <- matrix(rowMeans(PostFlare), dimnames=list(rownames(PostFlare), "Mean")))

(ControlMSD <- cbind(ControlMean, ControlSD))
(ControlMSD <- ControlMSD*100)
(BaselineMSD <- cbind(BaselineMean, BaselineSD))
(BaselineMSD <- BaselineMSD*100)
(FlareNTMSD <- cbind(FlareNTMean, FlareNTSD))
(FlareNTMSD <- FlareNTMSD*100)
(FlareTMSD <- cbind(FlareTMean, FlareTSD))
(FlareTMSD <- FlareTMSD*100)
(PostFlareMSD <- cbind(PostFlareMean, PostFlareSD))
(PostFlareMSD <- PostFlareMSD*100)

####Beta Dist before

#####
pMSD <- function(x){
  (nrow <- dim(x)[1])
  (ncol <- dim(x)[2])
  (y <- matrix(nrow=nrow,ncol=ncol, dimnames=list(rownames(x), colnames(x))))
  (total <- 100)
  (y[1,2] <- (x[1,2])/100)
  (y[1,1] <- (x[1,1])/100)
  for(i in 2:nrow){
    (y[i,2] <- x[i,2]/(total - x[i-1,1]))
    (total <- total - x[i-1,1])
    (y[i,1] <- x[i,1]/total)
  }
  return(y)
}

#####

#####
(nControl <- pMSD(ControlMSD))
(nBaseline <- pMSD(BaselineMSD))
(nFlareNT <- pMSD(FlareNTMSD))
(nFlareT <- pMSD(FlareTMSD))
(nPostFlare <- pMSD(PostFlareMSD))

####
pcontroltn <- urnorm(n=1000, mean=nControl[2,1], sd=nControl[2,2],lb=0,ub=1)
hist(pcontroltn,100)
sd(pcontroltn)
mean(pcontroltn)
####
numrow <- 25
numcol <- dim(ControlMSD)[1]
Cdata <- matrix(nrow=numrow,ncol=numcol)
colnames(Cdata) <- dimnames(ControlMSD)[[1]]
rownames(Cdata) <- rownames(Cdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(Cdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  Cdata[i,1] <- urnorm(n=1,mean=nControl[1,1],sd=nControl[1,2],lb=0,ub=1)
  total <- total - Cdata[i,1]
  
  ###Actinobacteria
  r <- urnorm(n=1, mean=nControl[2,1], sd=nControl[2,2], lb=0, ub=1)
  Cdata[i,2] <- total*r
  total <- total - Cdata[i,2]
  
  ###Proteobacteria
  r <- urnorm(n=1, mean=nControl[3,1], sd=nControl[3,2], lb=0, ub=1)
  Cdata[i,3] <- total*r
  total <- total - Cdata[i,3]
  
  ###Bacteriodetes
  r <- urnorm(n=1, mean=nControl[4,1], sd=nControl[4,2], lb=0, ub=1)
  Cdata[i,4] <- total*r
  total <- total - Cdata[i,4]
  
  ###Other
  Cdata[i,5] <- total
  
}


Cdata
colMeans(Cdata)
apply(Cdata,2,sd)
ControlMSD

#####
getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

########Parameters

(CparFirm <- getBetaParams(nControl[1,1], nControl[1,2]))
rpearson(n=1, params=CparFirm)
(CparActino <- getBetaParams(nControl[2,1], nControl[2,2]))
rpearson(n=1, params=CparActino)
(CparProteo <- getBetaParams(nControl[3,1], nControl[3,2]))
rpearson(n=1, params=CparProteo)
(CparBact <- getBetaParams(nControl[4,1], nControl[4,2]))
rpearson(n=1, params=CparBact)

(BparFirm <- getBetaParams(nBaseline[1,1], nBaseline[1,2]))
rpearson(n=1, params=BparFirm)
(BparActino <- getBetaParams(nBaseline[2,1], nBaseline[2,2]))
rpearson(n=1, params=BparActino)
(BparProteo <- getBetaParams(nBaseline[3,1], nBaseline[3,2]))
rpearson(n=1, params=BparProteo)
(BparBact <- getBetaParams(nBaseline[4,1], nBaseline[4,2]))
rpearson(n=1, params=BparBact)

(FNTparFirm <- getBetaParams(nFlareNT[1,1], nFlareNT[1,2]))
rpearson(n=1, params=FNTparFirm)
(FNTparActino <- getBetaParams(nFlareNT[2,1], nFlareNT[2,2]))
rpearson(n=1, params=FNTparActino)
(FNTparProteo <- getBetaParams(nFlareNT[3,1], nFlareNT[3,2]))
rpearson(n=1, params=FNTparProteo)
(FNTparBact <- getBetaParams(nFlareNT[4,1], nFlareNT[4,2]))
rpearson(n=1, params=FNTparBact)

(FTparFirm <- getBetaParams(nFlareT[1,1], nFlareT[1,2]))
rpearson(n=1, params=FTparFirm)
(FTparActino <- getBetaParams(nFlareT[2,1], nFlareT[2,2]))
rpearson(n=1, params=FTparActino)
(FTparProteo <- getBetaParams(nFlareT[3,1], nFlareT[3,2]))
rpearson(n=1, params=FTparProteo)
(FTparBact <- getBetaParams(nFlareT[4,1], nFlareT[4,2]))
rpearson(n=1, params=FTparBact)

(PFparFirm <- getBetaParams(nPostFlare[1,1], nPostFlare[1,2]))
rpearson(n=1, params=PFparFirm)
(PFparActino <- getBetaParams(nPostFlare[2,1], nPostFlare[2,2]))
rpearson(n=1, params=PFparActino)
(PFparProteo <- getBetaParams(nPostFlare[3,1], nPostFlare[3,2]))
rpearson(n=1, params=PFparProteo)
(PFparBact <- getBetaParams(nPostFlare[4,1], nPostFlare[4,2]))
rpearson(n=1, params=PFparBact)

####PostFlare
numrow <- 25
numcol <- dim(PostFlareMSD)[1]
PFdata <- matrix(nrow=numrow,ncol=numcol)
colnames(PFdata) <- dimnames(PostFlareMSD)[[1]]
rownames(PFdata) <- rownames(PFdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(PFdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  PFdata[i,1] <- rpearson(n=1, params=PFparFirm)
  total <- total - PFdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=PFparActino)
  PFdata[i,2] <- total*r
  total <- total - PFdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=PFparProteo)
  PFdata[i,3] <- total*r
  total <- total - PFdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=PFparBact)
  PFdata[i,4] <- total*r
  total <- total - PFdata[i,4]
  
  ###Other
  PFdata[i,5] <- total
  
}
colMeans(PFdata)
apply(PFdata,2,sd)
PostFlareMSD







####Control
numrow <- 25
numcol <- dim(ControlMSD)[1]
Cdata <- matrix(nrow=numrow,ncol=numcol)
colnames(Cdata) <- dimnames(ControlMSD)[[1]]
rownames(Cdata) <- rownames(Cdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(Cdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  Cdata[i,1] <- rpearson(n=1, params=CparFirm)
  total <- total - Cdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=CparActino)
  Cdata[i,2] <- total*r
  total <- total - Cdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=CparProteo)
  Cdata[i,3] <- total*r
  total <- total - Cdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=CparBact)
  Cdata[i,4] <- total*r
  total <- total - Cdata[i,4]
  
  ###Other
  Cdata[i,5] <- total
  
}
Cdata