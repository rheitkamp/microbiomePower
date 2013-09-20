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

######
(a <- matrix(Control[1,], dimnames=list(colnames(Control[1,]), rownames(Control[1,]))))
(b <- matrix(Control[2,], dimnames=list(colnames(Control[2,]), rownames(Control[2,]))))
(c <- matrix(Control[3,], dimnames=list(colnames(Control[3,]), rownames(Control[3,]))))
(d <- matrix(Control[4,], dimnames=list(colnames(Control[4,]), rownames(Control[4,]))))
(e <- matrix(Control[5,], dimnames=list(colnames(Control[5,]), rownames(Control[5,]))))

(control <- cbind(a,b,c,d,e))
mode(control) <- "numeric"

(a <- matrix(Baseline[1,], dimnames=list(colnames(Baseline[1,]), rownames(Baseline[1,]))))
(b <- matrix(Baseline[2,], dimnames=list(colnames(Baseline[2,]), rownames(Baseline[2,]))))
(c <- matrix(Baseline[3,], dimnames=list(colnames(Baseline[3,]), rownames(Baseline[3,]))))
(d <- matrix(Baseline[4,], dimnames=list(colnames(Baseline[4,]), rownames(Baseline[4,]))))
(e <- matrix(Baseline[5,], dimnames=list(colnames(Baseline[5,]), rownames(Baseline[5,]))))

(baseline <- cbind(a,b,c,d,e))
mode(baseline) <- "numeric"

(a <- matrix(FlareNT[1,], dimnames=list(colnames(FlareNT[1,]), rownames(FlareNT[1,]))))
(b <- matrix(FlareNT[2,], dimnames=list(colnames(FlareNT[2,]), rownames(FlareNT[2,]))))
(c <- matrix(FlareNT[3,], dimnames=list(colnames(FlareNT[3,]), rownames(FlareNT[3,]))))
(d <- matrix(FlareNT[4,], dimnames=list(colnames(FlareNT[4,]), rownames(FlareNT[4,]))))
(e <- matrix(FlareNT[5,], dimnames=list(colnames(FlareNT[5,]), rownames(FlareNT[5,]))))

(flarent <- cbind(a,b,c,d,e))
mode(flarent) <- "numeric"

(a <- matrix(FlareT[1,], dimnames=list(colnames(FlareT[1,]), rownames(FlareT[1,]))))
(b <- matrix(FlareT[2,], dimnames=list(colnames(FlareT[2,]), rownames(FlareT[2,]))))
(c <- matrix(FlareT[3,], dimnames=list(colnames(FlareT[3,]), rownames(FlareT[3,]))))
(d <- matrix(FlareT[4,], dimnames=list(colnames(FlareT[4,]), rownames(FlareT[4,]))))
(e <- matrix(FlareT[5,], dimnames=list(colnames(FlareT[5,]), rownames(FlareT[5,]))))

(flaret <- cbind(a,b,c,d,e))
mode(flaret) <- "numeric"

(a <- matrix(PostFlare[1,], dimnames=list(colnames(PostFlare[1,]), rownames(PostFlare[1,]))))
(b <- matrix(PostFlare[2,], dimnames=list(colnames(PostFlare[2,]), rownames(PostFlare[2,]))))
(c <- matrix(PostFlare[3,], dimnames=list(colnames(PostFlare[3,]), rownames(PostFlare[3,]))))
(d <- matrix(PostFlare[4,], dimnames=list(colnames(PostFlare[4,]), rownames(PostFlare[4,]))))
(e <- matrix(PostFlare[5,], dimnames=list(colnames(PostFlare[5,]), rownames(PostFlare[5,]))))

(postflare <- cbind(a,b,c,d,e))
mode(postflare) <- "numeric"

#######
(controlmean <- colMeans(control))
(baselinemean <- colMeans(baseline))
(flarentmean <- colMeans(flarent))
(flaretmean <- colMeans(flaret))
(postflaremean <- colMeans(postflare))

(controlSD <- apply(control,2,sd))
(baselineSD <- apply(baseline,2,sd))
(flarentSD <- apply(flarent,2,sd))
(flaretSD <- apply(flaret,2,sd))
(postflareSD <- apply(postflare,2,sd))

#######
(ControlSD <- matrix(apply(Control,1,sd),dimnames=list(rownames(Control), "SD")))
(BaselineSD <- matrix(apply(Baseline,1,sd), dimnames=list(rownames(Baseline), "SD")))
(FlareNTSD <- matrix(apply(FlareNT,1,sd), dimnames=list(rownames(FlareNT), "SD")))
(FlareTSD <- matrix(apply(FlareT,1,sd), dimnames=list(rownames(FlareT), "SD")))
(PostFlareSD <- matrix(apply(PostFlare,1,sd), dimnames=list(rownames(PostFlare), "SD")))

(ControlMean <- matrix(rowMeans(Control), nrow=1, ncol=5, dimnames=list("Mean", rownames(Control))))
(BaselineMean <- matrix(rowMeans(Baseline), nrow=1, ncol=5, dimnames=list("Mean", rownames(Baseline))))
(FlareNTMean <- matrix(rowMeans(FlareNT), nrow=1, ncol=5, dimnames=list("Mean", rownames(FlareNT))))
(FlareTMean <- matrix(rowMeans(FlareT), nrow=1, ncol=5, dimname=list("Mean", rownames(FlareT))))
(PostFlareMean <- matrix(rowMeans(PostFlare), nrow=1, ncol=5, dimnames=list("Mean", rownames(PostFlare))))

#######
Premainder <- function(x) {
  rowsize <- dim(x)[1]
  colsize <- dim(x)[2]
  dims <- dimnames(x)[[2]]
  y <- matrix(nrow=rowsize, ncol=colsize)
  for(i in 1:rowsize){
    total <- apply(X=x,MARGIN=1,sum)[i]
    y[i,1] <- x[i,1]/100
    total <- total - x[i,1]
    for(k in 2:colsize){
      y[i,k] <- x[i,k]/ total
      total <- total - x[i,k]
    }
  } 
  return(matrix(y,nrow=colsize, ncol=rowsize, dimnames=list(dims, "Mean")))
}
######
Premainder <- function(x) {
  rowsize <- dim(x)[1]
  colsize <- dim(x)[2]
  y <- matrix(nrow=rowsize, ncol=colsize, dimnames=list(rownames(x),colnames(x)))
  for(i in 1:rowsize){
    total <- apply(X=x,MARGIN=1,sum)[i]
    y[i,1] <- x[i,1]/100
    total <- total - x[i,1]
    for(k in 2:colsize){
      y[i,k] <- x[i,k]/ total
      total <- total - x[i,k]
    }
  } 
  return(y)
}

########
(pcontrol <- Premainder(control*100))
(pbaseline <- Premainder(baseline*100))
(pflarent <- Premainder(flarent*100))
(pflaret <- Premainder(flaret*100))
(ppostflare <- Premainder(postflare*100))

(pcontrolmean <- colMeans(pcontrol))
(pbaselinemean <- colMeans(pbaseline))
(pflarentmean <- colMeans(pflarent))
(pflaretmean <- colMeans(pflaret))
(ppostflaremean <- colMeans(ppostflare))

rbind(pcontrolmean, pbaselinemean, pflarentmean, pflaretmean, ppostflaremean)

(pcontrolSD <- apply(pcontrol,2,sd))
(pbaselineSD <- apply(pbaseline,2,sd))
(pflarentSD <- apply(flarent,2,sd))
(pflaretSD <- apply(flaret,2,sd))
(ppostflareSD <- apply(ppostflare,2,sd))

rbind(pcontrolSD, pbaselineSD, pflarentSD, pflaretSD, ppostflareSD)

########
(ControlMP <- Premainder(ControlMean*100))
(BaselineMP <- Premainder(BaselineMean*100))
(FlareNTMP <- Premainder(FlareNTMean*100))
(FlareTMP <- Premainder(FlareTMean*100))
(PostFlareMP <- Premainder(PostFlareMean*100))

(ControlMSD <- cbind(ControlMP, ControlSD))
(BaselineMSD <- cbind(BaselineMP, BaselineSD))
(FlareNTMSD <- cbind(FlareNTMP, FlareNTSD))
(FlareTMSD <- cbind(FlareTMP, FlareTSD))
(PostFlareMSD <- cbind(PostFlareMP, PostFlareSD))

#######
getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

########
CparFirm <- getBetaParams(ControlMSD[1,1], ControlMSD[1,2])
CparActino <- getBetaParams(ControlMSD[2,1], ControlMSD[2,2])
CparProteo <- getBetaParams(ControlMSD[3,1], ControlMSD[3,2])
CparBact <- getBetaParams(ControlMSD[4,1], ControlMSD[4,2])

BparFirm <- getBetaParams(BaselineMSD[1,1], BaselineMSD[1,2])
BparActino <- getBetaParams(BaselineMSD[2,1], BaselineMSD[2,2])
BparProteo <- getBetaParams(BaselineMSD[3,1], BaselineMSD[3,2])
BparBact <- getBetaParams(BaselineMSD[4,1], BaselineMSD[4,2])

FNTparFirm <- getBetaParams(FlareNTMSD[1,1], FlareNTMSD[1,2]) 
FNTparActino <- getBetaParams(FlareNTMSD[2,1], FlareNTMSD[2,2])
FNTparProteo <- getBetaParams(FlareNTMSD[3,1], FlareNTMSD[3,2])
FNTparBact <- getBetaParams(FlareNTMSD[4,1], FlareNTMSD[4,2])

FTparFirm <- getBetaParams(FlareTMSD[1,1], FlareTMSD[1,2])
FTparActino <- getBetaParams(FlareTMSD[2,1], FlareTMSD[2,2])
FTparProteo <- getBetaParams(FlareTMSD[3,1], FlareTMSD[3,2])
FTparBact <- getBetaParams(FlareTMSD[4,1], FlareTMSD[4,2])

PFparFirm <- getBetaParams(PostFlareMSD[1,1], PostFlareMSD[1,2])
PFparActino <- getBetaParams(PostFlareMSD[2,1], PostFlareMSD[2,2])
PFparProteo <- getBetaParams(PostFlareMSD[3,1], PostFlareMSD[3,2])
PFparBact <- getBetaParams(PostFlareMSD[4,1], PostFlareMSD[4,2])

###Control###
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

apply(Cdata,2,sd)
apply(Control,1,sd)

colMeans(Cdata)
rowMeans(Control)

###Baseline###

numrow <- 25
numcol <- dim(BaselineMSD)[1]
Bdata <- matrix(nrow=numrow,ncol=numcol)
colnames(Bdata) <- dimnames(BaselineMSD)[[1]]
rownames(Bdata) <- rownames(Bdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(Bdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  Bdata[i,1] <- rpearson(n=1, params=BparFirm)
  total <- total - Bdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=BparActino)
  Bdata[i,2] <- total*r
  total <- total - Bdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=BparProteo)
  Bdata[i,3] <- total*r
  total <- total - Bdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=BparBact)
  Bdata[i,4] <- total*r
  total <- total - Bdata[i,4]
  
  ###Other
  Bdata[i,5] <- total
  
}

apply(Bdata,2,sd)
apply(Baseline,1,sd)

###FlareNT###

numrow <- 25
numcol <- dim(FlareNTMSD)[1]
FNTdata <- matrix(nrow=numrow,ncol=numcol)
colnames(FNTdata) <- dimnames(FlareNTMSD)[[1]]
rownames(FNTdata) <- rownames(FNTdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(FNTdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  FNTdata[i,1] <- rpearson(n=1, params=FNTparFirm)
  total <- total - FNTdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=FNTparActino)
  FNTdata[i,2] <- total*r
  total <- total - FNTdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=FNTparProteo)
  FNTdata[i,3] <- total*r
  total <- total - FNTdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=FNTparBact)
  FNTdata[i,4] <- total*r
  total <- total - FNTdata[i,4]
  
  ###Other
  FNTdata[i,5] <- total
  
}

apply(FNTdata,2,sd)
apply(FlareNT,1,sd)

###FlareT###

numrow <- 25
numcol <- dim(FlareTMSD)[1]
FTdata <- matrix(nrow=numrow,ncol=numcol)
colnames(FTdata) <- dimnames(FlareTMSD)[[1]]
rownames(FTdata) <- rownames(FTdata, do.NULL= FALSE, prefix= "Sample")


size <- dim(FTdata)[1]

for (i in 1:size){
  total <- 1
  ###Firmicutes
  FTdata[i,1] <- rpearson(n=1, params=FTparFirm)
  total <- total - FTdata[i,1]
  
  ###Actinobacteria
  r <- rpearson(n=1, params=FTparActino)
  FTdata[i,2] <- total*r
  total <- total - FTdata[i,2]
  
  ###Proteobacteria
  r <- rpearson(n=1, params=FTparProteo)
  FTdata[i,3] <- total*r
  total <- total - FTdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=FTparBact)
  FTdata[i,4] <- total*r
  total <- total - FTdata[i,4]
  
  ###Other
  FTdata[i,5] <- total
  
}

apply(FTdata,2,sd)
apply(FlareT,1,sd)

###PostFlare###


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

apply(PFdata,2,sd)
apply(PostFlare,1,sd)

write.table(apply(PFdata,2,sd), "PFdataSD.csv", sep=",")
#write.table(apply(PostFlare,1,sd), "PostFlareSD.csv", sep=",")