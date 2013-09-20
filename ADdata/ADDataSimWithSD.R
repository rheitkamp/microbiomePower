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

######
#(controlmean <- matrix(colMeans(control), dimnames=list(colnames(control), "Mean")))
#(baselinemean <- matrix(colMeans(baseline), dimnames=list(colnames(baseline), "Mean")))
#(flarentmean <- matrix(colMeans(flarent), dimnames=list(colnames(flarent), "Mean")))
#(flaretmean <- matrix(colMeans(flaret), dimnames=list(colnames(flaret), "Mean")))
#(postflaremean <- matrix(colMeans(postflare), dimnames=list(colnames(postflare), "Mean")))

(controlSD <- matrix(apply(control,2,sd), dimnames=list(colnames(control), "SD")))
(baselineSD <- matrix(apply(baseline,2,sd), dimnames=list(colnames(baseline), "SD")))
(flarentSD <- matrix(apply(flarent,2,sd), dimnames=list(colnames(flarent), "SD")))
(flaretSD <- matrix(apply(flarent,2,sd), dimnames=list(colnames(flaret), "SD")))
(postflareSD <- matrix(apply(postflare,2,sd), dimnames=list(colnames(postflare), "SD")))

######
Premainder <- function(x) { ###data is formated in percentages 
  rowsize <- dim(x)[1] ###gets the dimensions for row
  colsize <- dim(x)[2] ###gets the dimensions for col
  y <- matrix(nrow=rowsize, ncol=colsize, dimnames=list(rownames(x),colnames(x))) ###sets up empty matrix
  for(i in 1:rowsize){ ###for each row
    total <- apply(X=x,MARGIN=1,sum)[i] ###sums the total of the row
    y[i,1] <- x[i,1]/100 ### fills in the first column with the portions
    total <- total - x[i,1] ###number from first column is subtracted from the total 
    for(k in 2:colsize){ ###for each column
      y[i,k] <- x[i,k]/ total ###number from 2nd column is divided by the total and stored into new matrix
      total <- total - x[i,k] ###subtract number above from total...this continues 
    }
  } 
  return(y) ###returns new matrix
}

########
(pcontrol <- Premainder(control*100))
(pbaseline <- Premainder(baseline*100))
(pflarent <- Premainder(flarent*100))
(pflaret <- Premainder(flaret*100))
(ppostflare <- Premainder(postflare*100))

#####
pcontrolmean <- colMeans(pcontrol)
pbaselinemean <- colMeans(pbaseline)
pflarentmean <- colMeans(pflarent)
pflaretmean <- colMeans(pflaret)
ppostflare <- colMeans(ppostflare)

pcontrolSD <- apply(pcontrol,2,sd)
pbaselineSD <- apply(pbaseline,2,sd)
pflarentSD <- apply(pflarent,2,sd)
pflaretSD <- apply(pflaret,2,sd)
ppostflare <- apply(ppostflare,2,sd)

#####
(pcontrolmean <- matrix(colMeans(pcontrol), dimnames=list(colnames(pcontrol), "Mean")))
(pbaselinemean <- matrix(colMeans(pbaseline), dimnames=list(colnames(pbaseline), "Mean")))
(pflarentmean <- matrix(colMeans(pflarent), dimnames=list(colnames(pflarent), "Mean")))
(pflaretmean <- matrix(colMeans(pflaret), dimnames=list(colnames(pflaret), "Mean")))
(ppostflaremean <- matrix(colMeans(ppostflare), dimnames=list(colnames(ppostflare), "Mean")))

(pcontrolSD <- matrix(apply(pcontrol,2,sd), dimnames=list(colnames(pcontrol), "SD")))
(pbaselineSD <- matrix(apply(pbaseline,2,sd), dimnames=list(colnames(pbaseline), "SD")))
(pflarentSD <- matrix(apply(flarent,2,sd), dimnames=list(colnames(pflarent), "SD")))
(pflaretSD <- matrix(apply(flaret,2,sd), dimnames=list(colnames(pflaret), "SD")))
(ppostflareSD <- matrix(apply(ppostflare,2,sd), dimnames=list(colnames(ppostflare), "SD")))

######Before changing the matrix to a 5x1
(controlmsd <- rbind(controlmean, pcontrolmean, controlSD, pcontrolSD))
(baselinemsd <- rbind(baselinemean, pbaselinemean, baselineSD, pbaselineSD))
(flarentmsd <- rbind(flarentmean, pflarentmean, flarentSD, pflarentSD))
(flaretmsd <- rbind(flaretmean, pflaretmean, flaretSD, pflaretSD))
(postflaremsd <- rbind(postflaremean, ppostflaremean, postflareSD, ppostflareSD))

######After changing to a 5x1 and binding Mean and SD together in to one matrix
(pcontrolmsd <- cbind(pcontrolmean, pcontrolSD))
(pbaselinemsd <- cbind(pbaselinemean, pbaselineSD))
(pflarentmsd <- cbind(pflarentmean, pflarentSD))
(pflaretmsd <- cbind(pflaretmean, pflaretSD))
(ppostflaremsd <- cbind(ppostflaremean, ppostflareSD))

#######
getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)*(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}

#######Parameter for pcontrols
(pcparfirm <- getBetaParams(pcontrolmsd[1,1], pcontrolmsd[1,2]))
(pcparactino <- getBetaParams(pcontrolmsd[2,1], pcontrolmsd[2,2]))
(pcparproteo <- getBetaParams(pcontrolmsd[3,1], pcontrolmsd[3,2]))
(pcparbact <- getBetaParams(pcontrolmsd[4,1], pcontrolmsd[4,2]))

######Parameters for pbaseline
(pbparfirm <- getBetaParams(pbaselinemsd[1,1], pbaselinemsd[1,2]))
(pbparactino <- getBetaParams(pbaselinemsd[2,1], pbaselinemsd[2,2]))
(pbparproteo <- getBetaParams(pbaselinemsd[3,1], pbaselinemsd[3,2]))
(pbparbact <- getBetaParams(pbaselinemsd[4,1], pbaselinemsd[4,2]))

######Parameters for pflarent
(pfntparfirm <- getBetaParams(pflarentmsd[1,1], pflarentmsd[1,2]))
(pfntparactino <- getBetaParams(pflarentmsd[2,1], pflarentmsd[2,2]))
(pfntparproteo <- getBetaParams(pflarentmsd[3,1], pflarentmsd[3,2]))
(pfntparbact <- getBetaParams(pflarentmsd[4,1], pflarentmsd[4,2]))

######Parameter for pflaret
(pftparfirm <- getBetaParams(pflaretmsd[1,1], pflaretmsd[1,2]))
(pftparactino <- getBetaParams(pflaretmsd[2,1], pflaretmsd[2,2]))
(pftparproteo <- getBetaParams(pflaretmsd[3,1], pflaretmsd[3,2]))
(pftparbact <- getBetaParams(pflaretmsd[4,1], pflaretmsd[4,2]))

######Parameter for pbaseline
(ppfparfirm <- getBetaParams(ppostflaremsd[1,1], ppostflaremsd[1,2]))
(ppfparactino <- getBetaParams(ppostflaremsd[2,1], ppostflaremsd[2,2]))
(ppfparproteo <- getBetaParams(ppostflaremsd[3,1], ppostflaremsd[3,2]))
(ppfparbact <- getBetaParams(ppostflaremsd[4,1], ppostflaremsd[4,2]))

######pcontrol
numrow <- 25 ###number of participants/subjects
numcol <- dim(pcontrolmsd)[1] ###determines the number of columns
pcdata <- matrix(nrow=numrow,ncol=numcol) ###create new empty matrix
colnames(pcdata) <- dimnames(pcontrolmsd)[[1]] ###assigning colnames on empty matrix
rownames(pcdata) <- rownames(pcdata, do.NULL= FALSE, prefix= "Sample") ###assigning rownames for empty matrix


size <- dim(pcdata)[1] ###dimension of new matrix, number of rows

for (i in 1:size){
  total <- 1 ###starting total (can change this to reflect total in starting dataset)
  ###Firmicutes
  pcdata[i,1] <- rpearson(n=1, params=pcparfirm) ###simulate datapoint for Firm and store in empty matrix
  total <- total - pcdata[i,1] ###subtract number simulated from total
  
  ###Actinobacteria
  r <- rpearson(n=1, params=pcparactino) ###number is simulated and stored as r
  pcdata[i,2] <- total*r ###number simulated (r) is multiplied by the total to get actual proportion and store in new matrix
  total <- total - pcdata[i,2] ###number stored in new matrix is subtracted from total
  
  ###Proteobacteria
  r <- rpearson(n=1, params=pcparproteo)
  pcdata[i,3] <- total*r
  total <- total - pcdata[i,3]
  
  ###Bacteriodetes
  r <- rpearson(n=1, params=pcparbact)
  pcdata[i,4] <- total*r
  total <- total - pcdata[i,4]
  
  ###Other
  pcdata[i,5] <- total ###whatever is leftover is stored as others
  
}

pcdata
colMeans(pcdata)
apply(pcdata,2,sd)
controlmean
controlSD