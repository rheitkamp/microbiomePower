set.seed(1234) #always use the same seed during testing
library(PearsonDS) #use pearsons distribution library (beta distribution)
library(HMP) #use HMP package

#####################
#READ DATA FROM FILES
#####################

(Control <- read.table("ADControls.txt"))
(Baseline <- read.table("ADBaseline.txt"))
(FlareNT <- read.table("ADFlareNT.txt"))
(FlareT <- read.table("ADFlareT.txt"))
(PostFlare <- read.table("ADPostFlare.txt"))

colSums(Control) #sum of all columns, should equal 1
colSums(Baseline)
colSums(FlareNT)
colSums(FlareT)
colSums(PostFlare)

#######apply built-in sd function to rows of table (1) and cast as matrix with rows as 
#######table rows and column as "SD"
(ControlSD <- matrix(apply(Control,1,sd),dimnames=list(rownames(Control), "SD")))
(BaselineSD <- matrix(apply(Baseline,1,sd), dimnames=list(rownames(Baseline), "SD")))
(FlareNTSD <- matrix(apply(FlareNT,1,sd), dimnames=list(rownames(FlareNT), "SD")))
(FlareTSD <- matrix(apply(FlareT,1,sd), dimnames=list(rownames(FlareT), "SD")))
(PostFlareSD <- matrix(apply(PostFlare,1,sd), dimnames=list(rownames(PostFlare), "SD")))

#######create matrix of with built-in function rowMeans on table with 1 row, 5 columns, 
#######rows as "Mean" and columns as table rows
(ControlMean <- matrix(rowMeans(Control), nrow=1, ncol=5, dimnames=list("Mean", rownames(Control))))
(BaselineMean <- matrix(rowMeans(Baseline), nrow=1, ncol=5, dimnames=list("Mean", rownames(Baseline))))
(FlareNTMean <- matrix(rowMeans(FlareNT), nrow=1, ncol=5, dimnames=list("Mean", rownames(FlareNT))))
(FlareTMean <- matrix(rowMeans(FlareT), nrow=1, ncol=5, dimname=list("Mean", rownames(FlareT))))
(PostFlareMean <- matrix(rowMeans(PostFlare), nrow=1, ncol=5, dimnames=list("Mean", rownames(PostFlare))))

#######declare a function that creates a matrix of space-filled means (ie percent remainder)
######
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

########
#run Premainter on percentages in proportions-of-total matrices
########
(ControlMP <- Premainder(ControlMean*100))
(BaselineMP <- Premainder(BaselineMean*100))
(FlareNTMP <- Premainder(FlareNTMean*100))
(FlareTMP <- Premainder(FlareTMean*100))
(PostFlareMP <- Premainder(PostFlareMean*100))

################
#make table of percent-remainder means with proportions-of-total standard deviations
###############
(ControlMSD <- cbind(ControlMP, ControlSD))
(BaselineMSD <- cbind(BaselineMP, BaselineSD))
(FlareNTMSD <- cbind(FlareNTMP, FlareNTSD))
(FlareTMSD <- cbind(FlareTMP, FlareTSD))
(PostFlareMSD <- cbind(PostFlareMP, PostFlareSD))

#######
#declare a function to calculate beta parameters from a mean and standard deviation
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
#calculate the beta parameters for each dataset for each taxon
#######
collectParameters <- function(datamatrix) {
  taxa <- rownames(datamatrix)
  parameterNames <- list("type","a","b","location","scale")
  BetaParameters <- matrix(data = NA, nrow = length(taxa), ncol = length(parameterNames), dimnames=list(taxa, parameterNames))
  
  for (i in 1:nrow(BetaParameters)) {
    parameters <- getBetaParams(datamatrix[i,1],datamatrix[i,2]);
    BetaParameters[i,] <- c(parameters$type, parameters$a, parameters$b, parameters$location, parameters$scale)
  }
  
  return(BetaParameters)
}

controlPar <- collectParameters(ControlMSD)
baselinePar <- collectParameters(BaselineMSD)
flareNTPar <- collectParameters(FlareNTMSD)
flareTPar <- collectParameters(FlareTMSD)
postFlarePar <- collectParameters(PostFlareMSD)

###############
#calculate percent-remainder data for each sample
###############
spaceFill <- function (dataMatrix, distParameters, subjects){
  subjects <- subjects #number of subjects for rows
  taxa <- dim(dataMatrix)[1] #number of taxa for columns 
  
  Cdata <- matrix(nrow=subjects,ncol=taxa)
  taxaNames <- dimnames(dataMatrix)[[1]]
  colnames(Cdata) <- taxaNames #list taxa names in columns
  rownames(Cdata) <- rownames(Cdata, do.NULL= FALSE, prefix= "Sample") #call rows samples
  
  for (i in 1:subjects){
    total <- 1
    for (j in 1:(taxa-1)){
      r <- rpearson(n=1, params=distParameters[j,]);
      Cdata[i,j] <- total * r
      total <- total - Cdata[i,j];
    }
    Cdata[i,taxa] <- total;
  }
  return(Cdata)
}

Cdata <- spaceFill(ControlMSD,controlPar,4)
Bdata <- spaceFill(BaselineMSD, baselinePar,6)
Fdata <- spaceFill(FlareNTMSD,flareNTPar,10)
Tdata <- spaceFill(FlareTMSD,flareTPar,2)
Pdata <- spaceFill(PostFlareMSD,postFlarePar,15)

######
#display data as a barchart
######

Barchart.data(Cdata, title="ADControl")
Barchart.data(Bdata, title="ADBaseline")
Barchart.data(Fdata, title="ADFlare")
Barchart.data(Tdata, title="ADTreatment")
Barchart.data(Pdata, title="ADPost")
#barchart(x=Cdata,horizontal=FALSE, col=rainbow(5),xlab=NULL)

####
#get standard deviation of simulated samples
####
apply(Cdata,2,sd)
apply(Control,1,sd)

colMeans(Cdata)
rowMeans(Control)

###Baseline### (repeat)

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

###FlareNT###  (repeat)

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
Barchart.data(FNTdata, title="ADFlareNT")
apply(FNTdata,2,sd)
apply(FlareNT,1,sd)

###FlareT###  (repeat)

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

###PostFlare### (repeat)


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