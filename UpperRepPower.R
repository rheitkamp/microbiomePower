set.seed(123)
library(HMP)

#Pilot Data
#********************************************
#Loading data 
Fall <- read.csv("pUpperRepFall.csv")
Fallreads <- Fall*11000
Fallreads <- Fallreads[1:12]

colnamesFall <- colnames(Fall)


UpperRepS1 <- read.csv("UpperRepS1.csv")
S1reads <- UpperRepS1*11000
S1reads <- S1reads[,1:12]

colnames(UpperRepS1) <- colnames(Fall)

UpperRepS2 <- read.csv("UpperRepS2.csv")
S2reads <- UpperRepS2*11000
S2reads <- S2reads[,1:12]

colnames(UpperRepS2) <- colnames(Fall)

UpperRepS3 <- read.csv("UpperRepS3.csv")
S3reads <- UpperRepS3*11000
S3reads <- S3reads[,1:12]

colnames(UpperRepS3) <- colnames(Fall)

UpperRepS4 <- read.csv("UpperRepS4.csv")
S4reads <- UpperRepS4*11000
S4reads <- S4reads[,1:12]

colnames(UpperRepS4) <- colnames(Fall)

UpperRepS5 <- read.csv("UpperRepS5.csv")
S5reads <- UpperRepS5*11000
S5reads <- S5reads[,1:12]

colnames(UpperRepS5) <- colnames(Fall)

################################################################################
#############################Siglev=0.01#Pooled#Taxa############################
################################################################################

######
##S1##
######

Gp1=Fallreads
Gp2=S1reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S1')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS1_1P5_175by5s.csv", sep=",")

#######
##S2### <---
#######

Gp1=Fallreads
Gp2=S2reads

#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S2')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS2_1P5_175by5s.csv", sep=",")

######
##S3##
######


Gp1=Fallreads
Gp2=S3reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S3')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS3_1P5_175by5s.csv", sep=",")

######
##S4##
######

Gp1=Fallreads
Gp2=S4reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S4')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS4_1P105_175by5s.csv", sep=",")

######
##S5##
######

Gp1=Fallreads
Gp2=S5reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S5')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS5_1P105_175by5s.csv", sep=",")

################################################################################
##############################Siglev=0.05#Pooled#Taxa###########################
################################################################################

######
##S1##
######

Gp1=Fallreads
Gp2=S1reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S1')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS1_5P5_175by5s.csv", sep=",")

#######
##S2### 
#######

Gp1=Fallreads
Gp2=S2reads

#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S2')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS2_5P5_175by5s.csv", sep=",")

######
##S3##
######


Gp1=Fallreads
Gp2=S3reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S3')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS3_5P5_175by5s.csv", sep=",")

######
##S4##
######

Gp1=Fallreads
Gp2=S4reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S4')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS4_5P5_175by5s.csv", sep=",")

######
##S5##
######

Gp1=Fallreads
Gp2=S5reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S5')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
taxaselection=which(pio>0.01)
Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
nc=dim(Gp1sel)[2]
colnames(Gp1sel)[nc]<-'Pooled taxa'
colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1sel)
C.alpha.multinomial(Gp2sel)


#Fitting DM
fit1=DM.MoM((Gp1sel))
fit2=DM.MoM((Gp2sel))
fit1$reads=colSums(Gp1sel)
fit2$reads=colSums(Gp2sel)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS5_5P5_175by5s.csv", sep=",")


################################################################################
#############################Siglev=0.01#Unpooled#Taxa##########################
################################################################################

######
##S1##
######

Gp1=Fallreads
Gp2=S1reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S1')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS1_1P5_175by5s_NP.csv", sep=",")

#######
##S2###
#######

Gp1=Fallreads
Gp2=S2reads

#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S2')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS2_1P5_175by5s_NP.csv", sep=",")

######
##S3##
######


Gp1=Fallreads
Gp2=S3reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S3')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS3_1P5_175by5s_NP.csv", sep=",")

######
##S4##
######

Gp1=Fallreads
Gp2=S4reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S4')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS4_1P5_175by5s_NP.csv", sep=",")

######
##S5##
######

Gp1=Fallreads
Gp2=S5reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S5')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.01)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS5_1P5_175by5s_NP.csv", sep=",")


################################################################################
#############################Siglev=0.05#Unpooled#Taxa##########################
################################################################################

######
##S1##
######

Gp1=Fallreads
Gp2=S1reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S1')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS1_5P5_175by5s_NP.csv", sep=",")

#######
##S2###
#######

Gp1=Fallreads
Gp2=S2reads

#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S2')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS2_5P5_175by5s_NP.csv", sep=",")

######
##S3##
######


Gp1=Fallreads
Gp2=S3reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S3')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS3_5P5_175by5s_NP.csv", sep=",")

######
##S4##
######

Gp1=Fallreads
Gp2=S4reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S4')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS4_5P5_175by5s_NP.csv", sep=",")

######
##S5##
######

Gp1=Fallreads
Gp2=S5reads


#Preprocessing 1:removing taxa that are zeroes in both groups
#If the data loaded consists of columns given by taxa units and rows subjects (usually called species data),
#then the following preprocessing should be done. Note that you are using the data of the package (ranked data),
#then there is not need to do the following.

comb=rbind(Gp1,Gp2)
index=which(colSums(comb)==0)
if(length(index)>0){
  Gp1=Gp1[,-index]
  Gp2=Gp2[,-index]
}



#Computing and visualizing taxa frequency means for each group
mygroup <- list(Gp1,Gp2)
names(mygroup)<-c('Fall','S5')

Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoM)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)


#Selecting taxa with weigthed average across groups above 1%
#On one part this is to facilitate convergence of the test statistics to the asymptotic
#distribution and being able to use the properties of this distribution to answer the hypothesis

#On a second part, this is also motivated by the uncertainty regarding the abundance of rare taxa, 
#meaning, we do not know if their abundances are correct or not.


fit1=DM.MoM(Gp1)
fit2=DM.MoM(Gp2)
fit1$reads=apply(Gp1,1,sum)
fit2$reads=apply(Gp2,1,sum)

group.data <- list(Gp1, Gp2)
pio <- pioest(group.data)

#pio=pioest(fit1,fit2)
#taxaselection=which(pio>0.01)
#Gp1sel=cbind(Gp1[,taxaselection],rowSums(Gp1[,-taxaselection]))
#Gp2sel=cbind(Gp2[,taxaselection],rowSums(Gp2[,-taxaselection]))
#nc=dim(Gp1sel)[2]
#colnames(Gp1sel)[nc]<-'Pooled taxa'
#colnames(Gp2sel)[nc]<-'Pooled taxa'

#Hypothesis testing;

#Testing DM versus M
C.alpha.multinomial(Gp1)
C.alpha.multinomial(Gp2)


#Fitting DM
fit1=DM.MoM((Gp1))
fit2=DM.MoM((Gp2))
fit1$reads=colSums(Gp1)
fit2$reads=colSums(Gp2)

#Ploting
Sites.pi.MoMsel=rbind(fit1$pi,fit2$pi)
Sites.pi.MoMsel=Sites.pi.MoMsel[,order(abs(Sites.pi.MoMsel[1,]-Sites.pi.MoMsel[2,]),decreasing=T)]
rownames(Sites.pi.MoMsel)<-rownames(Sites.pi.MoM)
matplot(t(Sites.pi.MoMsel),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,0.6))
y=colnames(Sites.pi.MoMsel)
nc=dim((Sites.pi.MoMsel))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,0.6,0.05))
legend("topright",legend=c(rownames(Sites.pi.MoMsel)), cex=0.8, col=c('blue','red'), bty="n",pch=19)
text(-1.5,0.01,srt=90,adj=-1,labels='Taxa Abundance',xpd=T, cex=1)
#text(21,-0.01,srt=0,adj=0,labels='Taxa Abundance',cex=1)

#Comparing the taxa frequency mean across groups
pvalue<-Xmcupo.sevsample(mygroup,dim(Gp2)[2])
pvalue


#Power calculation procedure
pi_2gr=rbind(fit1$pi, fit2$pi) 
group.theta=c(fit1$theta,fit2$theta)

MC=10000 #number of Monte Carlo experiments

ptm <- proc.time()

#Defining range of subjects and sample size 
subjects=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,
           160,165,170,175) #number of subjects per group
nreads=c(11000) #same number of reads per subject (see package to define different number of reads)    
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=list(rep(x[2], x[1]), rep(x[2], x[1])), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table

proc.time() - ptm

write.table(x=power.sample.size.table,file="UpperRepFallvsS5_5P5_175by5s_NP.csv", sep=",")
