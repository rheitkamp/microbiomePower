#Loading Library
#**********************************************
library(HMP)

#Read CSV files
#**********************************************
data=read.csv(file="miceRDP_phylum.csv", header = TRUE)
metadata=read.csv(file="metadata.csv", header = TRUE) 


#Obtaining the groups of interest
#**********************************************
ObObpopulation=metadata[metadata$Genotype=="ob/ob",2]
Leanpopulation=metadata[metadata$Genotype=="+/+",2]

dataObOb=data[,as.vector(ObObpopulation)]
rownames(dataObOb)<-as.vector(data[,1])
dataLean=data[,as.vector(Leanpopulation)]
rownames(dataLean)<-as.vector(data[,1])

#Multinomial versus DM
#**********************************************
#The two most prominent
pvalue1=C.alpha.multinomial(t(dataObOb)[,c(1:2)])
pvalue2=C.alpha.multinomial(t(dataLean)[,c(1:2)])

weirMoM(t(dataObOb), se=TRUE) #overdispersion with standard deviation
weirMoM(t(dataObOb), se=TRUE) #overdispersion with standard deviation

#Two sample RAD probability-mean comparison 
#************************************
#The two most abundant taxa
mygroup <- list(t(dataObOb)[,c(1:2)],t(dataLean)[,c(1:2)])
pvalue<-Xmcupo.sevsample(mygroup, 2)

#The least most abundant taxa
mygroup <- list(t(dataObOb)[,c(3:7)],t(dataLean)[,c(3:7)])
pvalue<-Xmcupo.sevsample(mygroup, 5)
DM.MoM(t(dataObOb)[,c(3:7)])
DM.MoM(t(dataLean)[,c(3:7)])


#All taxa 
pvalue<-Xmcupo.sevsample(mygroup,7)
DM.MoM(t(dataObOb))
DM.MoM(t(dataLean))

 
#Sample size calculation using MC procedures
#********************************************
fit1=DM.MoM(t(dataObOb))
fit2=DM.MoM(t(dataLean))

#Choosing subtaxa to perform comparison guidelines
#To compare common taxa, namely, pi>0
#To compare taxa that combines have taxa probability > 0.01


#Selecting common taxa with total counts >5 and colapsing the rest into one taxa
dataObObsel=rbind(dataObOb[which(rowSums(dataLean)>5),],colSums(dataObOb[-which(rowSums(dataLean)>5),]))
dataLeansel=rbind(dataLean[which(rowSums(dataLean)>5),],colSums(dataLean[-which(rowSums(dataLean)>5),]))


#Subselection
dataObObsel=rbind(dataObObsel[1:4,],colSums(dataObObsel[-c(1:4),]))
dataLeansel=rbind(dataLeansel[1:4,],colSums(dataLeansel[-c(1:4),]))

fit1=DM.MoM(t(dataObObsel))
fit2=DM.MoM(t(dataLeansel))

reads=200;subjects=20
fit1$reads=rep(reads,subjects)
fit2$reads=rep(reads,subjects)
subjects*reads*pioest(fit1,fit2)


#Bar plots
taxa.prob=cbind(fit1$pi,fit2$pi); colnames(taxa.prob)=c("Obese","Lean")
barplot(taxa.prob,beside=TRUE, col=rainbow(length(fit1$pi)))
legend("topright",rownames(taxa.prob), cex=0.6,  bty="n", fill=rainbow(length(fit1$pi)));
pdiff=(fit2$pi-fit1$pi)/fit2$pi

pi_2gr=rbind(fit1$pi, fit2$pi) #RAD probability-mean
group.theta=c(fit1$theta,fit2$theta) #overdispersion

MC=1000 #number of Monte Carlo experiments

#Defining range of subjects and sample size 
subjects=c(4,6,8,9,10,11,12) #number of subjects per group
nreads=c(1000,2000,4000,5000) #same number of reads per subject (see package to define different number of reads)		
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

#MC procedure to compute power for range of subjects and reads per sample
power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=rep(x[2], x[1]), MC, pi0=fit1$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)

#MC procedures for computing the probability that all taxa are greater than 5
MC=1000

realization=apply(MCr,1,function(x,reads,P,gammas){
dirmult_data <- Dirichlet.multinomial(Nrs=rep(reads,P),gammas)
Kest=sum(colSums(dirmult_data)> 5)== sum(gammas>0)},reads=reads,P=P,gammas=fit2$gamma)
prob=sum(realization)/MC
prob


Kest=prod()

dirmult_data <- Dirichlet.multinomial(Nrs=rep(reads*P,100), fitfit1$gamma)
question = dirmult_data>5

mult_data <- Multinomial(Nrs=(1000*10), fitfit1$pi)
mult_data



#Read CSV files
#**********************************************
data=read.csv(file="miceRDP_genus.csv", header = TRUE)
data=cbind(data[,1],data[rowSums(data[,-1])!=0,2:dim(data)[2]])

metadata=read.csv(file="metadata.csv", header = TRUE) 




#Obtaining the groups of interest
#**********************************************
ObObpopulation=metadata[metadata$Genotype=="ob/ob",2]
Leanpopulation=metadata[metadata$Genotype=="+/+",2]

dataObOb=data[,as.vector(ObObpopulation)]
rownames(dataObOb)<-as.vector(data[,1])
dataLean=data[,as.vector(Leanpopulation)]
rownames(dataLean)<-as.vector(data[,1])


#Multinomial versus DM
#**********************************************

fit1=DM.MoM(t(dataObOb))
fit2=DM.MoM(t(dataLean))
fit1$reads=colSums(dataObOb)
fit2$reads=colSums(dataLean)
effectsize(fit1,fit2)


#Selecting taxa > 0.01

reads=1000;subjects=20
fit1$reads=rep(reads,subjects)
fit2$reads=rep(reads,subjects)
pio=pioest(fit1,fit2)


dataObObsel=rbind(dataObOb[which(pio>0.01),],colSums(dataObOb[-which(pio>0.01),]))
dataLeansel=rbind(dataLean[which(pio>0.01),],colSums(dataLean[-which(pio>0.01),]))


fit1=DM.MoM(t(dataObObsel))
fit2=DM.MoM(t(dataLeansel))


fit1$reads=colSums(dataObObsel)
fit2$reads=colSums(dataLeansel)
effectsize(fit1,fit2)

#Bar plots
taxa.prob=cbind(fit1$pi,fit2$pi); colnames(taxa.prob)=c("Obese","Lean")
barplot(taxa.prob,beside=TRUE, col=rainbow(length(fit1$pi)))

legend("topright",rownames(taxa.prob), cex=0.6,  bty="n", fill=rainbow(length(fit1$pi)));

x11()
prop=fit1$pi

plot(prop,axes="FALSE", ann="FALSE",ylim=c(0,0.7),pch=19,col="red",type="o")
title("Taxa Abundances at Phylum level for Vervet on TAD and Chow")
y=strsplit(as.character(names(fit1$pi)), " . ")
y=unlist(lapply(y, function(x) x[1]))
#y[c(16,17)]=c("Unclassified","Unclassified")
axis(1,at=1:length(fit1$pi), lab=F, hadj=1)
text(1:length(fit1$pi),-0.05, srt=90, adj=1,labels=y,xpd=T, cex=0.5)
axis(2, las=1, at=seq(0,0.7,0.05))
prop=fit2$pi
points(prop,pch=19,col="blue",type="o")
points(pio,pch=19,col="green",type="o")

legend("topright",c("Obese","Lean","weighted average"),pch=19,col= c("red","blue","green"))


pdiff=(fit2$pi-fit1$pi)/fit2$pi




pi_2gr=rbind(fit1$pi, fit2$pi) #RAD probability-mean
group.theta=c(fit1$theta,fit2$theta) #overdispersion



MC=1000 #number of Monte Carlo experiments

#Defining range of subjects and sample size 
subjects=c(4,6,8,9,10,11,12) #number of subjects per group
nreads=c(1000,2000,4000,5000) #same number of reads per subject (see package to define different number of reads)		

subjects=c(10) #number of subjects per group
nreads=c(1000) #same number of reads per subject (see package to define different number of reads)		

subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

#MC procedure to compute power for range of subjects and reads per sample
power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=rep(x[2], x[1]), MC, pi0=fit1$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)


prob1=apply(range,1,function(y,param){
realization=apply(as.matrix(rep(1,MC)),1,function(x,reads,P,gammas){
dirmult_data <- Dirichlet.multinomial(Nrs=rep(reads,P),gammas)
Kest=sum(colSums(dirmult_data)> 0)== sum(gammas>0)},reads=y[2],P=y[1],gammas=param)
prob=sum(realization)/MC
prob
},param=fit1$gamma)
prob.range=cbind(range,prob1)

prob2=apply(range,1,function(y,param){
realization=apply(as.matrix(rep(1,MC)),1,function(x,reads,P,gammas){
dirmult_data <- Dirichlet.multinomial(Nrs=rep(reads,P),gammas)
Kest=sum(colSums(dirmult_data)> 0)== sum(gammas>0)},reads=y[2],P=y[1],gammas=param)
prob=sum(realization)/MC
prob
},param=fit2$gamma)
prob2.range=cbind(range,prob2)

prob.k=prob1*prob2+prob1*(1-prob2)+prob2*(1-prob1)
prob.k.range=cbind(range,prob.k)
power.sample.size2=cbind(range,power*prob.k)
power.sample.size2





