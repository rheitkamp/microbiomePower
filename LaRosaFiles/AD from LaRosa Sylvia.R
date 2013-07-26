library(HMP)


#Pilot Data
#********************************************
#Loading data 
controls = read.table("ADcontrolReads.txt")
flares = read.table("FlareReads.txt")

Gp1=controls
Gp2=flares


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
names(mygroup)<-c('control','flares')

#
#
#
#start here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#
#


Sites.pi.MoM=do.call(rbind,lapply(mygroup,function(x) DM.MoM(x)$pi))
Sites.pi.MoM=Sites.pi.MoM[,order(abs(Sites.pi.MoM[1,]-Sites.pi.MoM[2,]),decreasing=T)]

#Ploting taxa frequency means
matplot(t(Sites.pi.MoM),type='o',pch=19,axes="FALSE", ann="FALSE",col=c('blue','red'),ylim=c(-0.25,1.0))
y=colnames(Sites.pi.MoM)
nc=dim((Sites.pi.MoM))[2]
axis(1,at=1:nc, lab=F, hadj=3,pos=-0.009)
text(1:nc,-0.02, srt=90,adj=1,labels=y,xpd=T, cex=1)
axis(2, las=1, at=seq(0,1,0.1))
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

#Defining range of subjects and sample size 
subjects=c(5,10,20,25) #number of subjects per group
nreads=c(50000) #same number of reads per subject (see package to define different number of reads)		
subject.reads = expand.grid(subj=subjects, reads=nreads) 
range=subject.reads

power=apply(range,1,function(x){MC.Xmcupo.statistics(Nrs=rep(x[2], x[1]), MC, pi0=fit2$pi, group.pi=pi_2gr, group.theta=group.theta, type = "ha", siglev = 0.05)})
power.sample.size=cbind(range,power)
power.sample.size.table=xtabs(power~subj+reads,data=power.sample.size)
power.sample.size.table


#Functions that need to be overwritten


#Effect Size function
#********************************************
effectsize=function(par1,par2){
	Kc=length(par1$pi)
	group.parameter.estimated <- list()
	group.parameter.estimated[[1]] = c(par1$pi, par1$theta, par1$reads)
	group.parameter.estimated[[2]] = c(par2$pi, par2$theta, par2$reads)
	Xmcupo = Xmcupo.statistics(group.parameter.estimated, K=Kc)


	#Maximum distance between taxa proportions
	pi1=c(1,0);pi2=c(0,1)
	gp.estimated <- list()
	gp.estimated[[1]] = c(pi1, par1$theta, par1$reads)
	gp.estimated[[2]] = c(pi2, par2$theta, par2$reads)
	Xmcupo.gp = Xmcupo.statistics(gp.estimated, K=2)

	#Phi=sqrt(Xmcupo/(sum(c(par1$reads,par2$reads))))
	#Phi=sqrt(Xmcupo/(length(par1$reads)+length(par2$reads)))
	return(c(Xmcupo[1],Xmcupo.gp[1],sqrt(Xmcupo[1]/Xmcupo.gp[1])))
}

#Finding the common pi0 for two sample comparison
#********************************************
pioest= function (par1,par2){
    
    K=length(par1$pi)
    group.parameter <- list()
    group.parameter[[1]] = c(par1$pi, par1$theta, par1$reads)
    group.parameter[[2]] = c(par2$pi, par2$theta, par2$reads)

    n.groups = length(group.parameter)
    index = as.matrix(seq(1:n.groups))
    Xscg = apply(index, 1, function(x) {
        pi = group.parameter[[x]][1:K]
        theta = group.parameter[[x]][K + 1]
        P = length(group.parameter[[x]])
        nreads.data = group.parameter[[x]][(K + 2):P]
        N_1Cj = ((theta * (sum(nreads.data^2) - sum(nreads.data)) + 
            sum(nreads.data))/(sum(nreads.data))^2)
        Out1 = c(pi, N_1Cj)
        Out1
    })
    Xscg=Xscg[rowSums(Xscg)!=0,];
    K=dim(Xscg)[1]-1
    pi0 = (Xscg[1:K, ] %*% as.matrix(1/Xscg[K + 1, ]))/sum(1/Xscg[K +1, ])
    return(pi0)
}



 weirMoM= function (data, se = FALSE) 
{
    K <- ncol(data)
    J <- nrow(data)
	
	colSumsdata=apply(data,2,function(z){sum(as.numeric(z))})
	rowSumsdata=apply(data,1,function(z){sum(as.numeric(z))})
	totalsample=sum(apply(data,2,function(z){sum(as.numeric(z))}))
    
	
    MoM <-colSumsdata/totalsample
    Sn <- rowSumsdata
	
    MSP <- (J - 1)^(-1) * sum(rowSums((data/rowSumsdata - matrix(rep(MoM, 
        J), J, K, byrow = T))^2) * Sn)
    MSG <- (totalsample - J)^(-1) * sum(rowSums(data/rowSumsdata * 
        (1 - data/rowSumsdata)) * Sn)
    nc <- 1/(J - 1) * (sum(Sn) - sum(Sn^2)/sum(Sn))
    MoM.wh <- (MSP - MSG)/(MSP + (nc - 1) * MSG)
    if (se) {
        std.er <- sqrt(2 * (1 - MoM.wh)^2/(J - 1) * ((1 + (nc - 
            1) * MoM.wh)/nc)^2)
        list(theta = MoM.wh, se = std.er)
    }
    else MoM.wh
}



Xmcupo.statistics.H = function(K,group.parameter){                       
                n.groups=length(group.parameter)
                index=as.matrix(seq(1:n.groups))

                group.parameter.estimated<-list()
                for(x in index){
                                # Dirichlet-Multinomial
                                pi=group.parameter[[x]][1:K];
                                theta=group.parameter[[x]][K+1];
                                P=length(group.parameter[[x]])
                                nreads.data=as.matrix(group.parameter[[x]][(K+2):P]);
                                data=Dirichlet.multinomial(nreads.data,shape=pi*(1-theta)/theta)

                                #MoM of theta and Unbiased estimator of pi
								totalsample=sum(apply(data,2,function(z){sum(as.numeric(z))}))
                                pi.MoM=apply(data,2,function(z){sum(as.numeric(z))})/totalsample
                                ## replace zeros 
                                pi.MoMb=pi.MoM
                                
                                q=pi.MoM
                                r=length(as.vector(which(q==0)))
                                rr=length(q)-r
                                q[which(q!=0)]=q[which(q!=0)]-r/(rr*2*(totalsample+1))
                                q[which(q==0)]=1/(2*(totalsample+1))
                                pi.MoMb=q
                                                
                                

                                theta.MoM=weirMoM(data)
                                                                                                
                                group.parameter.estimated[[x]]=c(pi.MoMb,theta.MoM,t(nreads.data))                                           
                }
                                                                                
                #Xmcupo.statistics
                Xmcupo=Xmcupo.statistics(group.parameter.estimated,K)                                                                                                                         
}


MC.Xmcupo.statistics=function (Nrs, MC, pi0, group.pi, group.theta, type = "hnull", 
    siglev = 0.05) 
{
    if (missing(Nrs)) {
        stop("Nrs missing")
    }
    if (missing(MC)) {
        stop("MC missing")
    }
    if (missing(group.theta)) {
        stop("group.theta missing")
    }
    if (missing(pi0)) {
        stop("pi0 missing")
    }
    if (missing(group.pi) && type == "ha") {
        stop("group.pi missing")
    }
    #Nrs = t(t(Nrs))
    MCC = as.matrix(seq(1, 1, length = MC))
    n.groups = length(group.theta)
    index = as.matrix(seq(1:n.groups))
    group.parameter = list()
    K = length(pi0)
    if (type == "hnull") {
        for (i in index) {
            group.parameter[[i]] = c(pi0, group.theta[i], Nrs[[i]])
        }
    }
    else if (type == "ha") {
        for (i in index) {
            group.parameter[[i]] = c(group.pi[i, ], group.theta[i], 
                Nrs[[i]])
        }
    }
    else {
        stop(paste("Can't find type ", type, sep = ""))
    }
    Xmcupot = t(apply(MCC, 1, function(x) {
        Xmcupo.statistics.H(K, group.parameter)
    }))
    Xmcupo = t(Xmcupot)
    dgf = (length(group.theta) - 1) * (K - 1)
    Xmcupo_pval = apply(Xmcupo, 2, function(t) {
        q.alpha = qchisq(p = (1 - siglev), df = dgf, ncp = 0, 
            lower.tail = TRUE)
        sum((t[t != "NaN"] > q.alpha)/(sum(t != "NaN")))
    })
}


