##################################################################
#                           RAT DATA                             #
##################################################################
View(apply(LSCCED2, 2, sum))
LSCCED2COrder <- LSCCED2[,order(colSums(LSCCED2), decreasing=TRUE)]
LSCCD4COrder <- LSCCD4[,order(colSums(LSCCD4), decreasing=TRUE)]
LSCED4COrder <- LSCED4[,order(colSums(LSCED4), decreasing=TRUE)]


LSCCED2COrder <- LSCCED2COrder[,-which(colSums(LSCCED2COrder) == 0)]
LSCCD4COrder <- LSCCD4COrder[,-which(colSums(LSCCD4MSD) == 0)]
LSCED4COrder <- LSCED4COrder[,-which(colSums(LSCED4COrder) == 0)]


write.table(LSCCED2COrder, "RatWoundGenusBaselineD2_CbyS.txt", sep = ",")
write.table(LSCCD4COrder, "RatWoundGenusControlsD4_CbyS.txt", sep=",")
write.table(LSCED4COrder, "RatWoundGenusExperimentalD4_CbyS.txt", sep=",")

d2LSCsortedE <- d2LSCsorted[,-which(colSums(d2LSCsorted) == 0)]
d4LSCCsortedE <- d4LSCCsorted[,-which(colSums(d4LSCCsorted) == 0)]
d4LSCEsortedE <- d4LSCEsorted[,-which(colSums(d4LSCEsorted) == 0)]

mean <- apply(d2LSCsortedE,2,mean)
sd <- apply(d2LSCsortedE,2,sd)

LSCCED2MSD <- cbind(mean, sd)

write.table(LSCCED2MSD, "RatWoundGenusBaselineD2MeanSD.txt", sep=",")

mean <- apply(d4LSCCsortedE,2,mean)
sd <- apply(d4LSCCsortedE,2,sd)

LSCCD4MSD <- cbind(mean, sd)

write.table(LSCCD4MSD, "RatWoundGenusControlsD4MeanSD.txt", sep=",")

mean <- apply(d4LSCEsortedE,2,mean)
sd <- apply(d4LSCEsortedE,2,sd)

LSCED4MSD <- cbind(mean, sd)

write.table(LSCED4MSD, "RatWoundGenusExperiementalD4MeanSD.txt", sep=",")

Premainder(LSCCED2MSD)
Premainder(LSCCD4MSD)
Premainder(LSCED4MSD)

write.table(Premainder(LSCCED2MSD), "RatWoundGenusBaselineD2Premainder.txt", sep=",")
write.table(Premainder(LSCCD4MSD), "RatWoundGenusControlsD4Premainder.txt", sep=",")
write.table(Premainder(LSCED4MSD), "RatWoundGenusExperimentalD4Premainder.txt", sep=",")

collectParameters(Premainder(LSCCED2MSD))
collectParameters(Premainder(LSCCD4MSD))
collectParameters(Premainder(LSCED4MSD))

write.table(collectParameters(Premainder(LSCCED2MSD)), "RatWoundGenusBaselineD2Parameters.txt", sep=",")
write.table(collectParameters(Premainder(LSCCD4MSD)), "RatWOundGenusControlsD4Parameters.txt", sep=",")
write.table(collectParameters(Premainder(LSCED4MSD)), "RatWOundGenusExperimentalD4Parameters.txt", sep=",")

########################################################
#                     AD Phylum                        #
########################################################

mean <- apply(pbaseline,2,mean)
sd <- apply(pbaseline,2,sd)

ADPbaselineMSD <- cbind(mean,sd)

write.table(ADPbaselineMSD, "ADPbaselineMSD.txt", sep=",")

mean <- apply(pcontrol,2,mean)
sd <- apply(pcontrol,2,sd)

ADPcontrolMSD <- cbind(mean, sd)

write.table(ADPcontrolMSD, "ADPcontrolMSD.txt", sep=",")

mean <- apply(pflarent,2,mean)
sd <- apply(pflarent,2,sd)

ADPflarentMSD <- cbind(mean,sd)

write.table(ADPflarentMSD, "ADPflarentMSD.txt", sep=",")

mean <- apply(pflaret, 2, mean)
sd <- apply(pflaret, 2, sd)

ADPflaretMSD <- cbind(mean,sd)

write.table(ADPflaretMSD, "ADPflaretMSD.txt", sep=",")

mean <- apply(ppostflare,2,mean)
sd <- apply(ppostflare,2,sd)

ADPpostflareMSD <- cbind(mean,sd)

write.table(ADPpostflareMSD, "ADPpostflare.txt", sep=",")

Premainder(ADPbaselineMSD)
Premainder(ADPcontrolMSD)
Premainder(ADPflarentMSD)
Premainder(ADPflaretMSD)
Premainder(ADPpostflareMSD)

write.table(Premainder(ADPbaselineMSD), "ADPBaselinePremainder.txt", sep=",")
write.table(Premainder(ADPcontrolMSD), "ADPControlPremainder.txt", sep=",")
write.table(Premainder(ADPflarentMSD), "ADPFlarentPremainder.txt", sep=",")
write.table(Premainder(ADPflaretMSD), "ADPFlaretPremainder.txt", sep=",")
write.table(Premainder(ADPpostflareMSD), "ADPPostflarePremainder.txt", sep=",")

collectParameters(Premainder(ADPbaselineMSD))
collectParameters(Premainder(ADPcontrolMSD))
collectParameters(Premainder(ADPflarentMSD))
collectParameters(Premainder(ADPflaretMSD))
collectParameters(Premainder(ADPpostflareMSD))

write.table(collectParameters(Premainder(ADPbaselineMSD)), "ADPBaselineParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADPcontrolMSD)), "ADPControlParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADPflarentMSD)), "ADPFlarentParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADPflaretMSD)), "ADPFlaretParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADPpostflareMSD)), "ADPPostflaretParameters.txt", sep=",")

############################################################
#                       AD Species                         #
############################################################

mean <- apply(sbaseline,2,mean)
sd <- apply(sbaseline,2,sd)

ADSbaselineMSD <- cbind(mean,sd)

write.table(ADSbaselineMSD, "ADSBaselineMSD.txt", sep=",")

mean <- apply(scontrol,2,mean)
sd <- apply(scontrol,2,sd)

ADScontrolMSD <- cbind(mean,sd)

write.table(ADScontrolMSD, "ADSControlsMSD.txt", sep=",")

mean <- apply(sflarent,2,mean)
sd <- apply(sflarent,2,sd)

ADSflarentMSD <- cbind(mean,sd)

write.table(ADSflarentMSD, "ADSFlarentMSD.txt", sep=",")

mean <- apply(sflaret,2,mean)
sd <- apply(sflaret,2,sd)

ADSflaretMSD <- cbind(mean, sd)

write.table(ADSflaretMSD, "ADSFlaretMSD.txt", sep=",")

mean <- apply(spostflare,2,mean)
sd <- apply(spostflare,2,sd)

ADSpostflareMSD <- cbind(mean,sd)

write.table(ADSpostflareMSD, "ADSPostflareMSD.txt", sep=",")

Premainder(ADSbaselineMSD)
Premainder(ADScontrolMSD)
Premainder(ADSflarentMSD)
Premainder(ADSflaretMSD)
Premainder(ADSpostflareMSD)

write.table(Premainder(ADSbaselineMSD), "ADSBaselinePremainser.txt", sep=",")
write.table(Premainder(ADScontrolMSD), "ADSControlPremainder.txt", sep=",")
write.table(Premainder(ADSflarentMSD), "ADSFlarentPremainder.txt", sep=",")
write.table(Premainder(ADSflaretMSD), "ADFlaretPremainder.txt", sep=",")
write.table(Premainder(ADSpostflareMSD), "ADSPostflarePremainder.txt", sep=",")

collectParameters(Premainder(ADSbaselineMSD))
collectParameters(Premainder(ADScontrolMSD))
collectParameters(Premainder(ADSflarentMSD))
collectParameters(Premainder(ADSflaretMSD))
collectParameters(Premainder(ADSpostflareMSD))

write.table(collectParameters(Premainder(ADSbaselineMSD)), "ADSBaselineParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADScontrolMSD)), "ADSControlParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADSflarentMSD)), "ADSFlarentParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADSflaretMSD)), "ADSFlaretParameters.txt", sep=",")
write.table(collectParameters(Premainder(ADSpostflareMSD)), "ADSPostflareParameters.txt", sep=",")
