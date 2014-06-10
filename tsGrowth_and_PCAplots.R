#This script looks at the differences in tree swallow chick fatty acids across sites and makes comparisons to bolus fatty acids, bug FAs, and seston FAs - script looks at seasonal trends in PUFAs in the seston, in the bugs and in birds.  Also does NMDS to see if bugs, primary production, and birds all have similar FA profiles
#PTK 24 May 2014

library(vegan)

#load all fatty acid data
#load TS chick FA data
setwd('~/bugEmergenceProject')
chick.fa<-read.csv('tsFAs.csv')

#load bolus FAs
bolus.fa<-read.csv('bolusFAs.csv')

#load bug FAs
bug.fa<-read.csv('bugFAs_ug_mg.csv')

#load seston FAs
ses.fa<-read.csv('SES63_FAs.csv')

#load chick growth data
chickGrowth<-read.csv('chickGrowth_25Nov2013.csv')

#pcoa with seston data
#This looks at how the FA profiles across all sites and time points are different from each other
ses<-ses.fa[,6:52]
ses<-as.matrix(ses)

ses.dist<-vegdist(ses,method='bray')
ses.pcoa<-cmdscale(ses.dist)

plot(ses.pcoa,pch=19,cex=0,xlab='PCA1',ylab='PCA2',col=c(rep('red',8),rep('blue',8),rep('green',8),rep('orange',8),rep('violet',8),rep('cyan',8)))
text(ses.pcoa[,1],ses.pcoa[,2],ses.fa[,2],col=c(rep('red',8),rep('blue',8),rep('green',8),rep('orange',8),rep('violet',8),rep('cyan',8))) #It looks like there is strong grouping by river (ILR v UMR v Em.)

anosim(ses.dist,grouping=ses.fa[,2]) #sig. = 0.001

#Now try this with bird FA data
chickFA<-chick.fa[,as.numeric(c(17:63))]
chickFA<-as.matrix(chickFA)

chick.dist<-vegdist(chickFA,method='bray')
chick.pcoa<-cmdscale(chick.dist)

plot(chick.pcoa,cex=0,xlab='PCA1',ylab='PCA2') #No real grouping by site
text(chick.pcoa[,1],chick.pcoa[,2],chick.fa$Site_Abbrev)

summary(anosim(chick.dist,grouping=chick.fa$Site_Abbrev)) #significance = 0.001

#PCA with bolus data
bolusFA<-as.matrix(bolus.fa[,6:48])

bolus.dist<-vegdist(bolusFA)
bolus.pcoa<-cmdscale(bolus.dist)

plot(bolus.pcoa,cex=0)
text(bolus.pcoa[,1],bolus.pcoa[,2],bolus.fa$Site_Abbrev)

#Compare chick growth to fatty acids in bolus as well as in livers and seston FA concentration
#calculate averages for pufa, omega 3, epa, dha
#First, calculate average chick growth by site
chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,mean,na.rm=T)[-1] #calculate mean chick growth

#calculate average PUFA, epa, w3, and dha for bolus
avgBolus.pufa<-tapply(bolus.fa$PUFA,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.pufa<-tapply(bolus.fa$PUFA,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.epa<-tapply(bolus.fa$C22_5n3c,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.epa<-tapply(bolus.fa$C22_5n3c,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.w3<-tapply(bolus.fa$Omega3,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.w3<-tapply(bolus.fa$Omega3,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.dha<-tapply(bolus.fa$C22_6n3,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.dha<-tapply(bolus.fa$C22_6n3,bolus.fa$Site_Abbrev,std.error,na.rm=T)
plot(avgBolus.pufa,chick.growth)
summary(lm(chick.growth~avgBolus.pufa)) #r2=0.89 p=0.007

#do the same for TS chicks
avgChick.pufa<-tapply(chick.fa$PUFA,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.pufa<-tapply(chick.fa$PUFA,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.epa<-tapply(chick.fa$C22_5n3c,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.epa<-tapply(chick.fa$C22_5n3c,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.w3<-tapply(chick.fa$Omega3,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.w3<-tapply(chick.fa$Omega3,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.dha<-tapply(chick.fa$C22_6n3,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.dha<-tapply(chick.fa$C22_6n3,chick.fa$Site_Abbrev,std.error,na.rm=T)
plot(avgChick.pufa,chick.growth) #no relationship

#do the same for seston
avgSes.pufa<-tapply(ses.fa$PUFA,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.pufa<-tapply(ses.fa$PUFA,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.epa<-tapply(ses.fa$C22_5n3c,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.epa<-tapply(ses.fa$C22_5n3c,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.w3<-tapply(ses.fa$Omega3,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.w3<-tapply(ses.fa$Omega3,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.dha<-tapply(ses.fa$C22_6n3,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.dha<-tapply(ses.fa$C22_6n3,ses.fa$Site_Abbrev,std.error,na.rm=T)
plot(avgSes.pufa,chick.growth)
summary(lm(chick.growth~avgSes.pufa)) #r2=0.43 p=0.154 slope = -0.0088

#Story so far -> chick growth influenced strongly by PUFA availability, but PUFA availability in aquatic primary production or seston does not seem to be a driver of that availability in forgaing.  PUFAs in the diet moreso a function of the taxonomic composition of the prey, and the availaility of "high PUFA" taxa compared to lower PUFA taxa (i.e. mainly brachyceran diptera) - therefore fatty acid characteristics of a subsidy appear to be vitally important for tree swallow nestling growth.

#makes bar graphs of average PUFA, EPA, DHA, ALA, ARA, and LIN for seston, boli, and birds

#Barplot of comparisons among seston, bolus, and bird fatty acids
#start with PUFAs
pufa.table<-data.frame(avgSes.pufa,avgBolus.pufa,avgChick.pufa)
pufa.table<-t(pufa.table)
se.pufa.table<-data.frame(seSes.pufa,seBolus.pufa,seChick.pufa)
se.pufa.table<-t(se.pufa.table)
w3.table<-data.frame(avgSes.w3,avgBolus.w3,avgChick.w3)
w3.table<-t(w3.table)
se.w3.table<-data.frame(seSes.w3,seBolus.w3,seChick.w3)
se.w3.table<-t(se.w3.table)
epa.table<-data.frame(avgSes.epa,avgBolus.epa,avgChick.epa)
epa.table<-t(epa.table)
se.epa.table<-data.frame(seSes.epa,seBolus.epa,seChick.epa)
se.epa.table<-t(se.epa.table)
dha.table<-data.frame(avgSes.dha,avgBolus.dha,avgChick.dha)
dha.table<-t(dha.table)
se.dha.table<-data.frame(seSes.dha,seBolus.dha,seChick.dha)
se.dha.table<-t(se.dha.table)

par(mfrow=c(2,2))
table=epa.table
se.table=se.epa.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main='EPA') 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0)

table<-dha.table
se.table<-se.dha.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main='DHA') 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0)

table<-w3.table
se.table<-se.w3.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main=expression(paste(omega,'3'))) 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0)

barplot(pufa.table,beside=T,ylim=c(0,max(se.pufa.table+pufa.table)),main='PUFA') 
error.bar(1.5,pufa.table[1,1],upper=se.pufa.table[1,1],lower=0)
error.bar(2.5,pufa.table[2,1],upper=se.pufa.table[2,1],lower=0)
error.bar(3.5,pufa.table[3,1],upper=se.pufa.table[3,1],lower=0)
error.bar(5.5,pufa.table[1,2],upper=se.pufa.table[1,2],lower=0)
error.bar(6.5,pufa.table[2,2],upper=se.pufa.table[2,2],lower=0)
error.bar(7.5,pufa.table[3,2],upper=se.pufa.table[3,2],lower=0)
error.bar(9.5,pufa.table[1,3],upper=se.pufa.table[1,3],lower=0)
error.bar(10.5,pufa.table[2,3],upper=se.pufa.table[2,3],lower=0)
error.bar(11.5,pufa.table[3,3],upper=se.pufa.table[3,3],lower=0)
error.bar(13.5,pufa.table[1,4],upper=se.pufa.table[1,4],lower=0)
error.bar(14.5,pufa.table[2,4],upper=se.pufa.table[2,4],lower=0)
error.bar(15.5,pufa.table[3,4],upper=se.pufa.table[3,4],lower=0)
error.bar(17.5,pufa.table[1,5],upper=se.pufa.table[1,5],lower=0)
error.bar(18.5,pufa.table[2,5],upper=se.pufa.table[2,5],lower=0)
error.bar(19.5,pufa.table[3,5],upper=se.pufa.table[3,5],lower=0)
error.bar(21.5,pufa.table[1,6],upper=se.pufa.table[1,6],lower=0)
error.bar(22.5,pufa.table[2,6],upper=se.pufa.table[2,6],lower=0)
error.bar(23.5,pufa.table[3,6],upper=se.pufa.table[3,6],lower=0)


#load water quality data to compare primary production to FA availability in chicks/bolus/bugs, but with bugs do multiple regression with taxa and time as variables
#load water quality data
wq<-read.csv('SPPPWQ.csv')
wq$Date<-sub('/14','/10',wq$Date)

chl<-c() #add water chemistry data to bug.fa dat to look at influence of limnological conditions on fatty acids
tpk<-c()
tnk<-c()
tss<-c()
vss<-c()
for(i in 1:nrow(bug.fa)){
	samplei=bug.fa[i,]
	wqi=wq[wq$Site_Abrreviation==samplei$Site_Abbrev,]
	wqi=wqi[!is.na(wqi$tpk),]
	x=wqi[which(min(abs(as.Date(samplei$Collection_date,'%m/%d/%y')-as.Date(wqi$Date,'%m/%d/%y')))==abs(as.Date(samplei$Collection_date,'%m/%d/%y')-as.Date(wqi$Date,'%m/%d/%y'))),]
	chl[i]=wqi$chlorophyll
	tpk[i]=wqi$tpk
	tnk[i]=wqi$tnk
	tss[i]=wqi$tss
	vss[i]=wqi$vss
}
bug.fa$chl=chl
bug.fa$TP=tpk
bug.fa$TN=tnk
bug.fa$tss=tss
bug.fa$vss=vss

plot(bug.fa$chl,bug.fa$PUFA)
summary(lm(bug.fa$PUFA~bug.fa$chl)) #r2=0.06 p=0.003 slope=-0.08
summary(lm(bug.fa$PUFA~bug.fa$chl+as.factor(bug.fa$Family))) #r2=0.39 p=0.0000007


