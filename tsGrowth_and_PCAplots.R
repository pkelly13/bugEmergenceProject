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

plot(ses.pcoa,cex=0,xlab='PCA1',ylab='PCA2')
text(ses.pcoa[,1],ses.pcoa[,2],ses.fa[,2]) #It looks like there is strong grouping by river (ILR v UMR v Em.)

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
chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,mean,na.rm=T)[-1] #calculate mean chick growth

avgBolus.pufa<-tapply(bolus.fa$PUFA,bolus.fa$Site_Abbrev,mean,na.rm=T)
plot(avgBolus.pufa,chick.growth)
summary(lm(chick.growth~avgBolus.pufa)) #r2=0.89 p=0.007

avgChick.pufa<-tapply(chick.fa$PUFA,chick.fa$Site_Abbrev,mean,na.rm=T)
plot(avgChick.pufa,chick.growth) #no relationship

avgSes.pufa<-tapply(ses.fa$PUFA,ses.fa$Site_Abbrev,mean,na.rm=T)
plot(avgSes.pufa,chick.growth)
summary(lm(chick.growth~avgSes.pufa)) #r2=0.43 p=0.154 slope = -0.0088

#Story so far -> chick growth influenced strongly by PUFA availability, but PUFA availability in aquatic primary production or seston does not seem to be a driver of that availability in forgaing.  PUFAs in the diet moreso a function of the taxonomic composition of the prey, and the availaility of "high PUFA" taxa compared to lower PUFA taxa (i.e. mainly brachyceran diptera) - therefore fatty acid characteristics of a subsidy appear to be vitally important for tree swallow nestling growth.

#makes bar graphs of average PUFA, EPA, DHA, ALA, ARA, and LIN for seston, boli, and birds
pufa.table<-data.frame(avgSes.pufa,avgBolus.pufa,avgChick.pufa)
pufa.table<-t(pufa.table)

barplot(pufa.table,beside=T) #This should probably be in %total fatty acids

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


