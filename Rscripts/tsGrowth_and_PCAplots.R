#This script looks at the differences in tree swallow chick fatty acids across sites and makes comparisons to bolus fatty acids, bug FAs, and seston FAs - script looks at seasonal trends in PUFAs in the seston, in the bugs and in birds.  Also does PCA to see if bugs, primary production, and birds all have similar FA profiles
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

#load water quality data
wqData<-read.csv('SPPPWQ.csv')

#pcoa with seston data
#This looks at how the FA profiles across all sites and time points are different from each other
ses<-ses.fa[,6:52]
ses<-as.matrix(ses)

ses.dist<-vegdist(ses,method='bray')
ses.pcoa<-cmdscale(ses.dist)

anosim(ses.dist,grouping=ses.fa[,2]) #sig. = 0.001

#make data frame of environmental variables and specific fatty acids to run envfit
#start with specific fatty acids.  Just use the sites and the columns of fatty acids
fa<-data.frame(ses.fa[,10:(ncol(ses.fa)-6)])
fit<-envfit(ses.dist,fa)

#make data frame of environmental variables using the mean of wq data
TP<-tapply(wqData$tpk,wqData$Site_Abrreviation,mean,na.rm=T)
SRP<-tapply(wqData$srpk,wqData$Site_Abrreviation,mean,na.rm=T)
TN<-tapply(wqData$tnk,wqData$Site_Abrreviation,mean,na.rm=T)
NO3<-tapply(wqData$nox_nc,wqData$Site_Abrreviation,mean,na.rm=T)
NH4<-tapply(wqData$nhx_n,wqData$Site_Abrreviation,mean,na.rm=T)
tss<-tapply(wqData$tss,wqData$Site_Abrreviation,mean,na.rm=T)
vss<-tapply(wqData$vss,wqData$Site_Abrreviation,mean,na.rm=T)
chl<-tapply(wqData$chlorophyll,wqData$Site_Abrreviation,mean,na.rm=T)
wqTable<-data.frame(site=c('BN','BS','EN','ES','LL','SH'),TP,SRP,TN,NO3,NH4,tss,vss,chl)


#Make a table of these data to match to and make enviro. variable data frame
envTable<-c()
for(i in 1:nrow(ses.fa)){
	rowi=match(ses.fa$Site_Abbrev[i],wqTable$site)
	tp=wqTable$TP[rowi]
	srp=wqTable$SRP[rowi]
	tn=wqTable$TN[rowi]
	no3=wqTable$NO3[rowi]
	nh4=wqTable$NH4[rowi]
	tss=wqTable$tss[rowi]
	vss=wqTable$vss[rowi]
	chl=wqTable$chl[rowi]
	x=data.frame(tp,srp,tn,no3,nh4,tss,vss,chl)
	envTable=rbind(x,envTable)
}

fit.wq<-envfit(ses.dist,envTable)

plot(ses.pcoa,pch=19,cex=0,xlab='PCA1',ylab='PCA2',col=c(rep('red',8),rep('blue',8),rep('green',8),rep('orange',8),rep('violet',8),rep('cyan',8)))
text(ses.pcoa[,1],ses.pcoa[,2],ses.fa[,2],col=c(rep('red',8),rep('blue',8),rep('green',8),rep('orange',8),rep('violet',8),rep('cyan',8))) #It looks like there is strong grouping by river (ILR v UMR v Em.)
plot(fit,p.max=0.05,col='black',cex=0.5)
plot(fit.wq,p.max=0.05,col='red',cex=0.5)

#Now try this with bird FA data
chickFA<-chick.fa[,as.numeric(c(17:63))]
chickFA<-as.matrix(chickFA)

chick.dist<-vegdist(chickFA,method='bray')
chick.pcoa<-cmdscale(chick.dist)

chickfa<-data.frame(chick.fa[,c(6,10,21,23:27,30:32,35:38,40:49,51:54,57:60,62:63)])
chickfa.fit<-envfit(chick.dist,chickfa)

envTable<-c()
for(i in 1:nrow(chick.fa)){
	rowi=match(chick.fa$Site_Abbrev[i],wqTable$site)
	tp=wqTable$TP[rowi]
	srp=wqTable$SRP[rowi]
	tn=wqTable$TN[rowi]
	no3=wqTable$NO3[rowi]
	nh4=wqTable$NH4[rowi]
	tss=wqTable$tss[rowi]
	vss=wqTable$vss[rowi]
	chl=wqTable$chl[rowi]
	x=data.frame(tp,srp,tn,no3,nh4,tss,vss,chl)
	envTable=rbind(x,envTable)
}

fit.wq<-envfit(chick.dist,envTable)

plot(chick.pcoa,cex=0,xlab='PCA1',ylab='PCA2') #No real grouping by site
text(chick.pcoa[,1],chick.pcoa[,2],chick.fa$Site_Abbrev,col=c(rep('red',5),rep('blue',15),rep('green',12),rep('orange',10),rep('violet',10),rep('cyan',10)))
#plot(chickfa.fit,p.max=0.05,col='black')
plot(chickfa.fit,p.max=0.05,col='black')
summary(anosim(chick.dist,grouping=chick.fa$Site_Abbrev)) #significance = 0.001

#PCA with bolus data
bolusFA<-as.matrix(bolus.fa[,c(6:18,20:23,25:28,30,32:39,42:45,47:48)])

bolus.dist<-vegdist(bolusFA)
bolus.pcoa<-cmdscale(bolus.dist)
bolusFA.fit<-envfit(bolus.dist,bolusFA)

envTable<-c()
for(i in 1:nrow(bolus.fa)){
	rowi=match(bolus.fa$Site_Abbrev[i],wqTable$site)
	tp=wqTable$TP[rowi]
	srp=wqTable$SRP[rowi]
	tn=wqTable$TN[rowi]
	no3=wqTable$NO3[rowi]
	nh4=wqTable$NH4[rowi]
	tss=wqTable$tss[rowi]
	vss=wqTable$vss[rowi]
	chl=wqTable$chl[rowi]
	x=data.frame(tp,srp,tn,no3,nh4,tss,vss,chl)
	envTable=rbind(x,envTable)
}
bolusFA.wq.fit<-envfit(bolus.dist,envTable)

plot(bolus.pcoa,cex=0)
text(bolus.pcoa[,1],bolus.pcoa[,2],bolus.fa$Site_Abbrev,col=c(rep('red',nrow(bolus.fa[bolus.fa$Site_Abbrev=='BN',])),rep('blue',nrow(bolus.fa[bolus.fa$Site_Abbrev=='BS',])),rep('green',nrow(bolus.fa[bolus.fa$Site_Abbrev=='EN',])),rep('orange',nrow(bolus.fa[bolus.fa$Site_Abbrev=='ES',])),rep('violet',nrow(bolus.fa[bolus.fa$Site_Abbrev=='LL',])),rep('cyan',nrow(bolus.fa[bolus.fa$Site_Abbrev=='SH',]))))
plot(bolusFA.fit,p.max=0.05,col='black')
plot(bolusFA.wq.fit,p.max=0.05,col='red')

#Compare chick growth to fatty acids in bolus as well as in livers and seston FA concentration
#calculate averages for pufa, omega 3, epa, dha
#First, calculate average chick growth by site
chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,mean,na.rm=T)[-1] #calculate mean chick growth
se.chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,std.error,na.rm=T)[-1]

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
plot(avgSes.pufa,chick.growth,pch=19,ylim=c(min(chick.growth-se.chick.growth),max(chick.growth+se.chick.growth)))
errbar(avgSes.pufa,chick.growth,yplus=chick.growth+se.chick.growth,yminus=chick.growth-se.chick.growth,add=T)
summary(lm(chick.growth~avgSes.pufa)) #r2=0.43 p=0.154 slope = -0.0088


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


#calculate all of the above in % total fatty acids by dividing concentration by total fatty acids column - make a new data frame for each
#seston
ses.ptfa<-c()
for(i in 6:57){
	ses.ptfa<-cbind(ses.ptfa,(ses.fa[,i]/ses.fa[,58])*100)
}
colnames(ses.ptfa)=colnames(ses.fa[6:57])
ses.ptfa<-cbind(ses.fa[,1:5],ses.ptfa)
#bugs
bug.ptfa<-c()
for(i in 5:29){
	bug.ptfa<-cbind(bug.ptfa,(bug.fa[,i]/bug.fa[,30])*100)
}
colnames(bug.ptfa)<-colnames(bug.fa[5:29])
bug.ptfa<-cbind(bug.fa[,1:4],bug.ptfa)
#bolus
bolus.ptfa<-c()
for(i in 6:53){
	bolus.ptfa<-cbind(bolus.ptfa,(bolus.fa[,i]/bolus.fa[,54])*100)
}
colnames(bolus.ptfa)<-colnames(bolus.fa[6:53])
bolus.ptfa<-cbind(bolus.fa[,1:5],bolus.ptfa)
#TS chicks
chick.ptfa<-c()
for(i in 17:68){
	chick.ptfa<-cbind(chick.ptfa,(chick.fa[,i]/chick.fa[,69])*100)
}
colnames(chick.ptfa)<-colnames(chick.fa[17:68])
chick.ptfa<-cbind(chick.fa[,1:17],chick.ptfa)

#Calculate averages similar to above
#calculate average PUFA, epa, w3, and dha for bolus
avgBolus.pufa<-tapply(bolus.ptfa$PUFA,bolus.ptfa$Site_Abbrev,mean,na.rm=T)
seBolus.pufa<-tapply(bolus.ptfa$PUFA,bolus.ptfa$Site_Abbrev,std.error,na.rm=T)
avgBolus.epa<-tapply(bolus.ptfa$C22_5n3c,bolus.ptfa$Site_Abbrev,mean,na.rm=T)
seBolus.epa<-tapply(bolus.ptfa$C22_5n3c,bolus.ptfa$Site_Abbrev,std.error,na.rm=T)
avgBolus.w3<-tapply(bolus.ptfa$Omega3,bolus.ptfa$Site_Abbrev,mean,na.rm=T)
seBolus.w3<-tapply(bolus.ptfa$Omega3,bolus.ptfa$Site_Abbrev,std.error,na.rm=T)
avgBolus.dha<-tapply(bolus.ptfa$C22_6n3,bolus.ptfa$Site_Abbrev,mean,na.rm=T)
seBolus.dha<-tapply(bolus.ptfa$C22_6n3,bolus.ptfa$Site_Abbrev,std.error,na.rm=T)

#do the same for TS chicks
avgChick.pufa<-tapply(chick.ptfa$PUFA,chick.ptfa$Site_Abbrev,mean,na.rm=T)
seChick.pufa<-tapply(chick.ptfa$PUFA,chick.ptfa$Site_Abbrev,std.error,na.rm=T)
avgChick.epa<-tapply(chick.ptfa$C22_5n3c,chick.ptfa$Site_Abbrev,mean,na.rm=T)
seChick.epa<-tapply(chick.ptfa$C22_5n3c,chick.ptfa$Site_Abbrev,std.error,na.rm=T)
avgChick.w3<-tapply(chick.ptfa$Omega3,chick.ptfa$Site_Abbrev,mean,na.rm=T)
seChick.w3<-tapply(chick.ptfa$Omega3,chick.ptfa$Site_Abbrev,std.error,na.rm=T)
avgChick.dha<-tapply(chick.ptfa$C22_6n3,chick.ptfa$Site_Abbrev,mean,na.rm=T)
seChick.dha<-tapply(chick.ptfa$C22_6n3,chick.ptfa$Site_Abbrev,std.error,na.rm=T)

#do the same for seston
avgSes.pufa<-tapply(ses.ptfa$PUFA,ses.ptfa$Site_Abbrev,mean,na.rm=T)
seSes.pufa<-tapply(ses.ptfa$PUFA,ses.ptfa$Site_Abbrev,std.error,na.rm=T)
avgSes.epa<-tapply(ses.ptfa$C22_5n3c,ses.ptfa$Site_Abbrev,mean,na.rm=T)
seSes.epa<-tapply(ses.ptfa$C22_5n3c,ses.ptfa$Site_Abbrev,std.error,na.rm=T)
avgSes.w3<-tapply(ses.ptfa$Omega3,ses.ptfa$Site_Abbrev,mean,na.rm=T)
seSes.w3<-tapply(ses.ptfa$Omega3,ses.ptfa$Site_Abbrev,std.error,na.rm=T)
avgSes.dha<-tapply(ses.ptfa$C22_6n3,ses.ptfa$Site_Abbrev,mean,na.rm=T)
seSes.dha<-tapply(ses.ptfa$C22_6n3,ses.ptfa$Site_Abbrev,std.error,na.rm=T)

#Also do it for all bugs, irregardless of taxa
avgBug.pufa<-tapply(bug.ptfa$PUFA,bug.ptfa$Site_Abbrev,mean,na.rm=T)
seBug.pufa<-tapply(bug.ptfa$PUFA,bug.ptfa$Site_Abbrev,std.error,na.rm=T)
avgBug.epa<-tapply(bug.ptfa$C22_5n3c,bug.ptfa$Site_Abbrev,mean,na.rm=T)
seBug.epa<-tapply(bug.ptfa$C22_5n3c,bug.ptfa$Site_Abbrev,std.error,na.rm=T)
avgBug.dha<-tapply(bug.ptfa$DHAC22_6n3,bug.ptfa$Site_Abbrev,mean,na.rm=T)
seBug.dha<-tapply(bug.ptfa$DHAC22_6n3,bug.ptfa$Site_Abbrev,std.error,na.rm=T)
avgBug.w3<-tapply(bug.ptfa$Omega3,bug.ptfa$Site_Abbrev,mean,na.rm=T)
seBug.w3<-tapply(bug.ptfa$Omega3,bug.ptfa$Site_Abbrev,std.error,na.rm=T)

#Barplot of comparisons among seston, bolus, and bird fatty acids
#start with PUFAs ----- this time with percent fatty acids
pufa.table<-data.frame(avgSes.pufa,avgBug.pufa,avgBolus.pufa,avgChick.pufa)
pufa.table<-t(pufa.table)
se.pufa.table<-data.frame(seSes.pufa,seBug.pufa,seBolus.pufa,seChick.pufa)
se.pufa.table<-t(se.pufa.table)
w3.table<-data.frame(avgSes.w3,avgBug.w3,avgBolus.w3,avgChick.w3)
w3.table<-t(w3.table)
se.w3.table<-data.frame(seSes.w3,seBug.w3,seBolus.w3,seChick.w3)
se.w3.table<-t(se.w3.table)
epa.table<-data.frame(avgSes.epa,avgBug.epa,avgBolus.epa,avgChick.epa)
epa.table<-t(epa.table)
se.epa.table<-data.frame(seSes.epa,seBug.epa,seBolus.epa,seChick.epa)
se.epa.table<-t(se.epa.table)
dha.table<-data.frame(avgSes.dha,avgBug.dha,avgBolus.dha,avgChick.dha)
dha.table<-t(dha.table)
se.dha.table<-data.frame(seSes.dha,seBug.dha,seBolus.dha,seChick.dha)
se.dha.table<-t(se.dha.table)


####<------------------------------------------------------>####

#Look at transfer of bacterial fatty acids seston, bolus, and birds
#Seston first
avgSes.c15i<-tapply(ses.fa$C15_0i,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.c15i<-tapply(ses.fa$C15_0i,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.c15ai<-tapply(ses.fa$C15ai,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.c15ai<-tapply(ses.fa$C15ai,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.c15<-tapply(ses.fa$C15_0,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.c15<-tapply(ses.fa$C15_0,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.c151<-tapply(ses.fa$C15_1,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.c151<-tapply(ses.fa$C15_1,ses.fa$Site_Abbrev,std.error,na.rm=T)
avgSes.c17<-tapply(ses.fa$C17_0,ses.fa$Site_Abbrev,mean,na.rm=T)
seSes.c17<-tapply(ses.fa$C17_0,ses.fa$Site_Abbrev,std.error,na.rm=T)

#Bolus
avgBolus.c15i<-tapply(bolus.fa$C15_0i,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.c15i<-tapply(bolus.fa$C15_0i,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.c15ai<-tapply(bolus.fa$C15ai,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.c15ai<-tapply(bolus.fa$C15ai,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.c15<-tapply(bolus.fa$C15_0,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.c15<-tapply(bolus.fa$C15_0,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.c151<-tapply(bolus.fa$C15_1,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.c151<-tapply(bolus.fa$C15_1,bolus.fa$Site_Abbrev,std.error,na.rm=T)
avgBolus.c17<-tapply(bolus.fa$C17_0,bolus.fa$Site_Abbrev,mean,na.rm=T)
seBolus.c17<-tapply(bolus.fa$C17_0,bolus.fa$Site_Abbrev,std.error,na.rm=T)

#Birds
avgChick.c15i<-tapply(chick.fa$C15_0i,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.c15i<-tapply(chick.fa$C15_0i,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.c15ai<-tapply(chick.fa$C15ai,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.c15ai<-tapply(chick.fa$C15ai,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.c15<-tapply(chick.fa$C15_0,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.c15<-tapply(chick.fa$C15_0,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.c151<-tapply(chick.fa$C15_1,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.c151<-tapply(chick.fa$C15_1,chick.fa$Site_Abbrev,std.error,na.rm=T)
avgChick.c17<-tapply(chick.fa$C17_0,chick.fa$Site_Abbrev,mean,na.rm=T)
seChick.c17<-tapply(chick.fa$C17_0,chick.fa$Site_Abbrev,std.error,na.rm=T)

c15i.table<-data.frame(avgSes.c15i,avgBolus.c15i,avgChick.c15i)
c15i.table<-t(c15i.table)
se.c15i.table<-data.frame(seSes.c15i,seBolus.c15i,seChick.c15i)
se.c15i.table<-t(se.c15i.table)
c15ai.table<-data.frame(avgSes.c15ai,avgBolus.c15ai,avgChick.c15ai)
c15ai.table<-t(c15ai.table)
se.c15ai.table<-data.frame(seSes.c15ai,seBolus.c15ai,seChick.c15ai)
se.c15ai.table<-t(se.c15ai.table)
c15.table<-data.frame(avgSes.c15,avgBolus.c15,avgChick.c15)
c15.table<-t(c15.table)
se.c15.table<-data.frame(seSes.c15,seBolus.c15,seChick.c15)
se.c15.table<-t(se.c15.table)
c151.table<-data.frame(avgSes.c151,avgBolus.c151,avgChick.c151)
c151.table<-t(c151.table)
se.c151.table<-data.frame(seSes.c151,seBolus.c151,seChick.c151)
se.c151.table<-t(se.c151.table)
c17.table<-data.frame(avgSes.c17,avgBolus.c17,avgChick.c17)
c17.table<-t(c17.table)
se.c17.table<-data.frame(seSes.c17,seBolus.c17,seChick.c17)
se.c17.table<-t(se.c17.table)