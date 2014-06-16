#his script looks at other potential predictors of tree swallow growth including size of bolus, avg PUFA content of the bugs from each site, and site productivity.  
#PTK 14 June 2014

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

#load water chem data
wq<-read.csv('SPPPWQ.csv')

#load chick growth data
chickGrowth<-read.csv('chickGrowth_25Nov2013.csv')

#calculate average and SE for chick growth as well as EFA concentration for regressions
chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,mean,na.rm=T)[-1] #calculate mean chick growth
se.chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,std.error,na.rm=T)[-1]

#load bolus data
setwd('~/bugEmergenceProject')

bolus.data<-read.csv('bolusData_25Nov2013.csv')
#make uniqueID
bolus.data$bolusID<-paste(bolus.data$Box,bolus.data$Bolus,sep='.')

bolus.counts<-aggregate(as.numeric(bolus.data$Count),by=list(bolus.data$bolusID,bolus.data$Date,bolus.data$Site),FUN=sum,na.rm=T)
colnames(bolus.counts)<-c('bolusID','Date','Site','count')

avgCounts<-tapply(bolus.counts$count,bolus.counts$Site,mean,na.rm=T)[-1]

plot(avgCounts,chick.growth)
summary(lm(chick.growth~avgCounts)) #r2=0.013, p=0.83

#<---------Cant calculate biomass data for bolus because we do not have biomass data for all taxa that are represented in the diet----------->

#look at tree swallow growth vs site productivity
avgChl<-tapply(wq$chlorophyll,wq$Site_Abrreviation,mean,na.rm=T)
plot(avgChl,chick.growth) #no relationship

avgTP<-tapply(wq$tpk,wq$Site_Abrreviation,mean,na.rm=T)
plot(avgTP,chick.growth) #no relationship

#tree swallow growth vs seston fatty acids
avgSes.pufa<-tapply(ses.fa$PUFA,ses.fa$Site_Abbrev,mean,na.rm=T)
plot(avgSes.pufa,chick.growth)
summary(lm(chick.growth~avgSes.pufa)) #r2=0.44 p=0.15 slope=-.0088

#tree swallow growth vs %Brachycera
#source script that calculates % diet for each taxon
#use dietTable
source('TSdiet_25March2014.R')

plot(as.numeric(dietTable[3,2:7]),chick.growth)
summary(lm(chick.growth~as.numeric(dietTable[3,2:7]))) #r2=0.64 p=0.06 slope=-0.005364

#make plots of growth vs all taxa, see if there are any other paterns somewhere
for(i in 1:nrow(dietTable)){
	quartz()
	plot(as.numeric(dietTable[i,2:7]),chick.growth,pch=19,cex=1.2,main=dietTable[i,1])
} #no other evident patterns in diet composition in terms of % a certain taxa and growth

