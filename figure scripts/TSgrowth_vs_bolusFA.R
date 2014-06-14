#Script will make figures of tree swallow growth regressed against indices of diet fatty acid quality
#PTK 12 June 2014

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

#calculate average and SE for chick growth as well as EFA concentration for regressions
chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,mean,na.rm=T)[-1] #calculate mean chick growth
se.chick.growth<-tapply(chickGrowth$Growth..g.d.,chickGrowth$Site_Abreviation,std.error,na.rm=T)[-1]

#calculate average PUFA, epa, w3, w3:w6,dha, and ala for bolus
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

#make regression plots of EFAs vs tree swallow growth
setwd('~/bugEmergenceProject/potential figures')
pdf('chickGrowth_vs_FAs.pdf',height=11,width=8.5)
par(mfrow=c(2,2))
plot(avgBolus.epa,chick.growth,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,ylab='',xlab='',ylim=c(min(chick.growth-se.chick.growth),max(chick.growth+se.chick.growth)),main='EPA')
errbar(avgBolus.epa,chick.growth,yplus=chick.growth+se.chick.growth,yminus=chick.growth-se.chick.growth,add=TRUE)

plot(avgBolus.dha,chick.growth,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,ylab='',xlab='',ylim=c(min(chick.growth-se.chick.growth),max(chick.growth+se.chick.growth)),main='DHA')
errbar(avgBolus.dha,chick.growth,yplus=chick.growth+se.chick.growth,yminus=chick.growth-se.chick.growth,add=TRUE)
abline(lm(chick.growth~avgBolus.dha),lty=2)

plot(avgBolus.w3,chick.growth,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,ylab='',xlab='',ylim=c(min(chick.growth-se.chick.growth),max(chick.growth+se.chick.growth)),main=expression(paste(omega,' 3')))
errbar(avgBolus.w3,chick.growth,yplus=chick.growth+se.chick.growth,yminus=chick.growth-se.chick.growth,add=TRUE)
abline(lm(chick.growth~avgBolus.w3),lty=2)


plot(avgBolus.pufa,chick.growth,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,ylab='',xlab='',ylim=c(min(chick.growth-se.chick.growth),max(chick.growth+se.chick.growth)),main='PUFA')
errbar(avgBolus.pufa,chick.growth,yplus=chick.growth+se.chick.growth,yminus=chick.growth-se.chick.growth,add=TRUE)
abline(lm(chick.growth~avgBolus.pufa),lty=2)

mtext(expression(paste('Nestling growth (g day'^-1,')')),2,line=-1.6,cex=1.2,outer=T)
mtext(expression(paste('Mean bolus FA (',mu,'g mg DW'^-1,')')),side=1,line=-1.5,cex=1.2,outer=T)
dev.off()

#make similar graph showing regressions vs average bug PUFA, seston PUFA, and bolus size (weight)

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

#make graph of %Brachyceran vs chick growth
#source 'tsGrowth_vs_nonFApredictors.R' - use dietTable
setwd('~/bugEmergenceProject')
source('tsGrowth_vs_nonFApredictors.R')

#make pdf and save to potential figures folder
setwd('~/bugEmergenceProject/potential figures')
pdf('TSgrowth_vs_ptBrachycera.pdf',height=11,width=8.5)
plot(as.numeric(dietTable[3,2:7]),chick.growth,pch=19,cex=1.2,cex.lab=1.2,cex.axis=1.2,xlab='% Brachycera',ylab=expression(paste('Nestling growth (g day'^-1,')')),ylim=c(min(chick.growth-se.chick.growth),max(chick.growth+se.chick.growth)))
errbar(as.numeric(dietTable[3,2:7]),chick.growth,yplus=chick.growth+se.chick.growth,yminus=chick.growth-se.chick.growth,add=T)
dev.off()