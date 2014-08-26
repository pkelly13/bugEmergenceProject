#Analyze chick growth vs. specific FAs
#PTK 4 March 2014

#load chick growth data
setwd('~/Documents/Masters Thesis/Publication v2/data analyses/R files')
chick=read.csv('chickGrowth_25Nov2013.csv') #loads chick growth 

#load specific lipid data for bolus
bolus.fa<-read.csv('bolusFAs.csv')

#look for statistical difference among bolus FAs
#do in fixed effect model design
#need to add location to bolus data
river<-c()
for(i in 1:nrow(bolus.fa)){
	if(bolus.fa$Site_Abbrev[i]=='BN' | bolus.fa$Site_Abbrev[i]=='BS'){
		river[i]='ILR'
	}
	else if(bolus.fa$Site_Abbrev[i]=='EN' | bolus.fa$Site_Abbrev[i]=='ES'){
		river[i]='Em'
	}
	else{
		river[i]='UMR'
	}
}
bolus.fa$river=river


#calculate average chick growth
avg.growth<-tapply(chick$Growth..g.d.,chick$Site_Abreviation,mean,na.rm=T)[-1]
se.growth<-tapply(chick$Growth..g.d.,chick$Site_Abreviation,std.error,na.rm=T)

#check to see if chick growth is significantly different among sites
#remove BE
chick<-chick[chick$Site_Abreviation!='BE',]
#add River to chick data
river<-c()
for(i in 1:nrow(chick)){
	if(chick$Site_Abreviation[i]=='BN' | chick$Site_Abreviation[i]=='BS'){
		river[i]='ILR'
	}
	else if(chick$Site_Abreviation[i]=='EN' | chick$Site_Abreviation[i]=='ES'){
		river[i]='Em'
	}
	else{
		river[i]='UMR'
	}
}
chick$river=river
#add Jdate
chick$Jdate<-strptime(chick$Date,'%m/%d/%y')$yday+1
bolus.fa$Jdate<-strptime(as.Date(bolus.fa$Collection_date,'%m/%d/%y'),'%Y-%m-%d')$yday+1

#now do ANOVA on chick growth by site
growth.mod<-aov(chick$Growth..g.d.~chick$river+chick$river/chick$Site_Abreviation)
anova(growth.mod) #significant difference in chick growth among sites - may need to control for site
TukeyHSD(growth.mod)

#calculate avg bolus PUFA concentration
avg.pufa<-tapply(as.numeric(bolus.fa$PUFA),as.factor(bolus.fa$Site_Abbrev),mean,na.rm=T)

#make plot for paper/Bill's presentation
setwd('~/Documents/Masters thesis/publication v2/new figures')
pdf(file='chickGrowth_vs_PUFA.pdf')
par(mar=c(5,5,5,5))
plot(avg.pufa,avg.growth,cex=0,xlab=expression(paste('Bolus PUFA content (',mu,'g PUFA mg'^-1,')')),ylab=expression(paste('Tree swallow chick growth (g d'^-1,')')),cex.axis=1.2,cex.lab=1.5,ylim=c(min(avg.growth-se.growth[-1]),max(avg.growth+se.growth[-1])))
text(avg.pufa,avg.growth,labels=rownames(avg.pufa),cex=1.2)
errbar(avg.pufa,avg.growth,yplus=avg.growth+se.growth[-1],yminus=avg.growth-se.growth[-1],add=T,cex=0)
abline(lm(avg.growth~avg.pufa))
dev.off()
summary(lm(avg.growth~avg.pufa)) #r2=0.87 p=0.007

#calculate avg bolus EPA conc.
avg.epa<-tapply(bolus.fa$C20_5n3,bolus.fa$Site_Abbrev,mean,na.rm=T)
plot(avg.epa,avg.growth)
summary(lm(avg.growth~avg.epa)) #r2=0.80 p=0.016

#calculate avg bolus DHA conc.
avg.dha<-tapply(bolus.fa$C22_6n3,bolus.fa$Site_Abbrev,mean,na.rm=T)
plot(avg.dha,avg.growth)
summary(lm(avg.growth~avg.dha)) #r2=0.68 p=0.04

#calculate avg omega3 conc
avg.n3<-tapply(bolus.fa$Omega3,bolus.fa$Site_Abbrev,mean,na.rm=T)
plot(avg.n3,avg.growth)
summary(lm(avg.growth~avg.n3)) #r2=0.85 p<0.01

#add n3:n6 to bolus.fa data
bolus.fa$n3_n6<-bolus.fa$Omega3/bolus.fa$Omega6
#avg n3:n6
n3n6<-tapply(bolus.fa$n3_n6,bolus.fa$Site_Abbrev,mean,na.rm=T)
plot(n3n6,avg.growth)
summary(lm(avg.growth~n3n6)) #r2=0.82 p=0.01

#look at average chick growth by average size of bolus, both number of insects and mass
#load bolus data
bolus.data<-read.csv('bolusData_25Nov2013.csv')
#make uniqueID
bolus.data$bolusID<-paste(bolus.data$Box,bolus.data$Bolus,sep='.')

bolus.counts<-aggregate(as.numeric(bolus.data$Count),by=list(bolus.data$bolusID,bolus.data$Date,bolus.data$Site),FUN=sum,na.rm=T)
colnames(bolus.counts)<-c('bolusID','Date','Site','count')

avgCounts<-tapply(bolus.counts$count,bolus.counts$Site,mean,na.rm=T)[-1]

summary(lm(avg.growth~avgCounts)) #r2 = 0.01 p = 0.83