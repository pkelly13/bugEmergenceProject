#Script that will make barplots of trophic transfer of bacterial fatty acids across sites.  Uses seston, bolus and birds, no data for bugs.
#PTK 12 June 2014

#set working directory
setwd('~/bugEmergenceProject')
#source script
source('tsGrowth_and_PCAplots.R')

#save plot as a pdf to potential figures folder
setwd('~/bugEmergenceProject/potential figures')

pdf('BACfattyAcidTrophicTransfer_barplot.pdf',width=11,height=8.5)
par(mfrow=c(2,2))
table=c15i.table
se.table=se.c15i.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main='iC15:0') 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0,length=0.04)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0,length=0.04)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0,length=0.04)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0,length=0.04)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0,length=0.04)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0,length=0.04)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0,length=0.04)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0,length=0.04)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0,length=0.04)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0,length=0.04)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0,length=0.04)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0,length=0.04)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0,length=0.04)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0,length=0.04)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0,length=0.04)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0,length=0.04)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0,length=0.04)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0,length=0.04)



table<-c15ai.table
se.table<-se.c15ai.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main='aiC15:0') 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0,length=0.04)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0,length=0.04)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0,length=0.04)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0,length=0.04)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0,length=0.04)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0,length=0.04)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0,length=0.04)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0,length=0.04)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0,length=0.04)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0,length=0.04)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0,length=0.04)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0,length=0.04)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0,length=0.04)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0,length=0.04)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0,length=0.04)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0,length=0.04)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0,length=0.04)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0,length=0.04)


table<-c15.table
se.table<-se.c15.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main='C15:0',legend=c('Seston','Bolus','Nestlings'),args.legend=list(x=10.5,y=1.5,cex=1,bty='n')) 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0,length=0.04)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0,length=0.04)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0,length=0.04)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0,length=0.04)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0,length=0.04)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0,length=0.04)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0,length=0.04)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0,length=0.04)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0,length=0.04)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0,length=0.04)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0,length=0.04)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0,length=0.04)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0,length=0.04)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0,length=0.04)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0,length=0.04)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0,length=0.04)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0,length=0.04)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0,length=0.04)


table=c17.table
se.table=se.c17.table
barplot(table,beside=T,ylim=c(0,max(se.table+table)),main='C17:0') 
error.bar(1.5,table[1,1],upper=se.table[1,1],lower=0,length=0.04)
error.bar(2.5,table[2,1],upper=se.table[2,1],lower=0,length=0.04)
error.bar(3.5,table[3,1],upper=se.table[3,1],lower=0,length=0.04)
error.bar(5.5,table[1,2],upper=se.table[1,2],lower=0,length=0.04)
error.bar(6.5,table[2,2],upper=se.table[2,2],lower=0,length=0.04)
error.bar(7.5,table[3,2],upper=se.table[3,2],lower=0,length=0.04)
error.bar(9.5,table[1,3],upper=se.table[1,3],lower=0,length=0.04)
error.bar(10.5,table[2,3],upper=se.table[2,3],lower=0,length=0.04)
error.bar(11.5,table[3,3],upper=se.table[3,3],lower=0,length=0.04)
error.bar(13.5,table[1,4],upper=se.table[1,4],lower=0,length=0.04)
error.bar(14.5,table[2,4],upper=se.table[2,4],lower=0,length=0.04)
error.bar(15.5,table[3,4],upper=se.table[3,4],lower=0,length=0.04)
error.bar(17.5,table[1,5],upper=se.table[1,5],lower=0,length=0.04)
error.bar(18.5,table[2,5],upper=se.table[2,5],lower=0,length=0.04)
error.bar(19.5,table[3,5],upper=se.table[3,5],lower=0,length=0.04)
error.bar(21.5,table[1,6],upper=se.table[1,6],lower=0,length=0.04)
error.bar(22.5,table[2,6],upper=se.table[2,6],lower=0,length=0.04)
error.bar(23.5,table[3,6],upper=se.table[3,6],lower=0,length=0.04)


mtext('Site',side=1,cex=1.4,line=-1.5,outer=T)
mtext('% of total FAs',side=2,cex=1.4,line=-1.5,outer=T)
dev.off()