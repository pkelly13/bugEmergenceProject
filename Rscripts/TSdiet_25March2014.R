#comparison of diets for tree swallows

#load bolus data
setwd('~/Documents/Masters thesis/publication v2/data analyses/R files')
bolus<-read.csv('bolusData_25Nov2013.csv')

#get rid of non necessary information
bolus<-bolus[,-c(3,5,6,ncol(bolus))]

#aggregate bolus data, total number of insect collected for each site - use total number of insects for each site similar to what Quinney and Ankney 1985 did
tot.count<-aggregate(as.numeric(bolus$Count),by=list(bolus$Family,bolus$Order,bolus$Site),FUN=sum,na.rm=T)
colnames(tot.count)<-c('family','order','site','count')

#remove BE site
tot.count<-tot.count[tot.count$site!='BE',]

#calculate percent diet ---BY FAMILY--- for each site
sites<-unique(tot.count$site)
percentDiet<-c()
for(i in 1:length(sites)){
	sitei=tot.count[tot.count$site==sites[i],]
	sumi=sum(as.numeric(sitei$count),na.rm=T)
	x=c()
	for(j in 1:nrow(sitei)){
		x[j]=(as.numeric(sitei$count[j])/sumi)*100
	}
	y=cbind(sitei,percentDiet=x)
	percentDiet=rbind(percentDiet,y)
}

#BY ORDER
#CHANGE DIPTERA TO NEMATOCERA AND BRACHYCERA
diptera<-tot.count[tot.count$order=='Diptera',]
dipteraClass<-c()
for(i in 1:nrow(diptera)){
	if(diptera$family[i]=='Tipulidae' | diptera$family[i]=='Chironomidae'){
		dipteraClass[i]='Nematocera'
	}
	else{
		dipteraClass[i]='Brachycera'
	}
}
diptera$order=dipteraClass

tot.count=rbind(tot.count[tot.count$order!="Diptera",],diptera)
order.count<-aggregate(tot.count$count,by=list(tot.count$order,tot.count$site),FUN=sum,na.rm=T)
colnames(order.count)<-c('order','site','count')

#calculate percent diet
percentDiet<-c()
for(i in 1:length(sites)){
	sitei=order.count[order.count$site==sites[i],]
	sumi=sum(sitei$count,na.rm=T)
	x=c()
	for(j in 1:nrow(sitei)){
		x[j]=(sitei$count[j]/sumi)*100
	}
	y=cbind(sitei,percentDiet=x)
	percentDiet=rbind(percentDiet,y)
}
dietTable<-cast(percentDiet,order~site)