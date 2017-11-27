# Analyses and data compilation associated with Maranger et al. Stoichiometry of carbon, nitrogen, and phospohorus through the freshwater pipe
# SEJ

rm(list=ls())

#*******************
# Table 2 - summary of stoichiometrires from U.S. EPA National Lakes Assessment 2007 and U.S. Geological Survey stream data via their Water Quality Portal
#*******************

# stream/river data; see below for code to pull data
stream=read.table("USGSstreamCPNQ_11-26-17.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE,colClasses="character")

stream$P_mgL=as.numeric(stream$P_mgL)
stream$C_mgL=as.numeric(stream$C_mgL)
stream$N_mgL=as.numeric(stream$N_mgL)
stream$Q_ft3S=as.numeric(stream$Q_ft3S)

stream=stream[!is.na(stream$C_mgL),] # one site had an NA for TOC
stream=stream[stream$P_mgL>0,] # 3 sites had zero for TP
stream=stream[stream$C_mgL>0,] # 13 sites had zero for TOC

# calculate stoichiometric ratios, all are in mg per L
streamCP=(stream$C_mgL/12)/(stream$P_mgL/31)
streamCN=(stream$C_mgL/12)/(stream$N_mgL/14)
streamNP=(stream$N_mgL/14)/(stream$P_mgL/31)

# look at range of discharges
Q=stream[!is.na(stream[,5]),5]  # remove sites with out discharge; only 1517 have a reported discharge
Q=Q[Q>0]  # 15 sites had a reported discharge of zero and another two reported a negative discharge
length(Q)
range(Q)*0.0283168  #m3 s-1

length(streamCN)

summary(streamCN)
summary(streamCP)
summary(streamNP)

# lake and reservoir data
lakechem=read.csv("NLA2007_WaterQuality_20091123.csv",header=TRUE,stringsAsFactors=FALSE)
lakechem=lakechem[!duplicated(lakechem$SITE_ID),]

# pulling in more lake information and harmonizing data
lakeinfo=read.csv("NLA2007_SampledLakeInformation_20091113.csv",header=TRUE,stringsAsFactors=FALSE)
lakeinfo=lakeinfo[!duplicated(lakeinfo$SITE_ID),]
lakeinfo=lakeinfo[lakeinfo$SITE_ID%in%lakechem$SITE_ID,]
lakeinfo=lakeinfo[order(lakeinfo$SITE_ID),]
lakechem=lakechem[order(lakechem$SITE_ID),]
dim(lakechem)
dim(lakeinfo)
sum(lakeinfo$SITE_ID==lakechem$SITE_ID)

# molar stoichiometries
lakeCP=(lakechem$TOC/12)/(lakechem$PTL/1000/31)
lakeCN=(lakechem$TOC/12)/(lakechem$NTL/1000/14)
lakeNP=(lakechem$NTL/1000/14)/(lakechem$PTL/1000/31)

# range of surface areas
range(lakeinfo$AREA_HA)*0.01  #km2

# reservoir stoichiometries
sum(lakeinfo$LAKE_ORIGIN=='MAN-MADE')
summary(lakeCN[lakeinfo$LAKE_ORIGIN=='MAN-MADE'])
summary(lakeCP[lakeinfo$LAKE_ORIGIN=='MAN-MADE'])
summary(lakeNP[lakeinfo$LAKE_ORIGIN=='MAN-MADE'])

# lake stoichiometries
sum(lakeinfo$LAKE_ORIGIN=='NATURAL')
summary(lakeCN[lakeinfo$LAKE_ORIGIN=='NATURAL'])
summary(lakeCP[lakeinfo$LAKE_ORIGIN=='NATURAL'])
summary(lakeNP[lakeinfo$LAKE_ORIGIN=='NATURAL'])


# code for pull from USGS Water Quality Portal - completed 2017-11-26
library(dataRetrieval)

# using data from 2007-current
startDate="2007-01-01"
endDate="2017-11-30"

states=c('AL','AK','AZ','CA','CO','CT','DE','DC','FL','GA','HI','ID','IL','IN','IA','KS','AR','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VT','VA','WA','WV','WI','WY')

out=matrix(c(0,0,0,0,0),1,5)
# warnings on 5, 23, 26, 44, but the errors seem to be on the server side and we get data from those states
for(i in (1:length(states))[-8]){ #District of Columbia (i=8) had no sites with TOC, TP, and TN data
  print(i)
  siteListP=whatNWISsites(parameterCd="00665",stateCd=states[i])
  siteListC=whatNWISsites(parameterCd="00680",stateCd=states[i])
  siteListN=whatNWISsites(parameterCd="00600",stateCd=states[i])
  
  streamSitesP=siteListP[siteListP[,4]=="ST",]
  streamSitesC=siteListC[siteListC[,4]=="ST",]
  streamSitesN=siteListN[siteListN[,4]=="ST",]
  
  pairedSites=streamSitesP$site_no[streamSitesP$site_no%in%streamSitesC$site_no]
  pairedSites=pairedSites[pairedSites%in%streamSitesN$site_no]
  
  Pdata=readNWISqw(pairedSites,parameterCd="00665") # mg/L
  Cdata=readNWISqw(pairedSites,parameterCd="00680")	# mg/L
  Ndata=readNWISqw(pairedSites,parameterCd="00600") # mg/L
  Qdata=readNWISqw(pairedSites,parameterCd="00060")	# ft3/s	
  
  Pmeans=tapply(Pdata$result_va,Pdata$site_no,FUN=mean,na.rm=TRUE)
  Cmeans=tapply(Cdata$result_va,Cdata$site_no,FUN=mean,na.rm=TRUE)
  Nmeans=tapply(Ndata$result_va,Ndata$site_no,FUN=mean,na.rm=TRUE)
  Qmeans=tapply(Qdata$result_va,Qdata$site_no,FUN=mean,na.rm=TRUE)
  
  pairedMeans=pairedSites[pairedSites%in%names(Pmeans)]
  pairedMeans=pairedMeans[pairedMeans%in%names(Cmeans)]
  pairedMeans=pairedMeans[pairedMeans%in%names(Nmeans)]
  
  Pmeans=Pmeans[names(Pmeans)%in%pairedMeans]
  Cmeans=Cmeans[names(Cmeans)%in%pairedMeans]
  Nmeans=Nmeans[names(Nmeans)%in%pairedMeans]
  Qmeans=Qmeans[names(Qmeans)%in%pairedMeans]
  
  PCNQ=cbind(pairedMeans,Pmeans[match(pairedMeans,names(Pmeans))],Cmeans[match(pairedMeans,names(Cmeans))],Nmeans[match(pairedMeans,names(Nmeans))],NA)
  if(length(Qmeans)>0){
    PCNQ[match(names(Qmeans),PCNQ[,1]),5]=Qmeans
  }
  out=rbind(out,PCNQ)
}

out=out[-1,]
colnames(out)=c('site','P_mgL','C_mgL','N_mgL','Q_ft3S')
#write.table(out,"USGSstreamCPNQ_11-26-17.txt",row.names=FALSE,sep="\t")

#*******************
# Figure 2 - kernel density plots of lake and reservoir stoichiometries from the U.S. EPA National Lakes Assessment 2007
#*******************

# uses data loaded above for summary table - lakechem and lakeinfo dataframes
# uses stoichiometry vectors from above too - CP, CN, NP

# log10 transform stoichiometries
logCP=log10(CP)
logCN=log10(CN)
logNP=log10(NP)

# REF vs NOT by man-made vs natural - kernel density plots
lakeClass=paste(lakeinfo$LAKE_ORIGIN,lakeinfo$REF_NUTR,sep="|")
table(lakeClass)

# C:P
mmrefCPdens=density(logCP[lakeClass=='MAN-MADE|Y'])
mmnrCPdens=density(logCP[lakeClass=='MAN-MADE|N'])
natrefCPdens=density(logCP[lakeClass=='NATURAL|Y'])
natnrCPdens=density(logCP[lakeClass=='NATURAL|N'])

plot(natrefCPdens,lwd=3,col='grey70',type='l',xlab="C:P",ylab="Density",xlim=c(-0.5,4.5),ylim=c(0,1),lty=3,main="",axes=FALSE)
lines(natnrCPdens,lwd=3,col='grey70')
lines(mmrefCPdens,lwd=3,col='black',lty=3)
lines(mmnrCPdens,lwd=3,col='black')
axis(side=1,at=0:4,labels=c("1","10","100","1000","10000"))
axis(side=2,at=c(0,0.5,1),labels=c("0.0","0.5","1.0"))
box()
legend('topleft',c('lake-reference','lake-impacted','reservoir-reference','reservoir-impacted'),lwd=3,col=c('grey70','grey70','black','black'),lty=c(3,1,3,1),box.lty=0)

# C:N
mmrefCNdens=density(logCN[lakeClass=='MAN-MADE|Y'])
mmnrCNdens=density(logCN[lakeClass=='MAN-MADE|N'])
natrefCNdens=density(logCN[lakeClass=='NATURAL|Y'])
natnrCNdens=density(logCN[lakeClass=='NATURAL|N'])

plot(natrefCNdens,lwd=3,col='grey70',type='l',main="",xlab="C:N",ylab="Density",xlim=c(-0.5,4.5),ylim=c(0,2),lty=3,axes=FALSE)
lines(natnrCNdens,lwd=3,col='grey70')
lines(mmrefCNdens,lwd=3,col='black',lty=3)
lines(mmnrCNdens,lwd=3,col='black')
axis(side=1,at=0:4,labels=c("1","10","100","1000","10000"))
axis(side=2,at=c(0,1,2),labels=c("0.0","1.0","2.0"))
box()

# N:P
mmrefNPdens=density(logNP[lakeClass=='MAN-MADE|Y'])
mmnrNPdens=density(logNP[lakeClass=='MAN-MADE|N'])
natrefNPdens=density(logNP[lakeClass=='NATURAL|Y'])
natnrNPdens=density(logNP[lakeClass=='NATURAL|N'])

plot(natrefNPdens,lwd=3,col='grey70',type='l',main="",xlab="N:P",ylab="Density",xlim=c(-0.5,4.5),ylim=c(0,1.5),lty=3,axes=FALSE)
lines(natnrNPdens,lwd=3,col='grey70')
lines(mmrefNPdens,lwd=3,col='black',lty=3)
lines(mmnrNPdens,lwd=3,col='black')
axis(side=1,at=0:4,labels=c("1","10","100","1000","10000"))
axis(side=2,at=c(0,0.5,1,1.5),labels=c("0.0","0.5","1.0","1.5"))
box()