


setwd("C:/Users/nilu/Desktop/Transfer to 117/New folder")
invisible(utils::memory.limit(32000))



library (tidyr)
library(plyr)
library(stringr)
library(dplyr)
library(reshape2)
library(pivottabler)
library(data.table)
library(rgdal)
library(zoo)


library(rjags)
library(R2jags)
library(jagstools)
library(runjags)
library(mvc)

library(MCMCpack)
library(parallel)
library(mcmcplots)
library(superdiag)
library(abind)
library(snow)
library(coda)
library(R2WinBUGS)

library(maptools)
library(ggplot2)
library(corrplot)



#-------------- Model preparation and predictors--------------------- 


#--------------- Counts

df3=read.csv("data_Indextime_month2.csv",header=TRUE,stringsAsFactors = FALSE)


g2<- filter(df3,CZ_TYPE %in% c("C","Z"))

g2$CZ_FIPS[g2$CZ_FIPS==0] <- ""
g2$CZ_FIPS <- as.numeric(g2$CZ_FIPS)

g2$STATE_FIPS <- as.numeric(g2$STATE_FIPS)
g2$STATE_FIPS=formatC(g2$STATE_FIPS, width = 2, format = "d", flag = "0")

g2$CZ_FIPS=formatC(g2$CZ_FIPS, width = 3, format = "d", flag = "0") #injuri class NA ha avaz mishe va kolle column character mishe

g2$GEOid <- paste(g2$STATE_FIPS,g2$CZ_FIPS,sep = "")
g2$GEOid <- as.numeric(g2$GEOid)

df_c <- transform(g2, County_Index=match(GEOid, unique(GEOid)))
Tornado<- filter(df_c,EVENT_TYPE %in% c("Tornado"))



Tornado$TOR_F_SCALE[which(Tornado$TOR_F_SCALE=="EF0")]<-"F0"
Tornado$TOR_F_SCALE[which(Tornado$TOR_F_SCALE=="EF1")]<-"F1"
Tornado$TOR_F_SCALE[which(Tornado$TOR_F_SCALE=="EF2")]<-"F2"
Tornado$TOR_F_SCALE[which(Tornado$TOR_F_SCALE=="EF3")]<-"F3"
Tornado$TOR_F_SCALE[which(Tornado$TOR_F_SCALE=="EF4")]<-"F4"
Tornado$TOR_F_SCALE[which(Tornado$TOR_F_SCALE=="EF5")]<-"F5"


F01<- filter(Tornado,TOR_F_SCALE %in% c("F0","F1"))
F01<-dcast(F01, YEAR ~ GEOid,fun.aggregate =length, value.var = "EVENT_TYPE")

F23<- filter(Tornado,TOR_F_SCALE %in% c("F2","F3"))
F23<-dcast(F23, YEAR ~ GEOid,fun.aggregate =length, value.var = "EVENT_TYPE")

F45<- filter(Tornado,TOR_F_SCALE %in% c("F4","F5"))
F45<-dcast(F45, YEAR ~ GEOid,fun.aggregate =length, value.var = "EVENT_TYPE")
F45 <- merge(F45,expand.grid(YEAR=1950:2018),  all=TRUE) 
F45[is.na(F45)]<-0

###### Aggregate or Group_by and Cross-table #####

count<-dcast(Tornado, YEAR ~ GEOid,fun.aggregate =length, value.var = "EVENT_TYPE")
count<-count[,1:3047]



####################### Property and Crop Loss ###########

Loss<-dcast(Tornado, YEAR ~ GEOid,fun.aggregate =sum, value.var = "DAMAGE_PROPERTY_num",na.rm=TRUE)
Loss<-Loss[,1:3047]
Loss[is.na(Loss)]<-0



CLoss<-dcast(Tornado, YEAR ~ GEOid,fun.aggregate =sum, value.var = "DAMAGE_CROPS_num",na.rm=TRUE)
CLoss<-CLoss[,1:3047]
CLoss[is.na(CLoss)]<-0

tLoss<-CLoss[,2:3047]+Loss[,2:3047]
tLoss$year<-CLoss[,1]
tLoss<-tLoss[,c(3047,1:3046)]





################ CPI adjustment #############

CPI<-read.csv("CPI.csv",header = TRUE)
CPI<-as.data.frame(CPI)
CPI<-as.matrix(CPI)


inflation_factor<-as.matrix(251.095/CPI[1:69,2])

tloss_adj<-matrix(NA,69,3047)
tloss_adj<-as.data.frame(tloss_adj)
colnames(tloss_adj)<-colnames(tLoss)

for (i in 2:3047){
  tloss_adj[1:69,i]<-tLoss[1:69,i]*inflation_factor}
tloss_adj[,1]<-tLoss[,1]



#   Log scale  
tloss_adj[tloss_adj==0]<-1
tloss_adj[,2:3047]<-log(tloss_adj[,2:3047])

#####################################
GEOID<-data.frame(GEOID=as.numeric(colnames(count)[2:3047]))

pop2000<-read.csv("co-est00int-tot.csv",header=TRUE ,stringsAsFactors = FALSE)

pop2000<-pop2000[-which(pop2000$COUNTY==0),]
pop2000$STATE=formatC(pop2000$STATE, width = 2, format = "d", flag = "0")

pop2000$COUNTY=formatC(pop2000$COUNTY, width = 3, format = "d", flag = "0") 

pop2000$GEOID <- paste(pop2000$STATE,pop2000$COUNTY,sep = "")
pop2000$GEOID <- as.numeric(pop2000$GEOID )
pop2000<-pop2000[,c(21,1:20)]

pop2010<-read.csv("PEP_2017_PEPANNRES_with_ann.csv",header=TRUE,stringsAsFactors = FALSE)





pop2000<-pop2000[pop2000$GEOID %in% GEOID$GEOID,]
pop2010<-pop2010[pop2010$GEOID %in% pop2000$GEOID,]



pop1990<-read.csv("population counties 1990-2000.csv",header=TRUE ,stringsAsFactors = FALSE)
pop1990<-pop1990[-which(nchar(pop1990$GEOID)<3),]

pop1990<-pop1990[pop1990$GEOID %in% pop2000$GEOID,]

pop<-cbind(pop1990,pop2000,pop2010)



###########################     Land area ( data of 2010 )    #########

Area<-read.csv("DEC_10_SF1_GCTPH1.US05PR_with_ann.csv",header=TRUE ,stringsAsFactors = FALSE)
LandArea<-Area[,c(5,12)]
LandArea<-LandArea[-1,]
colnames(LandArea)<-c("GEOID","SquareMile")
LandArea<-LandArea[-which(nchar(LandArea$GEOID)<3),]
LandArea$GEOID=formatC(LandArea$GEOID, width = 4, format = "d", flag = "0")

LandArea$GEOID<-as.numeric(LandArea$GEOID)
LandArea$SquareMile<-as.numeric(LandArea$SquareMile)


LandArea<-LandArea[LandArea$GEOID %in% pop2000$GEOID,]



popDensity<-pop

for (i in 1:28){popDensity[i,2:3010]<- pop[i,2:3010]/LandArea[,2]}


#####################   Urban Area (data of 2010) #################

Urban<-read.csv("PctUrbanRural_County.csv",header=TRUE ,stringsAsFactors = FALSE)
Urban<-Urban[,c(1,2,3,4,9)]
Urban$AREA_URBAN<-(Urban$AREA_URBAN*3.86102e-7)
Urban$STATE<- formatC(Urban$STATE, width = 2, format = "d", flag = "0")
Urban$COUNTY<- formatC(Urban$COUNTY, width = 3, format = "d", flag = "0")

Urban$GEOID <- paste(Urban$STATE,Urban$COUNTY,sep = "")
Urban$GEOID <- as.numeric(Urban$GEOID )
Urban<-Urban[Urban$GEOID %in% pop2000$GEOID,]



################## Elevation ######

elevation<-read.delim("POP_PLACES_20181001.txt", sep="|",header=TRUE,stringsAsFactors = FALSE)

elev<-elevation[,c(5,7,17)]

elev$STATE_NUMERIC<- formatC(elev$STATE_NUMERIC, width = 2, format = "d", flag = "0")
elev$COUNTY_NUMERIC<- formatC(elev$COUNTY_NUMERIC, width = 3, format = "d", flag = "0")

elev$GEOID <- paste(elev$STATE_NUMERIC,elev$COUNTY_NUMERIC,sep = "")
elev$GEOID <- as.numeric(elev$GEOID )

El<-as.data.frame(ddply(elev,~GEOID,summarise,Eleveation_ft=mean(ELEV_IN_FT,na.rm = TRUE)))



El<-El[El$GEOID %in% pop2000$GEOID,]


######## AG LAnd- No of operations based on aggi#######

AGland<-as.data.frame(read.csv("AG land.csv",header=TRUE ,stringsAsFactors = FALSE))
AGland1997<-filter(AGland,Year %in% c(1997))
AGland2002<-filter(AGland,Year %in% c(2002))
AGland2012<-filter(AGland,Year %in% c(2012))



AG<-AGland[,c(7,11,20)]
AG$State.ANSI<- formatC(AG$State.ANSI, width = 2, format = "d", flag = "0")
AG$County.ANSI<- formatC(AG$County.ANSI, width = 3, format = "d", flag = "0")
AG$GEOID <- paste(AG$State.ANSI,AG$County.ANSI,sep = "")
AG$GEOID <- as.numeric(AG$GEOID )

AGland_ave<-as.data.frame(ddply(AG,~GEOID,summarise,No_Operation=mean(Value,na.rm = TRUE)))

AGland_ave<-AGland_ave[AGland_ave$GEOID %in% pop2000$GEOID,]


####### Farm operation : total No of farm operations #########


farmland<-as.data.frame(read.csv("total No of farm operations.csv",header=TRUE ,stringsAsFactors = FALSE))



farm<-farmland[,c(7,11,20)]
farm$State.ANSI<- formatC(farm$State.ANSI, width = 2, format = "d", flag = "0")
farm$County.ANSI<- formatC(farm$County.ANSI, width = 3, format = "d", flag = "0")
farm$GEOID <- paste(farm$State.ANSI,farm$County.ANSI,sep = "")
farm$GEOID <- as.numeric(farm$GEOID )

farm_ave<-as.data.frame(ddply(farm,~GEOID,summarise,No_Operation=mean(Value,na.rm = TRUE)))

farm_ave<-farm_ave[farm_ave$GEOID %in% pop2000$GEOID,]



############# some counties even did not show up in Farm and AGland data, So, That means

## they have zero No. of operation. Now add those GEOID from pop2000 GEOIDs with a zero value ###############

zero_AG <-as.data.frame(pop2000$GEOID[!(pop2000$GEOID %in% AGland_ave$GEOID)])
zero_farms <- as.data.frame(pop2000$GEOID[!(pop2000$GEOID %in% farm_ave$GEOID)])

AGland_ave <-  merge(AGland_ave, expand.grid( GEOID=zero_AG[,1]),  all=TRUE) 
farm_ave <- merge(farm_ave, expand.grid( GEOID=zero_farms[,1]),  all=TRUE) 
farm_ave[is.na(farm_ave)]<-0

##################### Housing units ###########################

HU2000<-read.csv("hu-est00int-tot.csv",header=TRUE ,stringsAsFactors = FALSE)

HU2000<-HU2000[-which(is.na(HU2000$COUNTY)),]
HU2000$STATE=formatC(HU2000$STATE, width = 2, format = "d", flag = "0")
HU2000$COUNTY=formatC(HU2000$COUNTY, width = 3, format = "d", flag = "0") 
HU2000$GEOID <- paste(HU2000$STATE,HU2000$COUNTY,sep = "")
HU2000$GEOID <- as.numeric(HU2000$GEOID )
HU2000<-HU2000[,c(18,1:17)]
HU2000<-HU2000[,c(1,3,7:16)]

HU2000<-HU2000[HU2000$GEOID %in% pop2000$GEOID,]



HU2010<-as.data.frame(read.csv("PEP_2017_PEPANNHU_with_ann.csv",header=TRUE ,stringsAsFactors = FALSE))
HU2010<-HU2010[,c(2,6:13)]
colnames(HU2010)[1]<-"GEOID"
HU2010<-HU2010[-1,]
HU2010[,2:9]<-as.numeric(as.matrix((HU2010[,2:9])))
HU2010$GEOID<-as.numeric(HU2010$GEOID)
HU2010<-HU2010[HU2010$GEOID %in% pop2000$GEOID,]

HU<-cbind(HU2000,HU2010)



###### match counties in F01, F23, F45 ##########


GEOID<-data.frame(GEOID=pop2000$GEOID)
F01<-F01[,names(F01) %in% pop2000$GEOID,]
F23<-F23[,names(F23) %in% pop2000$GEOID,]
F45<-F45[,names(F45) %in% pop2000$GEOID,]


new_F23 <-as.data.frame(matrix(NA,69,3009))

colnames(new_F23) <- pop2000$GEOID

keep <- as.matrix(colnames(F23))
new_F23[,(names(new_F23) %in% keep)]<- F23[,1:ncol(F23)]
new_F23[is.na(new_F23)]<- 0



new_F45 <-as.data.frame(matrix(NA,69,3009))

colnames(new_F45) <- pop2000$GEOID

keep <- as.matrix(colnames(F45))
new_F45[,(names(new_F45) %in% keep)]<- F45[,1:ncol(F45)]
new_F45[is.na(new_F45)]<- 0

new_F01 <-as.data.frame(matrix(NA,69,3009))

colnames(new_F01) <- pop2000$GEOID

keep <- as.matrix(colnames(F01))
new_F01[,(names(new_F01) %in% keep)]<- F01[,1:ncol(F01)]
new_F01[is.na(new_F01)]<- 0



############################ Latitude and  Longitude ###########


counties<-read.csv("UScounties_TableToExcel.csv",header=TRUE ,stringsAsFactors = FALSE)


lat<-data.frame(GEOID=as.numeric(colnames(count)),lat=NA)
long<-data.frame(GEOID=as.numeric(colnames(count)),lon=NA)
Region <-data.frame(GEOID=as.numeric(colnames(count)),Region=NA)

for (i in 1:3009)
{
  lat[i,2]<-counties$Latitude[which(counties$FIPS==colnames(count)[i])]
  print(i)
}


for (i in 1:3009)
{
  long[i,2]<-counties$Longitude[which(counties$FIPS==colnames(count)[i])]
  print(i)
}


for (i in 1:3009)
{
  Region[i,2]<-Tornado$Region[which(Tornado$GEOid==colnames(count)[i])[1]]
  print(i)
}




#----------------------------------- Model------------------------------

Y<-tloss_adj[11:69,]

moderate<-F01[11:69,]
severe<-new_F23[11:69,2:658]
devast<-new_F45[11:69,2:658]

nsites = ncol(Y)
nyears = nrow(Y)


FarmUnit<-log(farm_ave[,2])
UrbanArea<-log((Urban[,5])+1)
popD<-log((popDensity[,2:658])+1)
popD<-as.data.frame(colMeans(popD,na.rm=TRUE))
popD<-popD[,1]
HousingUnit<-log(HU[,2:658])
HousingUnit<-as.data.frame(colMeans(HousingUnit,na.rm=TRUE))
HousingUnit<-HousingUnit[,1]

jags.data <- list('moderate' = moderate,'severe' = severe, 'devast' = devast,'Y' = Y,'popD'=popD,'HousingUnit'=HousingUnit,'UrbanArea'=UrbanArea, 'FarmUnit'= FarmUnit,'nsites' = nsites, 'nyears'=nyears)
jags.params <- c('Y.pred','alpha','beta1','beta2','beta3','a0','a1','a2','a3',
                 'b01','b02','b03','b04',
                 'b11','b21','b31','b41',
                 'b12','b22','b32','b42',
                 'b13','b23','b33','b43',
                 'sigma_alpha','sigma_beta1','sigma_beta2','sigma_beta3','sigma_Y') 
jags.inits <- NULL

cat(
  
  "model
  
  {          
  
  for (i in 1:nsites)
  {
  for (t in 1:nyears)
  {
  Y[t,i] ~ dnorm(mu_Y[t,i],sigma_Y[i])
  Y.pred[t,i] ~ dnorm(mu_Y[t,i],sigma_Y[i])
  
  
  mu_Y[t,i] <- alpha[i] + (beta1[i]*moderate[t,i]) + (beta2[i]*severe[t,i]) + (beta3[i]*devast[t,i])
  
  }
  
  
  
  
  
  }
  
  for (j in 1:nsites)
  {
  
  alpha[j]~ dnorm(mu_alpha[j],sigma_alpha[j])
  mu_alpha[j]<-a0+b01*popD[j]+b02*HousingUnit[j]+b03*UrbanArea[j]+b04*FarmUnit[j]
  
  
  beta1[j] ~ dnorm(mu_beta1[j], sigma_beta1[j])
  mu_beta1[j] = a1+b11*popD[j]+b21*HousingUnit[j]+b31*UrbanArea[j]+b41*FarmUnit[j] 
  
  beta2[j] ~ dnorm(mu_beta2[j], sigma_beta2[j])
  mu_beta2[j] = a2+b12*popD[j]+b22*HousingUnit[j]+b32*UrbanArea[j]+b42*FarmUnit[j] 
  
  beta3[j] ~ dnorm(mu_beta3[j], sigma_beta3[j])
  mu_beta3[j] = a3+b13*popD[j]+b23*HousingUnit[j]+b33*UrbanArea[j]+b43*FarmUnit[j] 
  
  sigma_Y[j]~ dunif(0,10)
  sigma_alpha[j] ~ dunif(0,10)
  
  sigma_beta1[j]~ dunif(0,10)
  sigma_beta2[j] ~ dunif(0,10)
  sigma_beta3[j] ~ dunif(0,10)
  
  
  }
  
  
  a0 ~ dnorm(0,0.01)
  a1 ~ dnorm(0,0.01)
  a2 ~ dnorm(0,0.01)
  a3 ~ dnorm(0,0.01)
  
  b01~ dnorm(0,0.01)
  b02~ dnorm(0,0.01)
  b03~ dnorm(0,0.01)
  b04~ dnorm(0,0.01)
  b11~ dnorm(0,0.01)
  b21~ dnorm(0,0.01)
  b31~ dnorm(0,0.01)
  b41~ dnorm(0,0.01)
  b12~ dnorm(0,0.01)
  b22~ dnorm(0,0.01)
  b32~ dnorm(0,0.01)
  b42~ dnorm(0,0.01)
  b13~ dnorm(0,0.01)
  b23~ dnorm(0,0.01)
  b33~ dnorm(0,0.01)
  b43~ dnorm(0,0.01)
  
  }",file="Tor_new2.bug" )


set.seed(123)
jagsfit<-jags(data=jags.data, inits=jags.inits, jags.params, n.iter=1000,n.chains=3, model.file='Tor_new2.bug')


#results
attach.jags(jagsfit)