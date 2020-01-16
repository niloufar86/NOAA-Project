gc()
library(data.table)
library (tidyr)
library(plyr)
library(stringr)
library(dplyr)
library(reshape2)
library(pivottabler)
library(zoo)
library(lubridate)
library(data.table)
library(rgdal)


library(maps)
library(maptools)
library(ggplot2)
library(ggmap)
library(RColorBrewer)
library(corrplot)

library(rpca)
library(Matrix)


library(geosphere)
library(matrixStats)

library(fpc)
library(cluster)
library(factoextra)
library(sp)

library(rjags)
library(R2jags)
library(jagstools)
library(runjags)
library(mvc)
library(SparseM)
library(MCMCpack)
library(parallel)
library(mcmcplots)
library(superdiag)
library(abind)
library(snow)
library(coda)
library(R2WinBUGS)






#--------------------------------data and  filtering F2+ tornadoes --------------------------

tornado<-read.csv("all_tornadoes.csv",header=TRUE,stringsAsFactors = FALSE)

tor<-filter(tornado,magnitude %in% c(2,3,4,5))

tor$date<-parse_date_time(tor$date, orders = c("ymd", "dmy", "mdy"))
tor$date<-as.Date(as.character(tor$date))
tor<-tor[-which(tor$slon==0),]

count<- dcast(tor, date ~ GEOid ,fun.aggregate =length, value.var = "event")



#------------------------------------------ Affected Area -County area


Area<-read.csv("DEC_10_SF1_GCTPH1.US05PR_with_ann.csv",header=TRUE ,stringsAsFactors = FALSE)
LandArea<-Area[,c(5,12)]
LandArea<-LandArea[-c(1,2,3),]
colnames(LandArea)<-c("GEOID","SquareMile")
LandArea<-LandArea[-which(nchar(LandArea$GEOID)<3),]
LandArea$GEOID=formatC(LandArea$GEOID, width = 4, format = "d", flag = "0")

LandArea$GEOID<-as.numeric(LandArea$GEOID)
LandArea$SquareMile<-as.numeric(LandArea$SquareMile)


GEOID<-data.frame(GEOID=as.numeric(as.character(colnames(count)[2:2529])))
Area<- left_join(GEOID,LandArea,"GEOID")

Area[is.na(Area[,2]),2]<-1
Affected<-as.data.frame(array(NA,c(4518,2528)))
for ( i in 1:4518)
{Affected[i,]<-count[i,2:2529]*Area[,2]
}


Affected_total<-rowSums(Affected)
hist(Affected_total,col="gray80",xlab="Area,square mile",main="")




#------------------- Finding events which hit  n number of  counties in same day


for ( i in 2:2529)
{count[which(count[,i]>0),i]<-1}


count_total<-as.data.frame(rowSums(count[,2:2529]))

quantile(count_total[,1],0.95)

count8<-count[-(which(count_total<9)),]


#---------------------------------- Histograms

plot(density(count_total[,1]),xlab="No of Counties",main="",cex.lab=1.4,cex.axis=1.4)
hist(count_total[,1],col="gray80",xlab="No of Counties",main="")

df<-data.frame(date=tor$date,lat=tor$slat,long=tor$slon, GEOID=tor$GEOid)
df$event<-"tornado"
matrix2<-dcast(df,date~lat,fun.aggregate =length, value.var = "event")
total<-as.data.frame(rowSums(matrix2[,2:2637]))

plot(density(total[,1]),xlab="No of locations",main="",cex.lab=1.4,cex.axis=1.4)
hist(total[,1],col="gray80",xlab="No of locations",main="")


#-----Seasonal histogram

Events<-data.frame(Date=count$date, t=c(1:4518) )
Events$year<-year(Events$Date)
Events$month<-month(Events$Date)
Events$total_counties<-count_total[,1]


MAM<-Events[which(Events$month==3 | Events$month==4 | Events$month==5),]

JJA<-Events[which(Events$month==6 | Events$month==7 | Events$month==8),]


hist(MAM[,5],col="gray80",xlab="No of Counties",main="March-April-May")

hist(JJA[,5],col="gray80",xlab="No of Counties",main="June-July-Aug")



#------- all tornado events lat and long which hit at least 8 county in one day


df<-data.frame(date=tor$date,lat=tor$slat,long=tor$slon, GEOID=tor$GEOid)


Allhits <- semi_join(df,count8,"date")


length (unique(df$GEOID))
length (unique(Allhits$GEOID))


Allhits$event<-"tornado"

matrix<-dcast(Allhits,date~lat,fun.aggregate =length, value.var = "event")

total<-rowSums(matrix[,2:1296])



newmat<-dcast(Allhits,date~GEOID,fun.aggregate =length, value.var = "event")

A<-newmat[,-(which(colMeans(newmat[,2:1359]==0) >0.99))]
B<-rowSums(A[,2:691])
length(which(B>7))


#---------------------- counties centroids -----------------

Centers<-data.frame(GEOID=unique(Allhits$GEOID),lat=NA,long=NA)


for ( i in 1:1358)
{
  Centers$lat[i]<-mean(Allhits$lat[which(Allhits$GEOID==Centers$GEOID[i])],na.rm=TRUE)
  Centers$long[i]<-mean(Allhits$long[which(Allhits$GEOID==Centers$GEOID[i])],na.rm=TRUE)
}




#------- Clustering Tornado Events

#---------------------------Finding  Optimum No of Cluster ----------------------


tot<-as.data.frame(total)
tot <- rbind(c(0), tot)
Optimum<-as.data.frame(array(NA,c(199,1)))

set.seed(123)


for (t in 1:199)
{
  
  w=tot[,1] 
  
  
  salam  <- data.frame(long=Allhits$long[(1+sum(w[1:t])):sum(w[1:(t+1)])], lat=Allhits$lat[(1+sum(w[1:t])):sum(w[1:(t+1)])])
  
  
  n<-total[t]
  
  pamk.best <- pamk(salam,krange = 1:(n-1),alpha=0.01)
  
  cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
  
  
  Optimum[t,1]<-pamk.best$nc
  
}



length(which(Optimum[,1]>1))



#krange in pamk:	
#integer vector. Numbers of clusters which are to be compared by 
#the average silhouette width criterion. Note: average silhouette width and Calinski-Harabasz
#can't estimate number of clusters nc=1. If 1 is included, a Duda-Hart test is applied 
#and 1 is estimated if this is not significant




#-----------------------------------------------------------------------------------


Events<-data.frame(date=matrix$date, t=c(1:199) , No_Clusters=Optimum[,1])
Events$year<-year(Events$date)




############## Clustering ###############


#--------- functions

geo.dist = function(df) {
  require(geosphere)
  d <- function(i,z){         # z[1:2] contain long, lat
    dist <- rep(0,nrow(z))
    dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
    return(dist)
  }
  dm <- do.call(cbind,lapply(1:nrow(df),d,df))
  return(as.dist(dm))
}


clust.sd = function(i, dat, clusters) {
  ind = (clusters == i)
  colSds(as.matrix(dat[ind,]),na.rm=TRUE)
}


clust.var = function(i, dat, clusters) {
  ind = (clusters == i)
  colVars(as.matrix(dat[ind,]),na.rm=TRUE)
}


#---------------------- clustering using pamk -------------------------------

Coresrows<-sum(Optimum)
Cores_lat<-as.data.frame(array(NA, c(1,1)))
Cores_lon<-as.data.frame(array(NA, c(1,1)))
sd_lat <-as.data.frame(array(NA, c(1,1)))
sd_lon <-as.data.frame(array(NA, c(1,1)))


tot <- as.data.frame(total)
tot <- rbind(c(0), tot)

for (t in 1:199)
{
  
  w <- tot[,1]
  
  
  
  salam  <- data.frame(long=Allhits$long[(1+sum(w[1:t])):sum(w[1:(t+1)])], lat=Allhits$lat[(1+sum(w[1:t])):sum(w[1:(t+1)])])
  
  n<-total[t]
  pamk<- pamk(salam,krange = 1:(n-1),alpha=0.01)
  K<-pamk$nc
  clusters<-pamk$pamobject$clustering
  
  
  
  centers<- pamk$pamobject$medoids
  sdsss<- sapply(unique(clusters), clust.sd, salam, clusters)
  varsss<- sapply(unique(clusters), clust.var, salam, clusters)
  
  
  
  
  Cores_lon1<-data.frame(V1=centers[1:K,1])
  Cores_lon<-rbind(Cores_lon,Cores_lon1)
  
  Cores_lat1<-data.frame(V1=centers[1:K,2])
  Cores_lat<-as.data.frame(rbind(Cores_lat,Cores_lat1))
  
  
  sd_lon1<-data.frame(V1=sdsss[1,1:K])
  sd_lon <- as.data.frame(rbind(sd_lon,sd_lon1))
  
  sd_lat1<-data.frame(V1=sdsss[2,1:K])
  sd_lat <- as.data.frame(rbind(sd_lat,sd_lat1))
  
  
  
}
rownames(Cores_lat)<-c(1:276)
rownames(Cores_lon)<-c(1:276)


#---------------------------------------------------Mapping----------- ----------------


#-------- finding cluster Cores for events happened in March-April-May


Events<-data.frame(date=matrix$date, t=c(1:199) , No_Clusters=Optimum[,1])
Events$year<-year(Events$date)
Events$month<-month(Events$date)

MAM<-Events[which(Events$month==3 | Events$month==4 | Events$month==5),]






colors=brewer.pal(n = 7, name = "YlGnBu")
colors=brewer.pal(n = 7, name = "YlOrRd")

map<-map("state", col="gainsboro", fill=TRUE,type = "l",add=FALSE)



tot <- as.data.frame(total)
tot <- rbind(c(0), tot)



for ( t in m7$t)
{
  
  w <- tot[,1]
  
  K<-Optimum[t,1]
  
  salam  <- data.frame(long=Allhits$long[(1+sum(w[1:t])):sum(w[1:(t+1)])], lat=Allhits$lat[(1+sum(w[1:t])):sum(w[1:(t+1)])])
  
  n<-total[t]
  pamk<- pamk(salam,krange = 1:(n-1),alpha=0.01)
  K<-pamk$nc
  clusters<-pamk$pamobject$clustering
  
  
  
  centers<- pamk$pamobject$medoids
  
  Cores_lon1<-data.frame(V1=centers[1:K,1])
  Cores_lat1<-data.frame(V1=centers[1:K,2])
  
  
  points(x = Cores_lon1[,1], y = Cores_lat1[,1], pch=16, lwd= 4,cex = 1.5, col=colors[7])
  
}



legend("bottomleft",c("1950-1960","1961-1970","1971-1980","1981-1990","1991-2000","2001-2010","2011-2018")
       , col = colors[1:7],ncol=1, pch = 16,cex = 1.3)





#------------- ploting Standard deviation and cores----

#-----sd

sd_lat<-sd_lat[2:276,1]
sd_lon<-sd_lon[2:276,1]

sd<-data.frame(data=matrix$date,year=year(matrix$date),lat=NA,long=NA)

for (t in 1:199)
{
  
  w=rbind(0,Optimum)
  
  sd[t,3]<-mean(sd_lat[(1+sum(w[1:t,1])):sum(w[1:(t+1),1])],na.rm=TRUE)
  sd[t,4]<-mean(sd_lon[(1+sum(w[1:t,1])):sum(w[1:(t+1),1])],na.rm=TRUE)
  
}

sd$Area<-sd$lat*(abs(sd$long))

plot(sd$year,sd$Area,pch=16,cex=1.5,xlab="", ylab="Average Cluster Area")
plot(sd$Area,pch=16,cex=1.5,xlab="time", ylab="Average Cluster Area")


sd$event<-"tornado"
Area_ave<-as.data.frame(ddply(sd,year~event,summarise,mean=mean(Area)))

plot(Area_ave$year,Area_ave$mean,pch=16,cex=1.5,xlab="", ylab="Average Cluster Area",cex.lab=1.5,cex.axis=1.4)
lines(lowess(Area_ave$year,Area_ave$mean,f=1/8),lwd=4,col="hotpink3")



#------ cores lat and long

Cores_lat<-Cores_lat[2:276,1]
Cores_lon<-Cores_lon[2:276,1]

Cores<-data.frame(data=matrix$date,year=year(matrix$date),lat=NA,long=NA)

for (t in 1:199)
{
  
  w=rbind(0,Optimum)
  
  Cores[t,3]<-mean(Cores_lat[(1+sum(w[1:t,1])):sum(w[1:(t+1),1])],na.rm=TRUE)
  Cores[t,4]<-mean(Cores_lon[(1+sum(w[1:t,1])):sum(w[1:(t+1),1])],na.rm=TRUE)
  
}

plot(Cores$year,Cores$lat,pch=16,cex=1.5,xlab="", ylab="latitude of cluster centroid")
plot(Cores$lat,pch=16,cex=1.5,xlab="Event, day index", ylab="latitude of cluster centroid")

plot(Cores$year,Cores$long,pch=16,cex=1.5,xlab="", ylab="longitude of cluster centroid")



Cores$event<-"tornado"
lat_ave<-as.data.frame(ddply(Cores,year~event,summarise,mean=mean(lat)))
long_ave<-as.data.frame(ddply(Cores,year~event,summarise,mean=mean(long)))

plot(lat_ave$year,lat_ave$mean,pch=16,cex=1.5,xlab="", ylab="Cluster Centroid, latitude",cex.lab=1.5,cex.axis=1.4)
lines(lowess(lat_ave$year,lat_ave$mean,f=1/8),lwd=4,col="blue")
plot(long_ave$year,long_ave$mean,pch=16,cex=1.5,xlab="", ylab="Cluster Centroid, longitude",cex.lab=1.5,cex.axis=1.4)
lines(lowess(long_ave$year,long_ave$mean,f=1/8),lwd=4,col="blue")




#------------------ map  ---------------------------

Events<-data.frame(date=matrix$date, t=c(1:199) , No_Clusters=Optimum[,1])
Events$year<-year(Events$date)
Events$month<-month(Events$date)

hist(Events$month,main="",col="gray90",ylab="frequency",xlab="month")




#colors<-brewer.pal(n = 7, name = "YlGnBu")
#name = "YlOrRd"

colors<-colorRampPalette(c("white", "red"), space = "rgb")(14)

colors<-colorRampPalette(c("white", "navy"), space = "rgb")(69)



library(viridis)
colors<-viridis_pal(option = "B")(13)  #magma
colors<-viridis_pal(option = "D")(13)  #inferno

vacol<-as.numeric(c(1:13))
colors<-rgb(0,0,vacol,maxColorValue = 13) #blue
colors<-rgb(vacol,0,0,maxColorValue = 13)  #red


map<-map("state", col="gainsboro", fill=TRUE,type = "l",add=FALSE)



tot <- as.data.frame(total)
tot <- rbind(c(0), tot)


for (i in 1:13)
  #for (i in 1:64)
  
{ 
  
  m<-Events[which(Events$year>1950+5*i & Events$year<1955+(5*i+1)),]
  #m<-Events[which(Events$year==1954+i ),]
  
  for ( t in m$t)
    
  {
    
    w <- tot[,1]
    
    K<-Optimum[t,1]
    
    salam  <- data.frame(long=Allhits$long[(1+sum(w[1:t])):sum(w[1:(t+1)])], lat=Allhits$lat[(1+sum(w[1:t])):sum(w[1:(t+1)])])
    
    n<-total[t]
    pamk<- pamk(salam,krange = 1:(n-1),alpha=0.01)
    K<-pamk$nc
    clusters<-pamk$pamobject$clustering
    
    
    
    centers<- pamk$pamobject$medoids
    
    Cores_lon1<-data.frame(V1=centers[1:K,1])
    Cores_lat1<-data.frame(V1=centers[1:K,2])
    
    
    points(x = Cores_lon1[,1], y = Cores_lat1[,1], 
           pch=16, lwd= 4,cex = 1.5, col=colors[i+1])
    
  }
  
}

legend("bottomleft",as.character(c(1950:2018)),
       col = colors[1:69], pch = 16,cex = 1,ncol=4)



legend("bottomleft",c("1950-1955","1956-1960","1961-1965","1966-1970","1971-1975","1976-1980","1981-1985","1986-1990",
                      '1991-1995','1996-2000','2001-2005','2006-2010','2011-2015','2015-2018')
       , col = colors[1:14],ncol=1, pch = 16,cex = 1.4)


#---------------------------

dev.off()


Events$event<-"tornado"
annual<-dcast(Events,year~event,fun.aggregate =length, value.var = "event")



plot(annual[,1],annual[,2],type="b",pch=16,cex=1.1,cex.lab=1.4,cex.axis=1.4,xlab="year",ylab="No of events")

lines(lowess(annual[,1],annual[,2],f=1/8),lwd=2,col="blue")




#relationship with ENSO


enso<-as.data.frame(Nino34$annual)
enso$year<-Nino34$year
annual_new<-full_join(enso,annual,"year")
annual_new[is.na(annual_new[,3]),3]<-0
plot(annual_new[,1],annual_new[,3],pch=16,cex=1.3,cex.lab=1.4,cex.axis=1.4,xlab="ENSO-Nino3.4 index",ylab="No of events")
lines(lowess(annual_new[,1],annual_new[,3],f=1/7),lwd=2,col="red")




#-----------------------------------------------------------------------------------------------

###### plot pamk"

pom<-pamk(salam,2:15)
plot(salam, col=pom$pamobject$clustering, pch=19, xlab="Counties", ylab="Total Damage", main="pamk(data))")


###### Group_by and Cross-table 

Tor<-as.data.frame(ddply(Tornado,YEAR~GEOid,summarise,Count=length(EVENT_TYPE)))
Tor2<-spread(Tor,GEOid,Count)
Tor3<-spread(Tor,YEAR,Count)

#or
salam<-as.data.frame(Tornado %>%
                       group_by(GEOid, YEAR) %>%
                       summarise(Count = length(EVENT_TYPE)) %>%
                       spread(GEOid, Count))

#or
count<-dcast(Tornado, YEAR ~ GEOid,fun.aggregate =length, value.var = "EVENT_TYPE")
count<-count[,1:3047]


#----------------------------- Annual Aggregates 

annual<-dcast(tor,year~event,fun.aggregate =length, value.var = "event")
#or
annual<-as.data.frame((ddply(tor,year~event,summarise,Count=length(event))))



dt<-as.data.frame((ddply(tor,year~date,summarise,Count=length(date))))
dt$event<-"tornado"


dt22<-as.data.frame((ddply(dt,year~event,summarise,count=length(date))))

#or

dt33<-aggregate(date~year,data=hi ,FUN = NROW)



#---------------------------------------- Chapter2--------------------------------
#-------------------------------------
#---------------------------
#---------------
#-------
#---



Events<-data.frame(date=matrix$date, t=c(1:199) , No_Clusters=Optimum[,1])
Events$year<-year(Events$date)
Events$month<-month(Events$date)
Events$event<-"tornado"
Events$date<-as.character(Events$date)

D<-data.frame(date=seq(as.Date("1952/03/21"), by = "day", length.out = 24392))
D$date<-as.character(D$date)

time<-left_join(D,Events,"date")
#time$event[is.na(time$event)]<-0
time$event[time$event=="tornado"]<-1
time$event<-as.numeric(time$event)
time$index<-c(1:24392)

#------------------------- Stripchart

stripchart(index~event,data=time,vertical=FALSE,cex=0.5,xlab="Time index,day",ylab="Event occurrence",cex.lab=1.5,cex.axis=1.5)

#---time to next event

time$arrival<-NA

time$arrival[which(time$event==1)]<-time$index[which(time$event==1)]

arrival<-data.frame(arrival=time$arrival[is.finite(time$arrival)], ind=c(1:199))

arrival$month<-time$month[which(time$event==1)]
arrival$year<-time$year[which(time$event==1)]



#----------------- Events return time


arrival$return<-NA
for ( i in 2:199)
{arrival$return[i]<-arrival$arrival[i]-arrival$arrival[i-1]
}



plot(arrival$ind,arrival$return,pch=16,type="b",
     ylab="Return Period",xlab="Event ID",cex.lab=1.5,cex.axis=1.5,main="Time between two consecutive events")
lines(lowess(arrival$ind,arrival$return,f=1/10),lwd=3,col="red")

plot(density(arrival$return[2:199]),cex.lab=1.5,cex.axis=1.5,lwd=2,main="")
hist(arrival$return,col="gray70",ylab="Frequency",xlab="return time,days",cex.lab=1.5,cex.axis=1.5)
hist(arrival$return,col="gray70",freq=FALSE,xlab="return time,days",cex.lab=1.5,cex.axis=1.5)




library(locfit)
density<-locfit(~ arrival$return[2:199],alpha=c(0.8,0))
#or density<-locfit(~ lp(arrival$return[2:199],nn=0.5,h=0))
plot(density,get.data=TRUE,xlab="Arrival time,days",cex.lab=1.4)
points(density,pch=16,col="red",cex=1.2) #???????????


y=fitted(density)
x=arrival$return[2:199]
plot(arrival$return[2:199],fitted(density),pch=16)



f<-data.frame(x=x,y=y)
plot(f$x,f$y)


g<-data.frame(x=x,y=log(y))
#g<-data.frame(x=x,y=1/y)
plot(g$x,g$y)


mod1<-glm(y~x,f,family=gaussian(link='log'))
summary(mod1)


mod2<-glm(y~x,g,family=gaussian)     #same as mod2<-lm(y~x,g)
summary(mod2)



mod3<-glm(y~x,f,family=Gamma(link='log'))
summary(mod3)





#--------------------- Exponential Distribution function

model <- nls(y ~ b*exp( -b * x), data = f, start = list(b = 0.01),trace = T)
#OR
#model <- nls(y ~ I(b*exp(1)^( b * x)), data = f, start = list(b = 0.01),trace = T)
summary(model)$coefficients[,1]
plot(f$x,f$y, main = "Fitted exponential function", sub = "Blue: fit; black: observed")
points(f$x, predict(model),  col = "blue")




R1<-1-(mod1$deviance/mod1$null.deviance)
R2<-1-(mod2$deviance/mod2$null.deviance)
R3<-1-(mod3$deviance/mod3$null.deviance)

plot(mod1$fitted.values,pch=16,ylab="Fitted values",cex.axis=1.5,cex.lab=1.5,cex.main=1.4, main="GLM gaussian(link=log)")

plot(mod2$fitted.values,pch=16,ylab="Fitted values",cex.axis=1.5,cex.lab=1.5,cex.main=1.4, main="GLM gaussian")

plot(mod3$fitted.values,pch=16,ylab="Fitted values",cex.axis=1.5,cex.lab=1.5,cex.main=1.4, main="GLM Gamma(link=log)")






plot(mod$residuals)
acf(mod$residuals)


#------------------- Seasonal analysis on return time -----------


MAM<-arrival[which(arrival$month==3| arrival$month==4 | arrival$month==5 ),]
plot(density(MAM$return[2:133]),cex.lab=1.5,cex.axis=1.5,lwd=2,main="March-Apr-May")
hist(MAM$return,col="gray70",ylab="Frequency",xlab="return time,days",cex.lab=1.5,cex.axis=1.5, main= "March-Apr-May")

JJA<-arrival[which(arrival$month==6| arrival$month==7 | arrival$month==8 ),]
plot(density(JJA$return),cex.lab=1.5,cex.axis=1.5,lwd=2,main="Jun-Jul-Aug")

SON<-arrival[which(arrival$month==9| arrival$month==10 | arrival$month==11 ),]
plot(density(SON$return),cex.lab=1.5,cex.axis=1.5,lwd=2,main="Sep-Oct-Nov")

DJF<-arrival[which(arrival$month==12| arrival$month==1|arrival$month==2),]
plot(density(DJF$return),cex.lab=1.5,cex.axis=1.5,lwd=2,main="Dec-Jan-Feb")



#---------------------------------- aggregate annual for return time  and in each season---------------------

f1<-aggregate(return~year,data=arrival,FUN=mean,na.rm=TRUE)
f2<-aggregate(return~year,data=MAM,FUN=mean,na.rm=TRUE)

plot(f1$return,type="b",pch=16, main="Annual Average return time : All the seasons (199 events)",
     ylab="Inter arrival time",xlab="year",cex.lab=1.5,cex.axis=1.4)
lines(lowess(f1$return,f=0.5),lwd=2,col="blue")




plot(f2$return,type="b",pch=16, main="Annual Average return time : March-April-MAy",
     ylab="Inter arrival time",xlab="year",cex.lab=1.5,cex.axis=1.4)
lines(lowess(f2$return,f=0.5),lwd=2,col="red")
plot(density(f2$return),cex.lab=1.5,cex.axis=1.5,lwd=2,main="Annual Average return time density plot :March-Apr-May")



#-------------------------- Lambda changes over time -----------------

#if (is.na(a)==TRUE || is.na(b)==TRUE)


Lambda<-as.data.frame(array(NA, c(36,2)))
colnames(Lambda)<-c('Lambda','R')

for (i in 1:36)
{  
  
  a=which(arrival$year==1952+i)[1]
  
  
  
  if ( is.na(a)==TRUE)
  {Lambda [i,]<- NA}
  
  else
    
  {
    if (length(which(arrival$year==1952+i+30))==0 )
    {Lambda [i,]<- NA}
    
    else
      
    { b=which(arrival$year==1952+i+30)[length(which(arrival$year==1952+i+30))]
    
    

    
    
    den<-locfit(~ arrival$return[a:b],alpha=c(0.8,0))
    
    
    y=fitted(den)
    x=arrival$return[a:b]
    
    
    
    f<-data.frame(x=x,y=y)
    
    #mod<-glm(y~x,f,family=gaussian(link='log'))
    
    
    mod <- nls(y ~ b*exp( -b * x), data = f, start = list(b = 0.01),trace = T)
    
    
    Lambda[i,1]<-summary(mod)$coefficients[1]
    
    }
  }
  
  
}



Lambda$label<-c(1983:2018)

Lambda2<-Lambda[-which(is.na(Lambda$Lambda)==TRUE),]

plot(Lambda$Lambda,type="b",lwd=2,pch=16, ylab="Lambda",xlab="Index of 30-year period",cex.lab=1.4)
text(Lambda$Lambda, labels = Lambda$label, pos = 3,col="blue")





#------------------- Lambda gradual convergence to one Lambda ------------


Teta<-as.data.frame(array(NA, c(65,3)))
colnames(Teta)<-c('Beta','Intercept','R')

for (i in 1:65)
{  
  a=which(arrival$year==1952+i)[1]
  b=which(arrival$year==2018)
  
  if (is.na(a)==TRUE)
    
  {Teta [i,]<- NA}
  
  else
    
    
  {
    
    den<-density(arrival$return[a:b])
    
    
    x=den$x[which(den$x>1)]
    y=den$y[which(den$x>1)]
    
    
    f<-data.frame(x=x,y=y)
    
    mod<-glm(y~x,f,family=gaussian(link='log'))
    
    Teta[i,1]<-mod$coefficients[2]
    Teta[i,2]<-mod$coefficients[1]
    Teta[i,3]<-1-(mod$deviance/mod$null.deviance)
    
  }
  
}

Teta$label<-c(1953:2017)

Teta2<-Teta[-which(is.na(Teta$Beta)==TRUE),]

plot(Teta$label,Teta$Beta,type="b",lwd=2,pch=16, ylab="Lambda",cex.lab=1.4,xlab="start year")


#----------------- map of first and last 30 year period --------------


Events<-data.frame(date=matrix$date, t=c(1:199) , No_Clusters=Optimum[,1])
Events$year<-year(Events$date)
Events$month<-month(Events$date)

colors<-colorRampPalette(c("white", "navy"), space = "rgb")(31)




colors<-colorRampPalette(c("white", "red"), space = "rgb")(31)

map<-map("state", col="gainsboro", fill=TRUE,type = "l",add=FALSE)



tot <- as.data.frame(total)
tot <- rbind(c(0), tot)



for ( i in 1:31)
{  
  
  m<-Events[which(Events$year==1987+i ),]
  
  for ( t in m$t)
    
  {
    
    w <- tot[,1]
    
    
    
    salam  <- data.frame(long=Allhits$long[(1+sum(w[1:t])):sum(w[1:(t+1)])], lat=Allhits$lat[(1+sum(w[1:t])):sum(w[1:(t+1)])])
    
    
    
    
    points(x = salam$long, y = salam$lat, 
           pch=16, lwd= 4,cex = 1.5, col=colors[i+1])
    
  }
  
  
}

       , col = colors[1:31],ncol=2, pch = 16,cex = 1)


write.csv(newmat,"Event-Counties.CSV",row.names = FALSE)
write.csv(A,"selected counties.CSV",row.names = FALSE)
