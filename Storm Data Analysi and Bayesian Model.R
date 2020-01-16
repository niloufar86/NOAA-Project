
library(data.table)
library (tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(pivottabler)
library(stringr)
library(rgdal)
library(zoo)



library(maps)
library(maptools)
library(ggplot2)
library(corrplot)

require(colorspace)
require(colorRamps)
require(RColorBrewer)
require(grDevices)
library(fields)


library(rpca)
library(Matrix)



library(biwavelet)
library(modifiedmk)

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

#-------------------------- reading data 


setwd("C:/Users/niloufar/Desktop/Data folder")

pop_timeseries<-read.csv("pop_timeseries.csv",header=TRUE,stringsAsFactors = FALSE)
GDP<-read.csv("gdp.csv",header=TRUE,stringsAsFactors = FALSE)
GDP<-GDP[1:51,]


tornado<-read.csv("1950-2017_all_tornadoes.csv",header=TRUE,stringsAsFactors = FALSE)
tornado$event<-"Tornado"
length(which(tornado$f1==0))

length(which(tornado$loss==0))


tornado$FIPS=formatC(tornado$FIPS, width = 2, format = "d", flag = "0")

tornado$f1=formatC(tornado$f1, width = 3, format = "d", flag = "0") 
tornado$GEOid <- paste(tornado$FIPS,tornado$f1,sep = "")
tornado$GEOid <- as.numeric(tornado$GEOid)


#------------------ loss categories

n=min(which(tornado$loss<1 & tornado$loss>0))

tornado$loss[n:64335]<-tornado$loss[n:64335]*(1e6)
tornado$loss[which(tornado$loss==0)]<-NA

min(tornado$loss[n:64335],na.rm=TRUE)
class(tornado$loss)

#------------------Consistent loss categories-After 1996

tornado$loss[which(tornado$loss<=50 & tornado$loss>9)]<-"C1"
tornado$loss[which(tornado$loss<=500 & tornado$loss>50)]<-"C2"
tornado$loss[which(tornado$loss<=5000 & tornado$loss>500)]<-"C3"
tornado$loss[which(tornado$loss<=50000 & tornado$loss>5000)]<-"C4"
tornado$loss[which(tornado$loss<=500000 & tornado$loss>50000)]<-"C5"
tornado$loss[which(tornado$loss<=5000000 & tornado$loss>500000)]<-"C6"
tornado$loss[which(tornado$loss<=50000000 & tornado$loss>5000000)]<-"C7"
tornado$loss[which(tornado$loss<=500000000 & tornado$loss>50000000)]<-"C8"
tornado$loss[which(tornado$loss<=5000000000 & tornado$loss>500000000)]<-"C9"
tornado$loss[which(tornado$loss>5000000000)]<-"C10"


#------------------ Climate Regions



#Function

State_names<-state.abb

codes<-data.frame(state=state.abb,climate_code=NA)
 



finder <- function(i) {
  row_index <- which(tornado$state==State_names[i])

return(row_index)
}


#--------------------------------Counts and GEOID --------------------------



count<-dcast(tornado, year ~ FIPS,fun.aggregate =length, value.var = "event")
count<-count[,-c(1,53)]


#Filtering thoes states with more than 60% zero data

count2<-count[,-(which(colMeans(count==0) > 0.6))]


GEOID<-data.frame(GEOID=colnames(count2))
GEOID$GEOID<-as.character(GEOID$GEOID)
for (i in 1:47)
{GEOID$STATE[i]<-tornado$state[which(tornado$FIPS==GEOID$GEOID[i])[1]]
}

Region <-data.frame(GEOID=as.numeric(colnames(count2)),Region=NA)

for (i in 1:47)
{
  Region[i,2]<-tornado$Region[which(tornado$FIPS==colnames(count2)[i])[1]]
  print(i)
}

GDP2<-GDP[-(which(colMeans(count==0) > 0.6)),]



#---------------- -------------Boxplots and outliers ----------------


par(mar=c(5, 5, 4, 4))

colnames(count2)<-GEOID$STATE
boxplot(count2,col="gray80",pch=16,cex=1.6,cex.axis=1.6,cex.lab=1.7,las=2,staplewex = 1.2, outwex =1,ylab="Number of Tornado")




#------------------------------ Robust PCA --------------------------
options(warn=1)

Count<-as.matrix(count2[1:69,])
xcent = sweep(Count,2,colMeans(Count))

rob <-rpca(xcent)
summary(rob)



rankMatrix(rob$L)
rankMatrix(Count)

#---------------- Scores and Correlation loadings

score<-as.data.frame(rob$L.svd$u)
score$year<-c(1950:2018)
score<-score[,c(28,1:27)]


score1<-as.matrix(rob$L.svd$u[,1])
score2<-as.matrix(rob$L.svd$u[,2])
score3<-as.matrix(rob$L.svd$u[,3])




corload1<-cor(Count,score1)
corload2<-cor(Count,score2)
corload3<-cor(Count,score3)

#Finding Significant Loadings for each PC:

Pcorload<-as.data.frame(array(NA,c(47,3)))
for (i in 1:47){
  Pcorload[i,1]<-cor.test(Count[,i],score1)$p.value
  Pcorload[i,2]<-cor.test(Count[,i],score2)$p.value
  Pcorload[i,3]<-cor.test(Count[,i],score3)$p.value
}



dn1<-data.frame(Geoid=as.numeric(GEOID$GEOID)[which(Pcorload[,1]<0.05)],value=corload1[which(Pcorload[,1]<0.05),1])
dn2<-data.frame(Geoid=as.numeric(GEOID$GEOID)[which(Pcorload[,2]<0.05)],value=corload2[which(Pcorload[,2]<0.05),1])
dn3<-data.frame(Geoid=as.numeric(GEOID$GEOID)[which(Pcorload[,3]<0.05)],value=corload3[which(Pcorload[,3]<0.05),1])


#---------------- image plot of sparse matrix


rotate <- function(x) t(apply(x, 2, rev))

sparse<-as.matrix(rob$S)
sparse[which(sparse==0)]<-NA

max(sparse,na.rm = TRUE)
min(sparse,na.rm = TRUE)





par(mar=c(5,4.5,4,7))
image(rotate(sparse),col = blue2red(7),breaks=c(-50,-25,0,25,50,75,100,125),axes=F)
mtext(text=c(2018:1950), side=2, line=0.3, at=seq(0,1,(1/68)), las=1, cex=0.7)

image.plot(rotate(sparse),legend.only=T,col = blue2red(7),breaks=c(-50,-25,0,25,50,75,100,125))
title(main = "sparse component image", font.main = 4)



#------------- outliers in sparse matrix


sparse<-rob$S


myFunc <- function(x) sd(x) 

sd_sparse<-as.data.frame(sapply(sparse, myFunc))


outlier<-as.data.frame(array(NA, c(69,47)))
for ( i in 1:47)
  
{
  outlier[which(sparse[1:69,i]>(2*sd_sparse[i,1]) | sparse[1:69,i]<(-(2*sd_sparse[i,1]))),i]<-1
}

colnames(sparse)=colnames(count2)
boxplot(sparse,col="gray80",pch=16)


outlier$total<-rowSums(outlier,na.rm=TRUE)
plot(c(1950:2018),outlier$total,pch=16,xlab="Year",ylab='Number of outliers',cex=1.6,cex.lab=1.6,cex.axis=1.6)




#------------- PCA variance barplot 


percentvar =round(as.matrix((rob$L.svd$d)^2/sum((rob$L.svd$d)^2)),digits=2)

b=barplot(percentvar[1:6,1], ylab="Explained variance", ylim=c(0,1),
          names=c("PC1","PC2","PC3","PC4","PC5","PC6"),space=0.2,cex.names = 1.7, cex.lab=1.7,cex.axis=1.8)
text(x = b, y = percentvar[1:6,1], label = percentvar[1:6,1], pos = 3, cex = 1.6, col="tomato3")



#------------- PCA corload maps



data("state.fips")
st <- map_data("state")
colnames(st)[5]<-"polyname"

state.fips<-data.frame(state.fips)
state.fips$polyname[which(state.fips$abb=="NY")]="new york"
state.fips$polyname[which(state.fips$abb=="NC")]="north carolina"
state.fips$polyname[which(state.fips$abb=="WA")]="washington"
state.fips$polyname[which(state.fips$abb=="VA")]="virginia"
state.fips$polyname[which(state.fips$abb=="MA")]="massachusetts"
state.fips$polyname[which(state.fips$abb=="MI")]="michigan"


st2 <- st %>%
  left_join(state.fips, by="polyname")

df<-data.frame(Geoid=as.numeric(GEOID$GEOID),value=corload3[,1])

colnames(df)[2]<-"PC3"

st2.df <- full_join(st2, df, by=c("fips" = "Geoid") )
st2.dn <- inner_join(st2, dn3, by=c("fips" = "Geoid") )

ggplot(st2.df, aes(long, lat,group = group)) + 
  geom_polygon(aes(fill = PC3), color="black")  +
  geom_path(data = st2.dn[1:nrow(st2.dn),], aes(long, lat, group = group), 
            color = "black", size = 1.5) +
  ggtitle("")+ 
  
  scale_fill_gradient2(high = "#0056B7",low = "red1",mid="white",limits = c(-0.45, 0.45),name="",midpoint =0,space ="Lab")+
  
  theme(legend.text=element_text(size=20),legend.key.size = unit(0.8, "cm"))+
  
  coord_quickmap()




#-------------- Scores plots

par(mar=c(5, 5, 4, 4))

t=c(1950:2018)

plot(t,rob$L.svd$u[,1],pch=16,xlab='Year',ylab='Scores-PC1',col='black',cex=1.8,cex.lab=2,cex.axis=2)
lines(lowess(t,rob$L.svd$u[,1],f=1/9),lwd=4,col="hotpink3")

plot(t,rob$L.svd$u[,2],pch=16,xlab='Year',ylab='Scores-PC2',col='black',cex=1.8,cex.lab=2,cex.axis=2)
lines(lowess(t,rob$L.svd$u[,2],f=1/9),lwd=4,col="hotpink3")


#-------------------------------Man kendall test Trend anlysis on first PC---------------------

#Autocorrelation Check

acf(score1)
pacf(score1)

# we saw autocorrelation which means we can not use Man kendall test as it is;


#Modified Man Kendall test
mmkh(score1[,1],ci=0.95)

#---------------------------------- Wavelet--------------------------



par(mar=c(5,5,3,8))

score<-as.data.frame(rob$L.svd$u[,4])
score$year<-c(1950:2018)
score<-score[,2:1]
colnames(score)[2]<-"value"
options(digits=3)
wlt=wt(score,sig.level = 0.95)

plot(wlt, type="power.corr.norm", xlab="Year",plot.cb = TRUE,cex.axis=1.9,cex.lab=1.9);

options(digits=7)

#---------------------------------Cliamte indexes -------------------------------

setwd("C:/Users/niloufar/Desktop/climate indices")


options(digits=7)
climate=read.csv("Climate_Indices_1943_2017.csv",header=TRUE,stringsAsFactors = FALSE)



#-----------------ENSO-  Nino 3.4 

Nino34<- read.delim("inino5_2.dat",sep="",skip=8, header = FALSE)
colnames(Nino34)<-c("year","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
Nino34<- Nino34[95:163,]

min(Nino34[,2:13],na.rm=TRUE)
Nino34[Nino34==-999.9]<-NA


Nino34$annual<-rowMeans(Nino34[,2:13],na.rm=TRUE)


#------------------------------ AO   --------------------#

AO <- read.delim("AO.dat",sep="",skip=1, header = FALSE)
AO<-AO[1:69,]

min(AO[,2:13],na.rm=TRUE)


AO$annual<-rowMeans(AO[,2:13],na.rm=TRUE)


#----------------------------  NAO -----------------------------
NAO <- read.delim("nao_3dp.dat",sep="",skip=1, header = FALSE)
NAO<-NAO[129:197,1:13]

min(NAO[,2:13],na.rm=TRUE)

NAO[NAO==-99.99]<-NA

NAO$annual<-rowMeans(NAO[,2:13],na.rm=TRUE)

#---------------------------- SOI---------------------

SOI <- read.delim("icpc_soi.dat",sep="",skip=8, header = FALSE)
SOI<-SOI[65:133,]

min(SOI[,2:13])

SOI$annual<-rowMeans(SOI[,2:13],na.rm=TRUE)

#-------------------------------    PDO    ------------------------#

Pdata <- read.csv("PDO.csv",header=FALSE,stringsAsFactors = FALSE)

Pdata<-Pdata[1:828,]


PDO<-as.data.frame(array(NA,c(69,13)))
PDO[,1]<-c(1950:2018)

for ( i in 0:68)
  
{ PDO [i+1,2:13]<- Pdata[(12*i+1):(12*i+12),2]
}

min(PDO[,2:13])
PDO$annual<-rowMeans(PDO[,2:13],na.rm=TRUE)


#---------------------------- AMO ---------------------------------

AMO<- read.delim("iamo_ersst_ts_2.dat",sep="",skip=47, header = FALSE)
colnames(AMO)<-c("year","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
AMO<-AMO[97:165,]

min(AMO[,2:13],na.rm=TRUE)
AMO[AMO==-999.9]<-NA

AMO$annual<-rowMeans(AMO[,2:13],na.rm=TRUE)




#######


ao<-AO[,14]
nino34<-Nino34[,14]
soi<-SOI[,14]
amo<-AMO[,14]
nao<-NAO[,14]
pdo<-PDO[,14]


#------------------------ scatter plots of climate index correlation with scores--------------

model1<-lm(score1~ao)
pval=round(summary(model1)$coefficients[2,4],digits=3)
plot(ao,score1,pch=16,xlab="AO index", ylab="PC1",cex=1.5,cex.lab=2,cex.axis=2,cex.main=1.6)
lines(lowess(ao,score1,f=1/2),col="blue",lwd=4)



#------------------------------- correlation with population --------------

A=as.data.frame(t(cor(score1,pop2)))



pvalue<-as.data.frame(array(NA,c(47,1)))
for (i in 1:47){
  pvalue[i,1]<-cor.test(score1,pop2[,i])$p.value
  
}
A$GEOID<-GEOID$GEOID

#write.csv(A,file="A.csv",row.names=FALSE)

data("state.fips")
st <- map_data("state")
colnames(st)[5]<-"polyname"

state.fips<-data.frame(state.fips)
state.fips$polyname[which(state.fips$abb=="NY")]="new york"
state.fips$polyname[which(state.fips$abb=="NC")]="north carolina"
state.fips$polyname[which(state.fips$abb=="WA")]="washington"
state.fips$polyname[which(state.fips$abb=="VA")]="virginia"
state.fips$polyname[which(state.fips$abb=="MA")]="massachusetts"
state.fips$polyname[which(state.fips$abb=="MI")]="michigan"


st2 <- st %>%
  left_join(state.fips, by="polyname")



df<-data.frame(Geoid=as.numeric(GEOID$GEOID),value=A[,1])

colnames(df)[2]<-"val"

st2.df <- full_join(st2, df, by=c("fips" = "Geoid") )





class <- classIntervals(A[,1], 3, style="quantile")

ggplot(st2.df, aes(long, lat,group = group)) + 
  geom_polygon(aes(fill = val), color="black")  +
  scale_fill_gradient2(high = "#5A6CC6",low = "red1",mid="white",name="",midpoint =0,space ="Lab")+
  
  
  
  theme(legend.text=element_text(size=18),legend.key.size = unit(1.2, "cm"),
        legend.key.width = unit(1, "cm") )+
  coord_quickmap()





#----------------------------------- Bayesian Model ------------------------

gc()
setwd("C:/Users/nilu/Desktop/Transfer to 117/climate indices")


climate=read.csv("Climate_Indices_1943_2017.csv",header=TRUE,stringsAsFactors = FALSE)
lastyear<-c(2018,0.2029166667,-0.19333333,0.266666667,0.540300000,-0.0874661354,-4.566667e-01)
climate<-rbind(climate,lastyear)

ao<-climate[,2]
nino34<-climate[,3]
soi<-climate[,4]
amo<-climate[,6]
nao<-climate[,5]
pdo<-climate[,7]
CO2<-log(CO2_an[,2])
cross<-nao*nino34






Area=read.csv("landarea.csv",header=TRUE,stringsAsFactors = FALSE)

pop<-array(NA,c(69,51))
for (i in 1:51)
{pop[,i]<-(pop_timeseries[1:69,i+1])/Area[i,2] }

pop<-log(pop)
#pop<-scale(pop2)
pop2<-pop[,-(which(colMeans(count==0) > 0.6))]

period<-array(0,c(69,1))
period[42:69,1]<-1
period<-period[,1]
DI<-period[1:69]

gdp<-log(GDP2[,2])
Reg<-as.factor(Region[,2])

nsites<-47
nyears<-69
################  indices ##################


Y<-count2[1:69,1:47]
rownames(Y)<-c(1:69)
no=length(which(Y==0))

nonzero<-as.data.frame(which(Y>0,arr.ind = TRUE))
zero<-as.data.frame(which(Y==0,arr.ind = TRUE))


###### Bayesian model ########

gc()
invisible(utils::memory.limit(32000))



set.seed(123)
jags.data <- list( 'Y' = Y,'nsites'=nsites,'nyears'=nyears,'cross'=cross,
                   'pop2'=pop2,'DI'=DI,
                   'nino34'=nino34,
                   'amo'=amo,'nao'=nao,'pdo'=pdo,
                   'Reg'=Reg,'gdp'=gdp)

jags.params <- c('Y.pred','lambda','alpha','beta',
                 'a0','a1','a2','a3','a4','a5','a6','a7',
                 
                 
                 'b10','b11','b12','b13','b14','b15','b16','b17',
                 'b20','b21','b22','b23','b24','b25','b26','b27',
                 
                 'tau_alpha','tau_beta1','tau_beta2','tau_beta3','tau_beta4',
                 'tau_beta5','tau_beta6','tau_beta7',
                 
                 'sigma_alpha','sigma_beta1','sigma_beta2','sigma_beta3','sigma_beta4',
                 'sigma_beta5','sigma_beta6','sigma_beta7') 
jags.inits <- NULL




cat(
  "
  model  
  {          
  ### -------------------------------- Level 1 -------------------------------------### 
  
  
  for (i in 1:nsites)
  {
  for (t in 1:nyears)
  {
  
  Y[t,i] ~ dpois(lambda[t,i])
  
  log(lambda[t,i])  <- alpha[i]+
  
  (beta[i,1]*pop2[t,i])+(beta[i,2]*DI[t])+
  
  (beta[i,3]*nino34[t])+(beta[i,4]*amo[t])+
  
  (beta[i,5]*nao[t])+(beta[i,6]*cross[t])+
  
  (beta[i,7]*pdo[t])
  
  Y.pred[t,i] ~ dpois(lambda[t,i])
  
  }
  }
  
  
  ##------------------ Hierarchical Level for Covariance on alpha and beta ------------------ ##
  
  for (j in 1:nsites)
  {
  
  alpha[j]~ dnorm(mu_alpha[j],tau_alpha)
  mu_alpha[j]<- a0 + b10*gdp[j] + b20*Reg[j]
  
  beta[j,1] ~ dnorm(mu_beta1[j], tau_beta1)
  mu_beta1[j] <- a1 + b11*gdp[j] + b21*Reg[j]
  
  beta[j,2] ~ dnorm(mu_beta2[j], tau_beta2)
  mu_beta2[j] <- a2 + b12*gdp[j] + b22*Reg[j] 
  
  beta[j,3] ~ dnorm(mu_beta3[j], tau_beta3)
  mu_beta3[j] <- a3 + b13*gdp[j] + b23*Reg[j] 
  
  beta[j,4] ~ dnorm(mu_beta4[j], tau_beta4)
  mu_beta4[j] <- a4 + b14*gdp[j] + b24*Reg[j]
  
  beta[j,5] ~ dnorm(mu_beta5[j], tau_beta5)
  mu_beta5[j] <- a5 + b15*gdp[j] + b25*Reg[j]
  
  beta[j,6] ~ dnorm(mu_beta6[j], tau_beta6)
  mu_beta6[j] <- a6 + b16*gdp[j] + b26*Reg[j] 
  
  beta[j,7] ~ dnorm(mu_beta7[j], tau_beta7)
  mu_beta7[j] <- a7 + b17*gdp[j] + b27*Reg[j] 
  
  
  
  
  }
  
  
  tau_alpha <- pow(sigma_alpha,-2)
  tau_beta1 <- pow(sigma_beta1,-2)
  tau_beta2 <- pow(sigma_beta2,-2)
  tau_beta3 <- pow(sigma_beta3,-2)
  tau_beta4 <- pow(sigma_beta4,-2)
  tau_beta5 <- pow(sigma_beta5,-2)
  tau_beta6 <- pow(sigma_beta6,-2)
  tau_beta7 <- pow(sigma_beta7,-2)
  
  
  
  sigma_alpha ~ dunif(10,100)
  sigma_beta1 ~ dunif(10,100)
  sigma_beta2 ~ dunif(10,100)
  sigma_beta3 ~ dunif(10,100)
  sigma_beta4 ~ dunif(10,100)
  sigma_beta5 ~ dunif(10,100)
  sigma_beta6 ~ dunif(10,100)
  sigma_beta7 ~ dunif(10,100)
  
  
  
  
  ############## Priors for level 1 #####################
  a0 ~ dnorm(0,0.001)
  a1 ~ dnorm(0,0.001)
  a2 ~ dnorm(0,0.001)
  a3 ~ dnorm(0,0.001)
  a4 ~ dnorm(0,0.001)
  a5 ~ dnorm(0,0.001)
  a6 ~ dnorm(0,0.001)
  a7 ~ dnorm(0,0.001)
  
  b10~ dnorm(0,0.001)
  b20~ dnorm(0,0.001)
  b11~ dnorm(0,0.001)
  b21~ dnorm(0,0.001)
  b12~ dnorm(0,0.001)
  b22~ dnorm(0,0.001)
  b13~ dnorm(0,0.001)
  b23~ dnorm(0,0.001)
  b14~ dnorm(0,0.001)
  b24~ dnorm(0,0.001)
  b15~ dnorm(0,0.001)
  b25~ dnorm(0,0.001)
  b16~ dnorm(0,0.001)
  b26~ dnorm(0,0.001)
  b17~ dnorm(0,0.001)
  b27~ dnorm(0,0.001)
  
  
  
  
  
  }",file="July18.bug" )

gc()
jagsfit<-jags(data=jags.data, inits=jags.inits, jags.params, n.iter=10000, n.chains=3, model.file='July18.bug')


#results
attach.jags(jagsfit)



#-------------------------- Output

R_squared<-as.data.frame(array(NA,c(nsites,2)))
adj_R<-as.data.frame(array(NA,c(nsites,2)))
R_squared[,1]<-colnames(Y)
adj_R[,1]<-colnames(Y)
colnames(adj_R)=c("fips","val")


zarib=(nyears-1)/(nyears-7-1)

Y_mean<-matrix(NA,ncol=nsites,nrow=1)
salam<-matrix(NA,ncol=nsites,nrow=nyears)
aleik<-matrix(NA,ncol=nsites,nrow=nyears)

Y_p<-jagsfit$BUGSoutput$median$Y.pred
#Y_p<-Ypred_median

for (j in 1:nsites)
{
  for (i in 1:nyears)
  {
    
    Y_mean[1,j]<-mean(Y[,j],na.rm=TRUE)
    salam[i,j]<-(Y[i,j]- Y_p[i,j])^2
    aleik[i,j]<-(Y[i,j]-Y_mean[1,j])^2
    
  }
}

for (j in 1:nsites)
{
  R_squared[j,2]<-1-((sum(salam[,j]))/(sum(aleik[,j])))
  #adj_R[j,2]<-1-((1-R_squared[j,2])*zarib)
  adj_R[j,2]<-1-(((sum(salam[,j]))/(sum(aleik[,j])))*zarib)
}



#------------------ Map of R-Squred--------------------

library(RColorBrewer)
library(classInt)


data("state.fips")
st <- map_data("state")
colnames(st)[5]<-"polyname"

state.fips<-data.frame(state.fips)
state.fips$polyname[which(state.fips$abb=="NY")]="new york"
state.fips$polyname[which(state.fips$abb=="NC")]="north carolina"
state.fips$polyname[which(state.fips$abb=="WA")]="washington"
state.fips$polyname[which(state.fips$abb=="VA")]="virginia"
state.fips$polyname[which(state.fips$abb=="MA")]="massachusetts"
state.fips$polyname[which(state.fips$abb=="MI")]="michigan"


st2 <- st %>%
  left_join(state.fips, by="polyname")



df<-data.frame(Geoid=as.numeric(GEOID$GEOID),value=adj_R[,2])

colnames(df)[2]<-"R_Squared"

st2.df <- full_join(st2, df, by=c("fips" = "Geoid") )



#cols=brewer.pal(8,"Blues")
#class <- classIntervals(adj_R[,2], nc, style="quantile",dataPrecision = 0)

ggplot(st2.df, aes(long, lat,group = group)) + 
  geom_polygon(aes(fill = R_Squared), color="black")  +
  
  
  #scale_fill_gradientn(colours =cols , na.value = "gray50",breaks=class$brks,labels=round(class$brks,2))+
  scale_fill_gradient2(high = "#0056B7",low = "red1",mid="white",name="",midpoint =0,limits = c(-0.65, 0.65),space ="Lab")+
  
  guides(fill = guide_colourbar(barheight = 10, direction = "vertical",
                                title.position="top", title.hjust = 1,title.vjust = 1))+  
  
  theme(legend.text=element_text(size=20),legend.key.size = unit(0.8, "cm"))+
  coord_quickmap()


#-----------------------------Posteriors -----------------

Beta_median<-jagsfit$BUGSoutput$median$beta

Betaa<-jagsfit$BUGSoutput$sims.list$beta

P<-array(NA,c(47,7))

for  (j in 1:7)
{for (i in 1:47)
{P[i,j]<-(length(which(Betaa[,i,j]>0)))/3000
}}





dn1<-data.frame(Geoid=as.numeric(GEOID[which(P[,1]>0.9 | P[,1]<0.1),1]),value=Beta_median[which(P[,1]>0.9 | P[,1]<0.1),1])


#---------------------- betas standard deviation



Beta_sd<-jagsfit$BUGSoutput$sd$beta

plot(colMeans(count2),Beta_sd[,1],pch=16,cex=2,cex.lab=1.7,cex.axis=1.8,ylab="Posterior Std,?? population ",
     xlab="Average tornado counts")




#--------- map of significant Betas boxplot

par(mar=c(4,5,3,5))
nam=c("Population","Radar","ENSO-Nino3.4","AMO","NAO","NAO*Nino3.4","PDO")
boxplot(Beta_median,col="gray60",names=nam,pch=16,cex=1.5,cex.axis=2,staplewex = 0.3, outwex =3,ylab="Beta",cex.lab=2)

alley<-c(40,34,14,25,38,39,13,11,23,3,13,1,22,8)

par(mfrow=c(1,8))

i=1 
boxplot(Beta_median[which(P[,i]>0.9 | P[,i]<0.1),i],ylim=c(-2.5,4), col="skyblue3",main=nam[i],pch=16,cex=2,cex.axis=1.8,staplewex = 0.3, outwex =2,ylab="Beta",cex.lab=1.8,cex.main=1.8)
boxplot(Beta_median[alley,i],ylim=c(-2.5,4),col="darksalmon",main=nam[i],pch=16,cex=2,cex.axis=1.8,staplewex = 0.3, outwex =2,ylab="",cex.lab=1.8, cex.main=1.8)



#------------------------------------



Beta_median[which(P[,1]<0.9 & P[,1]>0.1),1]<-0





#--------- map of Betas ------------

data("state.fips")
st <- map_data("state")
colnames(st)[5]<-"polyname"

state.fips<-data.frame(state.fips)
state.fips$polyname[which(state.fips$abb=="NY")]="new york"
state.fips$polyname[which(state.fips$abb=="NC")]="north carolina"
state.fips$polyname[which(state.fips$abb=="WA")]="washington"
state.fips$polyname[which(state.fips$abb=="VA")]="virginia"
state.fips$polyname[which(state.fips$abb=="MA")]="massachusetts"
state.fips$polyname[which(state.fips$abb=="MI")]="michigan"

st2 <- st %>%
  left_join(state.fips, by="polyname")



min(Beta_median[,1])
max(Beta_median[,1])



df<-data.frame(Geoid=as.numeric(GEOID$GEOID),
               value=as.numeric(Beta_median[,1]))

colnames(df)[2]<-"Beta"

st2.df <- full_join(st2, df, by=c("fips" = "Geoid") )
st2.dn <- inner_join(st2, dn1, by=c("fips" = "Geoid") )


ggplot(st2.df, aes(long, lat,group = group)) + 
  geom_polygon(aes(fill = Beta), color="black")  +
  geom_path(data = st2.dn[1:nrow(st2.dn),], aes(long, lat, group = group), 
            color = "black", size = 2) +
  ggtitle("")+
  
  scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-4.1, 4.1),space ="Lab")+
  #scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-1.7, 1.7),space ="Lab")+
  #scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-0.4, 0.7),space ="Lab")+
  #scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-3.3, 1.63),space ="Lab")+
  #scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-0.55, 0.8),space ="Lab")+
  #scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-1.1, 0.47),space ="Lab")+
  #scale_fill_gradient2(high = "#0056B7",name="",low = "red1",mid="white",midpoint =0,limits = c(-0.7, 0.34),space ="Lab")+
  
  
  
  
  guides(fill = guide_colourbar(barheight = 10, direction = "vertical",
                                title.position="top", title.hjust = 1,title.vjust = 1))+  
  
  theme(legend.text=element_text(size=20),legend.key.size = unit(0.8, "cm"))+
  coord_quickmap()





#-------------- map of tornado categories

tornado$magnitude[which(tornado$magnitude==(-9))]<-NA
F01<-filter(tornado,magnitude %in% c(0,1))
F2<-filter(tornado,magnitude %in% 2)
F3<-filter(tornado,magnitude %in% 3)
F45<-filter(tornado,magnitude %in% c(4,5))

rang<-as.matrix(c("green3","gold1","chocolate1","red3"))



map<-map("state", col="gray90", fill=TRUE,type = "l",add=FALSE)





points(x = F01$slon, y = F01$slat, pch=16,cex=0.6, lwd= 4, col = rang[1])
points(x = F2$slon, y = F2$slat, pch=16,cex=0.7, lwd= 4, col = rang[2])
points(x = F3$slon, y = F3$slat, pch=16,cex=0.9, lwd= 4, col = rang[3])
points(x = F45$slon, y = F45$slat, pch=16,cex=1, lwd= 4, col = rang[4])



legend("bottomright",as.character(c('F0-F1','F2','F3','F4-F5')), col = rang[1:4,1], pch = 16, cex=1.7 )









#---------------------- Second level figures --------------------

#--------------- Regions

par(mfrow = c(2,2),mar=c(11, 5, 3, 4))

Regi<-as.data.frame(as.numeric(as.character(Reg)))
names<-c('Northwest','Northern Plains', 'Upper Midwest','Northeast',
         'West','Southwest','South','Ohio Valley','Southeast')



plot(Regi[,1],Beta_median[,1],xaxt = 'n',cex=1.7,cex.lab=1.8,cex.axis=1.8,ylab="Beta-Population",
     xlab="",pch=16)
axis(1,at=c(1:9), labels=names, las=2, cex.axis=1.6)
abline(h=0, col="red",lwd=4, lty=2)








#------------------- GDP


par(mfrow = c(2,2),mar=c(5,5,2,4))

m<-array(NA,c(nrow(dn1),1))
for (i in 1:nrow(dn1))
{m[i]=which(as.numeric(GEOID$GEOID)==dn1$Geoid[i])
}


alley<-c(40,34,14,25,38,39,13,11,23,3,13,1,22,8)
GEOID$STATE[alley]


plot(gdp[m],Beta_median[m,1],pch=19,cex.lab=1.8,cex.axis=1.7,ylab="Beta-Population",
     xlab="log(GDP),million dollars", cex=2,ylim=c(-3,5))
points(gdp[-m],Beta_median[-m,1],pch="o",cex=1.6,col="gray40")


with(text(Beta_median[m,1]~gdp[m], labels = GEOID$STATE[m], pos = 4,cex=1.3,col="gray30"))
with(text(Beta_median[alley,1]~gdp[alley], labels = GEOID$STATE[alley], pos = 4,cex=1.3,col="red1"))

abline(h=0, col="royalblue2", lwd=3, lty=2)
lines(lowess(gdp,Beta_median[,1],f=0.5),col="gray50",lwd=2)





#----------Trend anlaysis on residuals


Resi<-array(NA,c(69,47))

for ( i in 1:69 )
{for (j in 1:47)
{
  Resi[i,j]<-(Y[i,j]- Y_p[i,j])
}}




coef<-array(NA, c(47,3))
colnames(coef)<-c("fips","alpha","time")

pvalue<-array(NA, c(47,3))
colnames(pvalue)<-c("fips","alpha","time")
pvalue[,1]<-colnames(Y)

t=c(1:69)


R<-array(NA, c(47,2))
R<-as.data.frame(R)
R[,1]<-colnames(Y)
colnames(R)<-c("fips","val")


for (i in 1:47){
  mod<- glm(formula=Resi[,i]~t,family = gaussian())
  
  coef[i,2:3]<-as.numeric(mod$coefficients)  
  
  
  R[i,2]<-1-((mod$deviance)/(mod$null.deviance))
  summ <- summary(mod)
  str(summ, max = 1)
  pvalue[i,2:3]<- as.numeric(summ$coefficients[,4])
  #Residuals[,i]<-mod$residuals
}

#-------- map of p-val

st2 <- st %>%
  left_join(state.fips, by="polyname")



df<-data.frame(Geoid=as.numeric(GEOID$GEOID),value=as.numeric(pvalue[,3]))

colnames(df)[2]<-"P_value"

st2.df <- full_join(st2, df, by=c("fips" = "Geoid") )




ggplot(st2.df, aes(long, lat,group = group)) + 
  geom_polygon(aes(fill = P_value), color="black")  +
  
  
  
  scale_fill_gradient2(name="",high = "#0056B7",low = "white",limits = c(0.2, 1),space ="Lab")+
  
  
  guides(fill = guide_colourbar(barheight = 10, direction = "vertical",
                                title.position="top", title.hjust = 1,title.vjust = 1))+  
  
  theme(legend.text=element_text(size=20),legend.key.size = unit(0.8, "cm"))+
  coord_quickmap()




#--------------- Doppler radar data


radar<-as.data.frame(read.csv("radar.csv",header=FALSE,stringsAsFactors = FALSE))

radar<-as.data.frame(radar[!is.na(radar),])
colnames(radar)<-"STATE"

radar_No<-as.data.frame(radar %>%
                          group_by(STATE) %>%
                          summarise(Count = length(STATE)) %>%
                          spread(STATE, Count))
#or

radar_No<-summarise(group_by(radar,STATE),Count = length(STATE))


#or
radar_No<-ddply(radar,~STATE,summarise,No = length(STATE))
radar_No<-ddply(radar,"STATE",summarise,No = length(STATE))


#--------------------------Adj R-squared vs radar and population change

radar2<-right_join(radar_No,GEOID,by="STATE")
radar2[is.na(radar2[,2]),2]<-0
plot(radar2$Count,adj_R$val,pch=16,cex=2,cex.lab=1.7,cex.axis=1.8,ylab="Adjusted R-squared",
     xlab="Number of WSR88-Doppler ")

pop_change<-array(NA,c(47,1))
for ( i in 1:47)
{pop_change[i,1]<-(pop2[69,i]-pop2[1,i])*100/pop2[1,i]}



par(mar=c(5,5,2,4))
plot(pop_change[which(pop_change<500),1],adj_R[which(pop_change<500),2],pch=16,cex=2,cex.lab=1.7,cex.axis=1.8,ylab="Adjusted R-squared",
     xlab="% Population change since 1950")
with(text(adj_R[alley,2]~pop_change[alley,1], labels = GEOID$STATE[alley], pos = 3,cex=1.3,col="red1"))


plot(pop_change[,1],adj_R[,2],pch=16,cex=2,cex.lab=1.7,cex.axis=1.8,ylab="Adjusted R-squared",
     xlab="% Population change since 1950")
with(text(adj_R[alley,2]~pop_change[alley,1], labels = GEOID$STATE[alley], pos = 4,cex=1.2,col="red1"))

#------------------------------------- Beta vs Beta------------------------

plot.new()

labs<- c("Beta-Population","Beta-Radar proxy","Beta-ENSO","Beta-AMO","Beta-NAO","Beta-NAO*ENSO","Beta-PDO")


par(mfrow = c(5,7),mar=c(5, 5, 3, 4))


for ( i in 1:7)
{  
  for (j in 1:7)
  {
  plot(Beta_median[,i],Beta_median[,j],pch=16,cex=1.4,cex.lab=1.4,cex.axis=1.4,ylab=labs[i],
       xlab=labs[j])
  
  }
}

#--------------------------- Beta sd vs Average count --------------------------------------


par(mfrow = c(3,3),mar=c(5, 5, 3, 4))

Beta_sd<-jagsfit$BUGSoutput$sd$beta



m<-array(NA,c(nrow(dn1),1))
for (i in 1:nrow(dn1))
{m[i]=which(as.numeric(GEOID$GEOID)==dn1$Geoid[i])
}

plot(colMeans(count2)[m],Beta_sd[m,1],pch=16,cex=2,cex.lab=1.7,cex.axis=1.8,ylab="Posterior Std,Beta population ",
     xlab="Average tornado counts")
points(colMeans(count2)[-m],Beta_sd[-m,1],pch="o",cex=1.6,col="deepskyblue2")





#------------------------------- Beta median vs Average counts

par(mfrow = c(3,3),mar=c(5, 5, 3, 4))



m<-array(NA,c(nrow(dn1),1))
for (i in 1:nrow(dn1))
{m[i]=which(as.numeric(GEOID$GEOID)==dn1$Geoid[i])
}

plot(colMeans(count2)[m],Beta_median[m,1],pch=16,cex=2,cex.lab=1.7,cex.axis=1.8,ylab="Beta-Population ",
     xlab="Average tornado counts")
points(colMeans(count2)[-m],Beta_sd[-m,1],pch="o",cex=1.6,col="deepskyblue2")


