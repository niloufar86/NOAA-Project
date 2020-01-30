library(data.table)
library (tidyr)
library(plyr)
library(stringr)
library(dplyr)
library(reshape2)
library(pivottabler)
library(zoo)
library(rpca)
library(Matrix)
library(maps)
library(RColorBrewer)
library(ggplot2)
library(ggmap)
library(corrplot)
library(rgdal)

#----- data cleaning

files <- list.files(pattern="*.csv")
all_data <- do.call(rbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE)))

all_data$County_State <- paste(all_data$CZ_NAME,all_data$STATE,sep = ", ")
data<-all_data

data$DP_num <- str_extract(data$DAMAGE_PROPERTY, "[0-9.]+")                  
data$DP_char<-paste(str_extract(data$DAMAGE_PROPERTY, "[aA-zZ]+"))

char<-c("H","h","k","K","M","m","B","b")
num<-c(100,100,1000,1000,1e+06,1e+06,1e+09,1e+09)

for ( i in 1:8)
  
{
  data$DP_char<-gsub(char[i],num[i],data$DP_char,fixed = TRUE)
}



data$DP_char<-gsub("NA",1,data$DP_char,fixed = TRUE)


data$DAMAGE_PROPERTY_num<-(as.numeric(data$DP_char))*(as.numeric(data$DP_num))

data$BEGIN_YEARMONTH=as.numeric(data$BEGIN_YEARMONTH)
data$BEGIN_DAY=as.numeric(data$BEGIN_DAY)


data$BEGIN_DAY=formatC(data$BEGIN_DAY, width = 2, format = "d", flag = "0")

data$Date <- paste(data$BEGIN_YEARMONTH,data$BEGIN_DAY,sep = "")
data$Date <- as.numeric(data$Date)
data <- data[order(data$Date),]

data$Date <- gsub('^(.{4})(.*)$', '\\1-\\2', data$Date)
data$Date <- gsub('^(.{7})(.*)$', '\\1-\\2', data$Date)
data$Date <- as.Date(data$Date,"%Y-%m-%d")


data$STATE_FIPS=formatC(data$STATE_FIPS, width = 2, format = "d", flag = "0")
tornado<- filter(data,EVENT_TYPE %in% c("Tornado"))

#--------------------------------Counts and GEOID --------------------------
tornado$event<-"T"
count<-dcast(tornado, YEAR ~ STATE_FIPS,fun.aggregate =length, value.var = "event")
count<-count[,-1]

count2<-count[,-(which(colMeans(count==0) > 0.62))]


GEOID<-data.frame(GEOID=colnames(count2))

for (i in 1:47)
{GEOID$STATE[i]<-tornado$STATE[which(tornado$STATE_FIPS==GEOID$GEOID[i])[1]]
}

#---------------- -------------Boxplots and outliers ----------------

par(mar=c(5, 5, 4, 4))

boxplot(count2,col="gray80",pch=16,cex=1.6,cex.axis=1.6,cex.lab=1.7,las=2,staplewex = 1.2, outwex =1,ylab="Number of Tornado")


#----------------------------- population data ------------------------------------

pop5060<- read.delim("st5060ts.dat",sep="",skip=27, header = FALSE)
pop50_60<-cbind(pop5060[1:51,],pop5060[64:114,2:7])
pop50_60<-pop50_60[,-(c(2,13))]
colnames(pop50_60)<-c("STATE",1950:1959)
pop50_60$FIPS<-GEOID$GEOID[2:52]
pop50_60<-pop50_60[,c(12,1:11)]

for (i in 3:12)
{pop50_60[,i]<-as.numeric(gsub(",","",pop50_60[,i]))}
pop50_60[,3:12]<-pop50_60[,3:12]*1000


pop6070<- read.delim("st6070ts.dat",sep="",skip=23, header = FALSE)
pop60_70<-cbind(pop6070[1:51,],pop6070[60:110,2:7])
pop60_70<-pop60_70[,-(c(2,13))]
colnames(pop60_70)<-c("STATE",1960:1969)
pop60_70$FIPS<-GEOID$GEOID[2:52]
pop60_70<-pop60_70[,c(12,1:11)]
for (i in 3:12)
{pop60_70[,i]<-as.numeric(gsub(",","",pop60_70[,i]))}
pop60_70[,3:12]<-pop60_70[,3:12]*1000


pop7080<- read.delim("st7080ts.dat",sep="",skip=14, header = FALSE)
pop70_80<-cbind(pop7080[1:51,],pop7080[54:104,3:7])
pop70_80<-pop70_80[,-13]
colnames(pop70_80)<-c("FIPS","STATE",1970:1979)
pop70_80$FIPS<-GEOID$GEOID[2:52]

for (i in 3:12)
{pop70_80[,i]<-as.numeric(as.character(pop70_80[,i]))}


pop80_84<- read.delim("st8090ts.dat",sep="",skip=11, header = FALSE)
pop80_84<-pop80_84[1:51,]
pop85_89<- read.delim("st8090ts.dat",sep="",skip=70, header = FALSE)
pop85_89<-pop85_89[,-7]
pop80_90<-cbind(pop80_84[,1:6],pop85_89[,2:6])
pop80_90$FIPS<-GEOID$GEOID[2:52]
pop80_90<-pop80_90[,c(12,1:11)]
colnames(pop80_90)<-c("FIPS","STATE",1980:1989)

for (i in 3:12)
{pop80_90[,i]<-as.numeric(as.character(pop80_90[,i]))}


pop9099<-read.csv("population counties 1990-2000.csv",header=TRUE ,stringsAsFactors = FALSE)
pop9099<-pop9099[which(nchar(pop9099$GEOID)<3),]
pop90_99<-pop9099[,-c(3,14)]
pop90_99<-pop90_99[,c(1,2,12,11,10,9,8,7,6,5,4,3)]
colnames(pop90_99)<-c("FIPS","STATE",1990:1999)


pop200010<-read.csv("st-est00int-01.csv",header=TRUE,stringsAsFactors = FALSE)
pop2000_9<-pop200010[-c(1:8,60:69),-c(2,13,14)]
colnames(pop2000_9)<-c("STATE",c(2000:2009))
pop2000_9$FIPS<-GEOID$GEOID[2:52]
pop2000_9<-pop2000_9[,c(12,1:11)]

for (i in 3:12)
{pop2000_9[,i]<-as.numeric(gsub(",","",pop2000_9[,i]))}



pop1018<-read.csv("PEP_2018_PEPANNRES_with_ann.csv",header=TRUE,stringsAsFactors = FALSE)
pop10_18<-pop1018[-1,-c(1,4,5)]
colnames(pop10_18)<-c("FIPS","STATE",c(2010:2018))
for (i in 3:11)
{pop10_18[,i]<-as.numeric(pop10_18[,i])}

pop10_18<-pop10_18[-52,]



pop_timeseries<-cbind(pop50_60[,1:12],pop60_70[,3:12],pop70_80[,3:12],
                      pop80_90[,3:12],pop90_99[,3:12],pop2000_9[,3:12],pop10_18[,3:11])
pop_timeseries$STATE<-GEOID$STATE[2:52]


pop_timeseries<-t(pop_timeseries)

colnames(pop_timeseries)<-pop_timeseries[1,]
pop_timeseries<-pop_timeseries[-c(1,2),]
rownames(pop_timeseries)<-c(1:69)
pop_timeseries<-as.data.frame(pop_timeseries)
pop_timeseries$year<-c(1950:2018)
pop_timeseries<-pop_timeseries[,c(52,1:51)]

for (i in 2:52)
{pop_timeseries[,i]<-as.numeric(as.character(pop_timeseries[,i]))}


#------------------------- Population density------------------------

Area=read.csv("landarea.csv",header=TRUE,stringsAsFactors = FALSE)

pop<-array(NA,c(69,51))
for (i in 1:51)
{pop[,i]<-(pop_timeseries[1:69,i+1])/Area[i,2] }


pop2<-pop[,-(which(colMeans(count==0) > 0.6))]

#---- map all events 1950-2018 based on intensity of storm 

storm<-read.csv("1950-2018_all_events.csv",header=TRUE,stringsAsFactors = FALSE)
storm$magnitude[which(tornado$magnitude==(-9))]<-NA
F01<-filter(tornado,magnitude %in% c(0,1))
F2<-filter(tornado,magnitude %in% 2)
F3<-filter(tornado,magnitude %in% 3)
F45<-filter(tornado,magnitude %in% c(4,5))

clr<-as.matrix(c("green3","gold1","chocolate1","red3"))

map<-map("state", col="gray90", fill=TRUE,type = "l",add=FALSE)


points(x = F01$slon, y = F01$slat, pch=16,cex=0.6, lwd= 4, col = rang[1])
points(x = F2$slon, y = F2$slat, pch=16,cex=0.7, lwd= 4, col = rang[2])
points(x = F3$slon, y = F3$slat, pch=16,cex=0.9, lwd= 4, col = rang[3])
points(x = F45$slon, y = F45$slat, pch=16,cex=1, lwd= 4, col = rang[4])

legend("bottomright",as.character(c('F0-F1','F2','F3','F4-F5')), col = clr[1:4,1], pch = 16, cex=1.7 )


#------------------------ climate variables --------------------------------

#-----------------ENSO-Nino 

Nino34<- read.delim("inino5_2.dat",sep="",skip=8, header = FALSE)
colnames(Nino34)<-c("year","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
Nino34<- Nino34[95:163,]

min(Nino34[,2:13],na.rm=TRUE)
Nino34[Nino34==-999.9]<-NA


Nino34$annual<-rowMeans(Nino34[,2:13],na.rm=TRUE)

#------------------------------ AO

AO <- read.delim("AO.dat",sep="",skip=1, header = FALSE)
AO<-AO[1:69,]

min(AO[,2:13],na.rm=TRUE)


AO$annual<-rowMeans(AO[,2:13],na.rm=TRUE)

#----------------------------  NAO 
NAO <- read.delim("nao_3dp.dat",sep="",skip=1, header = FALSE)
NAO<-NAO[129:197,1:13]

min(NAO[,2:13],na.rm=TRUE)

NAO[NAO==-99.99]<-NA

NAO$annual<-rowMeans(NAO[,2:13],na.rm=TRUE)

#---------------------------- SOI

SOI <- read.delim("icpc_soi.dat",sep="",skip=8, header = FALSE)
SOI<-SOI[65:133,]

min(SOI[,2:13])

SOI$annual<-rowMeans(SOI[,2:13],na.rm=TRUE)

#-------------------------------PDO

Pdata <- read.csv("PDO.csv",header=FALSE,stringsAsFactors = FALSE)

Pdata<-Pdata[1:828,]


PDO<-as.data.frame(array(NA,c(69,13)))
PDO[,1]<-c(1950:2018)

for ( i in 0:68)
  
{ PDO [i+1,2:13]<- Pdata[(12*i+1):(12*i+12),2]
}

min(PDO[,2:13])
PDO$annual<-rowMeans(PDO[,2:13],na.rm=TRUE)

#---------------------------- AMO

AMO<- read.delim("iamo_ersst_ts_2.dat",sep="",skip=47, header = FALSE)
colnames(AMO)<-c("year","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
AMO<-AMO[97:165,]

min(AMO[,2:13],na.rm=TRUE)
AMO[AMO==-999.9]<-NA

AMO$annual<-rowMeans(AMO[,2:13],na.rm=TRUE)


#--All

climate_data<-cbind(AO[,14],Nino34[,14],SOI[,14],AMO[,14],NAO[,14],PDO[,14])
colnames(climate_data)<-c("AO", "ENSO-nino3.4", "ENSO-SOI", "AMO", "NAO", "PDO")

#------------------ Principle Component Analysis--------------

#------Robust PCA 

Count<-as.matrix(count2[1:69,])
xcent = sweep(Count,2,colMeans(Count))
rob <-rpca(xcent)

#------------- PCA explained variance barplot and scores plot
score<-as.data.frame(rob$L.svd$u)
score$year<-c(1950:2018)
score<-score[,c(28,1:27)]


percentvar =round(as.matrix((rob$L.svd$d)^2/sum((rob$L.svd$d)^2)),digits=2)

b=barplot(percentvar[1:6,1], ylab="Explained variance", ylim=c(0,1),
          names=c("PC1","PC2","PC3","PC4","PC5","PC6"),space=0.2,cex.names = 1.7, cex.lab=1.7,cex.axis=1.8)
text(x = b, y = percentvar[1:6,1], label = percentvar[1:6,1], pos = 3, cex = 1.6, col="tomato3")



plot(c(1950:2018),rob$L.svd$u[,1],pch=16,xlab='Year',ylab='Scores-PC1',col='black',cex=2,cex.lab=1.8,cex.axis=1.6)
lines(lowess(c(1950:2018),rob$L.svd$u[,1],f=1/9),lwd=4,col="hotpink3")

#--------------------------------- Correlation with poulation and climate variables ----------------

scores_matrix=rob$L.svd$u[,1:3]
colnames(scores_matrix)<-c("PC1","PC2","PC3")
M<-cor(scores_matrix,climate_data)

Pmat<-as.matrix(array(NA,c(3,6)))
colnames(Pmat)<-colnames(M)
for (i in 1:3){
  
  for (j in 1:6){
    Pmat[i,j]<-cor.test(scores_matrix[,i],climate_data[,j])$p.value
  }}


corrplot(M, method="square", p.mat=Pmat,mar=c(1,0.5,0.5,1),
         insig="label_sig", sig.level = 0.05,tl.cex=1,tl.col="black",cl.cex=1,tl.srt=45)

Correlation=as.data.frame(t(cor(score1,pop2)))
Correlation$GEOID<-GEOID$GEOID

#--------------------------- Hierarchical Model --------------

gc()
invisible(utils::memory.limit(32000))


ao<-climate_data[,1]
nino34<-climate_data[,2]
soi<-climate_data[,3]
amo<-climate_data[,4]
nao<-climate_data[,5]
pdo<-climate_data[,6]
cross<-nao*nino34




period<-array(0,c(69,1))
period[42:69,1]<-1
period<-period[,1]
DI<-period[1:69]

pop2<-log(pop2)
gdp<-log(capital$change)
Reg<-as.factor(Alley[,2])

nsites<-47
nyears<-69

Y<-count2[1:69,1:47]
rownames(Y)<-c(1:69)


###### Bayesian model ########



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
  
  
  
  
  
  }",file="Jan18.bug" )


jagsfit<-jags(data=jags.data, inits=jags.inits, jags.params, n.iter=10000, n.chains=3, model.file='Jan18.bug')


attach.jags(jagsfit)