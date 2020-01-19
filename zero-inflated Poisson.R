gc()
setwd("C:/Users/nilu/Desktop/Transfer to 117/climate indices")


#---------------- Climate predictors-------------


Nino12<- read.delim("inino2_2.dat",sep="",skip=8, header = FALSE)
colnames(Nino12)<-c("year","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
min(Nino12[,2:13],na.rm=TRUE)
Nino12[Nino12==-999.9]<-NA
Nino12$annual<-rowMeans(Nino12[,2:13],na.rm=TRUE)

######################### ENSO-  Nino 3 ###############################

Nino3<- read.delim("inino3_2.dat",sep="",skip=8, header = FALSE)
colnames(Nino3)<-c("year","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
min(Nino3[,2:13],na.rm=TRUE)
Nino3[Nino3==-999.9]<-NA
Nino3$annual<-rowMeans(Nino3[,2:13],na.rm=TRUE)


######################### Global Sea Surface Temperature Northern Hamisphere    ###################

SST <- read.delim("ihadcrut4_nh_2.dat",sep="",skip=10, header = FALSE)

colnames(SST)<-c("year","Month","Temp")
min(SST[,3],na.rm=TRUE)

SST_an<-data.frame(year=c(1850:2018),Temperature=array(NA,169))

for (j in 0:167)
{
  i=12*j
  
  SST_an$Temperature [j+1]<-mean(SST[(i+1):(i+12),3])
  
}
SST_an[169,2]<-mean(SST[2017:2023,3])



#######################   AO   ########################


AO <- read.delim("icpc_ao_2.dat",sep="",skip=5, header = FALSE)
colnames(AO)<-c("year","Month","Index AO")
min(AO[,3],na.rm=TRUE)


AO_an<-data.frame(year=c(1950:2018),Index=array(NA,69))

for (j in 0:67)
{
  i=12*j
  
  AO_an$Index [j+1]<-mean(AO[(i+1):(i+12),3])
  
}
AO_an[69,2]<-mean(AO[817:824,3])





################  indices ##################


new_F23<-read.csv("F23_tornado.csv",header=TRUE,stringsAsFactors = FALSE)

moderate<-new_F23[11:69,]
rownames(moderate)<-c(1:59)
no=length(which(moderate==0))

nonzero<-as.data.frame(which(moderate>0,arr.ind = TRUE))
zero<-as.data.frame(which(moderate==0,arr.ind = TRUE))

################   variables #################

Ind<-moderate
for (i in 1:59)
{ Ind[i,which(Ind[i,]>0)]<-1}

Y<-moderate


time<-c(1:59)
nino12<-Nino12[105:163,14]
nino3<-Nino3[105:163,14]
AO_an<-as.data.frame(AO_an)
AO<-AO_an[11:69,2]


Elevation<-El[,2]
Elevation<-Elevation*0.000189
Lat<-lat[,2]
Long<-long[,2]
Lat<-Lat*(pi/180)
Long<-Long*(pi/180)
Reg<-as.factor(Region[,2])






##---------------------- zero inflated Poisson and the Bayesian model----------------##

gc()
invisible(utils::memory.limit(32000))

set.seed(123)
jags.data <- list( 'Y' = Y,'Ind'=Ind, 'time'=time, 'nino12' = nino12,'nino3'=nino3,'AO'=AO, 
                   'Elevation'= Elevation,'Lat'= Lat,'Long'= Long,'Reg'=Reg,
                   'nonzero'=nonzero,'zero' = zero)
jags.params <- c('Y.pred','lambda','Ind.pred','P','alpha','beta','gama','teta',
                 'a0','a1','a2','a3','a4','c0','c1','c2','c3','c4',
                 'b01','b02','b03','b04','b11','b21','b31','b41','b12','b22','b32','b42',
                 'b13','b23','b33','b43','b14','b24','b34','b44',
                 'd01','d02','d03','d04','d11','d21','d31','d41','d12','d22','d32','d42',
                 'd13','d23','d33','d43','d14','d24','d34','d44',
                 'sigma_alpha','sigma_beta1','sigma_beta2','sigma_beta3','sigma_beta4',
                 'sigma_gama','sigma_teta1','sigma_teta2','sigma_teta3','sigma_teta4') 
jags.inits <- NULL




cat(
  "
  model  
  {          
  ### -------------------------------- Level 1 -------------------------------------### 
  #--------  Zero
  
  for (r in 1:167731)
  
  {
  
  
  Y[zero[r,1],zero[r,2]] ~ dpois(lambda[zero[r,1],zero[r,2]])
  lambda[zero[r,1],zero[r,2]] <- Ind[zero[r,1],zero[r,2]]
  Y.pred[zero[r,1],zero[r,2]] ~ dpois(lambda[zero[r,1],zero[r,2]])
  
  } 
  
  
  #------- Non-Zero
  
  for (r in 1:9800)
  {
  
  
  Y[nonzero[r,1],nonzero[r,2]] ~ dpois(lambda[nonzero[r,1],nonzero[r,2]])
  log(lambda[nonzero[r,1],nonzero[r,2]])  <- alpha[nonzero[r,2]]+(beta[nonzero[r,2],1]*time[nonzero[r,1]])+(beta[nonzero[r,2],2]*nino12[nonzero[r,1]])+(beta[nonzero[r,2],3]*nino3[nonzero[r,1]])+(beta[nonzero[r,2],4]*AO[nonzero[r,1]])
  Y.pred[nonzero[r,1],nonzero[r,2]] ~ dpois(lambda[nonzero[r,1],nonzero[r,2]])
  
  }
  
  
  
  
  
  
  
  
  
  ##------------------ Hierarchical Level for Covariance on alpha and beta ------------------ ##
  
  for (j in 1:3009)
  {
  
  alpha[j]~ dnorm(mu_alpha[j],sigma_alpha[j])
  mu_alpha[j]<- a0 + b01*Lat[j] + b02*Long[j] + b03*Elevation[j] + b04*Reg[j]
  
  beta[j,1] ~ dnorm(mu_beta1[j], sigma_beta1[j])
  mu_beta1[j] <- a1 + b11*Lat[j] + b21*Long[j] + b31*Elevation[j] + b41*Reg[j]
  
  beta[j,2] ~ dnorm(mu_beta2[j], sigma_beta2[j])
  mu_beta2[j] <- a2 + b12*Lat[j] + b22*Long[j] + b32*Elevation[j] + b42*Reg[j] 
  
  beta[j,3] ~ dnorm(mu_beta3[j], sigma_beta3[j])
  mu_beta3[j] <- a3 + b13*Lat[j] + b23*Long[j] + b33*Elevation[j] + b43*Reg[j] 
  
  beta[j,4] ~ dnorm(mu_beta4[j], sigma_beta4[j])
  mu_beta4[j] <- a4 + b14*Lat[j] + b24*Long[j] + b34*Elevation[j] + b44*Reg[j]
  
  
  sigma_alpha[j] ~ dunif(0,10)
  sigma_beta1[j] ~ dunif(0,10)
  sigma_beta2[j] ~ dunif(0,10)
  sigma_beta3[j] ~ dunif(0,10)
  sigma_beta4[j] ~ dunif(0,10)
  
  }
  
  
  ############## Priors for level 1 #####################



  a0 ~ dnorm(0,0.1)
  a1 ~ dnorm(0,0.1)
  a2 ~ dnorm(0,0.1)
  a3 ~ dnorm(0,0.1)
  a4 ~ dnorm(0,0.1)
  
  b01~ dnorm(0,0.1)
  b02~ dnorm(0,0.1)
  b03~ dnorm(0,0.1)
  b04~ dnorm(0,0.1)
  
  b11~ dnorm(0,0.1)
  b21~ dnorm(0,0.1)
  b31~ dnorm(0,0.1)
  b41~ dnorm(0,0.1)
  
  b12~ dnorm(0,0.1)
  b22~ dnorm(0,0.1)
  b32~ dnorm(0,0.1)
  b42~ dnorm(0,0.1)
  
  b13~ dnorm(0,0.1)
  b23~ dnorm(0,0.1)
  b33~ dnorm(0,0.1)
  b43~ dnorm(0,0.1)
  
  b14~ dnorm(0,0.1)
  b24~ dnorm(0,0.1)
  b34~ dnorm(0,0.1)
  b44~ dnorm(0,0.1)
  
  ### -------------------------------- Level 2 -------------------------------------###
  
  for (k in 1:3009)
  {
  for (m in 1:59)
  
  {
  Ind[m,k] ~ dbinom(P[m,k],1)
  logit(P[m,k]) <- gama[k]+(teta[k,1]*time[m])+(teta[k,2]*nino12[m])+ (teta[k,3]*nino3[m]) + (teta[k,4]*AO[m])
  Ind.pred[m,k] ~ dbern(P[m,k])
  }
  }
  
  
  ##------------------ Hierarchical Level for Covariance on gamma and tetas ------------------ ##
  
  for (n in 1:3009)
  {
  
  gama[n]~ dnorm(mu_gama[n],sigma_gama[n])
  mu_gama[n]  <- c0 + d01*Lat[n] + d02*Long[n] + d03*Elevation[n] + d04*Reg[n]
  
  teta[n,1] ~ dnorm(mu_teta1[n], sigma_teta1[n])
  mu_teta1[n] <- c1 + d11*Lat[n] + d21*Long[n] + d31*Elevation[n] + d41*Reg[n]
  
  teta[n,2] ~ dnorm(mu_teta2[n], sigma_teta2[n])
  mu_teta2[n] <- c2 + d12*Lat[n] + d22*Long[n] + d32*Elevation[n] + d42*Reg[n]
  
  teta[n,3] ~ dnorm(mu_teta3[n], sigma_teta3[n])
  mu_teta3[n] <- c3 + d13*Lat[n] + d23*Long[n] + d33*Elevation[n] + d43*Reg[n]
  
  teta[n,4] ~ dnorm(mu_teta4[n], sigma_teta4[n])
  mu_teta4[n] <- c4 + d14*Lat[n] + d24*Long[n] + d34*Elevation[n] + d44*Reg[n]
  
  
  
  sigma_gama[n] ~ dunif(0,10)
  sigma_teta1[n] ~ dunif(0,10)
  sigma_teta2[n] ~ dunif(0,10)
  sigma_teta3[n] ~ dunif(0,10)
  sigma_teta4[n] ~ dunif(0,10)
  }
  
  ################# Priors for level 2 #####################
  
  
  
  c0 ~ dnorm(0,0.01)
  c1 ~ dnorm(0,0.01)
  c2 ~ dnorm(0,0.01)
  c3 ~ dnorm(0,0.01)
  c4 ~ dnorm(0,0.01)
  
  d01~ dnorm(0,0.01)
  d02~ dnorm(0,0.01)
  d03~ dnorm(0,0.01)
  d04~ dnorm(0,0.01)
  
  d11~ dnorm(0,0.01)
  d21~ dnorm(0,0.01)
  d31~ dnorm(0,0.01)
  d41~ dnorm(0,0.01)
  
  d12~ dnorm(0,0.01)
  d22~ dnorm(0,0.01)
  d32~ dnorm(0,0.01)
  d42~ dnorm(0,0.01)
  
  d13~ dnorm(0,01)
  d23~ dnorm(0,01)
  d33~ dnorm(0,01)
  d43~ dnorm(0,01)
  
  d14~ dnorm(0,01)
  d24~ dnorm(0,01)
  d34~ dnorm(0,01)
  d44~ dnorm(0,01)
  
  
  }",file="trial2_severe.bug" )

gc()
jagsfit<-jags(data=jags.data, inits=jags.inits, jags.params, n.iter=1000,n.thin=3, n.chains=3, model.file='trial2_severe.bug')

gc()
#results
attach.jags(jagsfit)



