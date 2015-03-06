###
### Simulate Data
###

beta.0=-100 #Intercept
beta.1=50 #Slope
W=matrix(c(.6,.3,.3),1,3) # Offset Matrix for seed trap size
y.no.seed=rep(0,10) #Sites with no seeds

ht.obs=matrix(0,1100,3) #storage matrix
ht.obs[,1]=c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),rep(11,100),rep(12,100),rep(13,100),rep(14,100),rep(15,100),rep(16,100),rep(17,100),rep(18,100),rep(19,100),rep(20,100)) #Indexing plants at each site

# Giving heights to each plant at each site. 
ht.obs[ht.obs[,1]==1,3]=rpois(10,10)
ht.obs[ht.obs[,1]==2,3]=rpois(10,20)
ht.obs[ht.obs[,1]==3,3]=rpois(10,30)
ht.obs[ht.obs[,1]==4,3]=rpois(10,40)
ht.obs[ht.obs[,1]==5,3]=rpois(10,50)
ht.obs[ht.obs[,1]==6,3]=rpois(10,60)
ht.obs[ht.obs[,1]==7,3]=rpois(10,70)
ht.obs[ht.obs[,1]==8,3]=rpois(10,80)
ht.obs[ht.obs[,1]==9,3]=rpois(10,90)
ht.obs[ht.obs[,1]==10,3]=rpois(10,100)

ht.obs[ht.obs[,1]==11,3]=rpois(100,110)
ht.obs[ht.obs[,1]==12,3]=rpois(100,120)
ht.obs[ht.obs[,1]==13,3]=rpois(100,130)
ht.obs[ht.obs[,1]==14,3]=rpois(100,140)
ht.obs[ht.obs[,1]==15,3]=rpois(100,150)
ht.obs[ht.obs[,1]==16,3]=rpois(100,160)
ht.obs[ht.obs[,1]==17,3]=rpois(100,170)
ht.obs[ht.obs[,1]==18,3]=rpois(100,180)
ht.obs[ht.obs[,1]==19,3]=rpois(100,190)
ht.obs[ht.obs[,1]==20,3]=rpois(100,200)

x = cbind(1,as.vector(by(ht.obs[,3]>100,ht.obs[,1],sum)))
y.seeds=beta.0+beta.1*x[,2] #Sites with seeds
y.seed=c(y.no.seed,y.seeds[11:20])%*%W #Simulated seed collection data


###
### Run MCMC #Kind of slow...
###

mm.gamma=function(mean,var){
  shape=mean^2/var^2
  rate=mean/var^2
  
  list(shape=shape,rate=rate)
}

mm.gamma3=mm.gamma(100,50)
a3=mm.gamma3$shape
a4=mm.gamma3$rate
#par(mfrow=c(1,1))
#plot(density(rgamma(10000,a3,a4)))

mm.gamma4=mm.gamma(100,35)
a1=mm.gamma4$shape
a2=mm.gamma4$rate
#par(mfrow=c(1,1))
#hist(rgamma(10000,a1,a2),freq=FALSE)
#hist(rlnorm(10000,log(200),1/1.2^2),freq=FALSE)

n.mcmc=100000
setwd("~/Documents/Willows")
source("change.pt.pois.regMCMC6.R")
output6=change.pt.pois.regMCMC6(n.mcmc=n.mcmc,X=ht.obs,y.seed=y.seed,opt.ht.strt=200,hshape1=log(200),hshape2=1/1.1^2,opt.ht.tune=100,beta=c(-1,1),beta.mns=c(0,0),beta.var=25,beta.tune=c(1,.5),N=4,N.tune=5,N.shape=log(200),N.rate=1/1.2^2,x.strt=10)

output6$p.value

output6$accept.N/n.mcmc
output6$accept.opt.ht/(n.mcmc)
output6$accept.beta/(n.mcmc)

mean(output6$x.save)


par(mfrow=c(2,2))
hist(output6$opt.ht.save,breaks=100,freq=FALSE)
lines(density(rlnorm(n.mcmc,log(200),1/1.1^2)))

hist(output6$N.save,breaks=100,freq=FALSE)
lines(density(rlnorm(1000,log(200),1/1.2^2)))

hist(output6$beta.save[,1],breaks=100,freq=FALSE)
lines(density(rnorm(1000,0,25)))

hist(output6$beta.save[,2],breaks=100,freq=FALSE)
lines(density(rnorm(1000,0,25)))

beta.var=.1
beta0=rnorm(1000,0,beta.var)
beta1=rnorm(1000,0,beta.var)
test=notExp(x%*%t(cbind(beta0,beta1)))
par(mfrow=c(1,1))
hist(test,breaks=100,freq=FALSE)

## it worked with opt.ht gamma informed priors mean 100 var 35 and a gamma proposal for opt.ht
    