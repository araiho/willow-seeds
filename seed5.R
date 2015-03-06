rm(list=ls())

setwd("/Users/raiho/Documents/Willows")
library(reshape)
library(segmented)
raw_seeds=read.csv("raw_seeds.csv")
melt.raw.seeds=melt(raw_seeds)#,id=c("Site","Trap.Type","X..Seeds"),na.rm=TRUE)
cast.raw.seeds=cast(melt.raw.seeds, Site ~ Trap.Type , sum)
#cast.raw.seeds1=cast(melt.raw.seeds, Site ~ Trap.Type , sum)
remove=c( 11, 12, 20) #,22 to remove WBT1
cast.raw.seeds.remove=cast.raw.seeds[-remove,]

remove.cols=c(3,4,5,6,7,8)
seed.count.traps=cast.raw.seeds.remove[,-remove.cols]
as.numeric(seed.count.traps[,1])

y.seed=as.matrix(seed.count.traps[,2:4])

seed.count=rowSums(cast.raw.seeds.remove[,2:10])

ht.obs=(read.csv("sex_obs_sites.csv"))
seed1=(read.csv("seedrain1.csv"))
ht.obs$sex<-NULL
ht.obs$species<-NULL
head(ht.obs)
ht.obs=ht.obs[ht.obs[,1]!="Rose",] #Takes Rose out because it has no heights
ht.obs[,1]=as.numeric(factor(ht.obs[,1])) #Giving index

ht.obs3=matrix(0,1,3)
ht.obs=matrix(rbind(as.matrix(ht.obs),ht.obs3),2229,3)

ht.obs1=matrix(0,90,3) #storage matrix
ht.obs1[,1]=c(rep(22,10),rep(23,10),rep(24,10),rep(25,10),rep(26,10),rep(27,10),rep(28,10),rep(29,10),rep(30,10)) #Indexing plants at each site

ht.obs=rbind(as.matrix(ht.obs),ht.obs1,names=NA)
ht.obs2=matrix(ht.obs,2319,3)
seed.count=c(seed.count,rep(0,9))

plot(sort(ht.obs[,3]))

ht.obs=cbind(ht.obs,rep(0,length(ht.obs[,1])))
for(i in 1:21){
ht.obs[ht.obs[,1]==i,4]=seed.count[i]
}
ht.obs=cbind(ht.obs,rep(0,length(ht.obs[,1])))

for(i in 1:length(ht.obs[,1])){
if(ht.obs[i,4]>0){
	ht.obs[i,5]=1
}
}

plot(ht.obs[ht.obs[,5]==0,3],seq(1,length(ht.obs[ht.obs[,5]==0,3]),1),col=2,xlim=c(0,450),ylim=c(0,2500))
points(ht.obs[ht.obs[,5]==1,3],seq(1,length(ht.obs[ht.obs[,5]==1,3]),1))

ht.obs=ht.obs[order(ht.obs[,4]),]
boxplot(ht.obs[,3]~ht.obs[,4]+ht.obs[,1])
for(i in 2:21){
points(ht.obs[ht.obs[,1]==i,3],col=i)
}

max.ht=round(max(ht.obs[,3])-1)
x.mat=matrix(0,length(unique(ht.obs[,1])),max.ht)
for (i in 1:max.ht){
	x.mat[,i] = as.vector(by(ht.obs[,3]>i,ht.obs[,1],sum))
}

seed.sim=matrix(0,length(unique(ht.obs[,1])),max.ht)
beta=c(-1,1)
W=matrix(c(.5,.25,.25),1,3)
for(i in 1:max.ht){
	seed.sim[,i]=rpois(20,exp(cbind(1,log(x.mat[,i]))%*%beta))
}
par(mfrow=c(1,2))
matplot(seed.sim,type="l",main="Simulated Seed Production")
lines(rowSums(y.seeds),lwd=3)

value=rep(0,length(unique(ht.obs[,1])))
for(i in 1:length(value)){
value[i] = which(seed.sim[i,] == y.seeds[i,])[1]
}
hist(value,breaks=21,col=8,xlim=c(0,max.ht),xlab=expression(paste(tilde(h))),main="beta=c(-1,1)")


matplot(x.mat,type="l")


par(mfrow=c(1,2))
matplot(seq.x,x.mat[1,],type="l",ylim=c(0,600),xlab=expression(paste(tilde(h))),ylab="x")
for(i in 2:21){
	matlines(seq.x,x.mat[i,],col=i)
}
matplot(seq.x,x.mat[1,],type="l",ylim=c(0,5),xlab=expression(paste(tilde(h))),ylab="x")
for(i in 2:21){
	matlines(seq.x,x.mat[i,],col=i)
}
par(mfrow=c(1,1))
seq.site=seq(1,21,1)
matplot(seq.site,x.mat[,1],type="l",xlab="Sites",ylab="x")
for(i in 2:440){
	matlines(seq.site,x.mat[,i],col=i)
}

###Putting in order by seeds observed

ht.obs1=cbind(ht.obs,matrix(0,2228,3))
for(r in 1:21){
	ht.obs1[ht.obs1[,1]==r,4:6]=seed.count.traps[r,2:4]	
}
plot(ht.obs1[,3],ht.obs1[,4])

ht.obs1=ht.obs1[order(ht.obs1[,6]),]

f=unique(ht.obs1[,1])


x.mat1=matrix(0,21,440)
seq.x1=seq(1,max(ht.obs1[,3]),1)
for (i in 1:440){
	x.mat1[,i] = as.vector(by(ht.obs1[,3]>i,ht.obs1[,1],sum))
}

x.mat3=x.mat1[order(seed.count),]

seed.count

par(mfrow=c(2,2))
seq.site=seq(1,21,1)
matplot(seq.site,x.mat3[,1],type="l",xlab="Sites",ylab="x")
for(i in 2:100){
	matlines(seq.site,x.mat3[,i],col=i)
}
matplot(seq.site,x.mat3[,101],type="l",xlab="Sites",ylab="x")
for(i in 102:200){
	matlines(seq.site,x.mat3[,i],col=i)
}
matplot(seq.site,x.mat3[,201],type="l",xlab="Sites",ylab="x")
for(i in 202:300){
	matlines(seq.site,x.mat3[,i],col=i)
}
matplot(seq.site,x.mat3[,301],type="l",xlab="Sites",ylab="x")
for(i in 302:440){
	matlines(seq.site,x.mat3[,i],col=i)
}

par(mfrow=c(2,2))
seq.site=seq(1,21,1)
matplot(sort(seed.count),x.mat3[,1],type="l",xlab="Seeds",ylab="x")
for(i in 2:100){
	matlines(sort(seed.count),x.mat3[,i],col=i)
}


par(mfrow=c(1,1))
seq.site=seq(1,21,1)
matplot(seq.site,x.mat3[,1],type="l",xlab="Sites",ylab="x")
for(i in 2:440){
	matlines(seq.site,x.mat3[,i],col=i)
}



abline(v=130)
c=rep(0,21)
c1=rep(0,21)
for(i in 1:21){
	c[i]=sum(ht.obs[ht.obs[,1]==i,3] > 200)
	#c1[i]=length(ht.obs[ht.obs[,1]==i,3])
}
plot(c1,seed.count)

ht.obs2=data.frame(matrix(0,1,4))
names(ht.obs2)=names(ht.obs1)

for(i in 1:21){
	if (sum(ht.obs1[ht.obs1[,1]==i,3] > 100) > 100 ){
	ht.obs2=rbind(ht.obs2,ht.obs1[ht.obs1[,1]==i,1:4])
	} 
	}
	
plot(ht.obs2[,3],ht.obs2[,4])
#find percentage above 2m
#find where the willows are mostly in prod sites or no

find.sum1<-function(x,i,ht){
	c(1,sum(x[x[,1]==i,3] > ht))
}
p=matrix(0,21,2)
for(j in 1:21){
	p[j,]=find.sum(ht.obs,j,130)
}

quantile(ht.obs1[ht.obs1[,4]==0,3],c(.05,.5,.975))
quantile(ht.obs1[ht.obs1[,4]>0,3],c(.05,.5,.975))

par(mfrow=c(1,1))
plot(density(ht.obs[ht.obs[,4]==0,3]),xlim=c(0,500),lwd=2,main=NA,xlab="Willow Height (cm)")
lines(density(ht.obs[ht.obs[,4]>0,3]),col="green",lwd=2)
abline(v=137,lty=2,lwd=2)
abline(v=85,lty=2,lwd=2,col="blue")
legend("topright",c("Sites without seeds","Sites with seeds","Average Height over all sites","Height threshold estimation"),lty=c(1,1,2,2),col=c(1,"green",1,"blue"),lwd=c(2,2,2,2))

ht.obs1[,4]==0 

max.ht=rep(0,21)
mean.ht=rep(0,21)
for(j in 1:21){
	max.ht[j]=max(ht.obs[ht.obs[,1]==j,3])
	mean.ht[j]=mean(ht.obs[ht.obs[,1]==j,3])
}
par(mfrow=c(1,2))
plot(max.ht,seed.count,xlab="Maximum Willow Height by Site (cm)",ylab="Seed Count")
plot(mean.ht,seed.count)

points(mean.ht,seed.count)

sapply(ht.obs,function(x) find.sum(ht.obs,))
i=as.list(seq(1,21,1))
aggregate(x=ht.obs,FUN=find.sum,by=i,simplify=FALSE)
aggregate(ht.obs,list(Height=ht.obs[,3] > 110,Site=ht.obs[,1]),length)

length(ht.obs,ht.obs[ht.obs==1,3]>100)



obs.melt=melt(ht.obs,id=c("site","sex_index"),na.rm=TRUE)
obs.ht=cast(obs.melt,site~variable,fun.aggregate=c(mean,var,sd))
obs.ht[,1]=as.numeric(obs.ht[,1])
obs.sd=signif(obs.ht[,4],3)
s.obs=cbind(seed1[,1],seed1[,3])

#plot(obs.ht,s.obs)

library(rjags)

obs.matrix=matrix(0,21,5)
obs.matrix=cbind(as.matrix(obs.ht),y.seed)
obs.matrix=obs.matrix[order(obs.matrix[,1]),] #It's very important they are in the right order!

seed.count=obs.matrix[,5]
obs.sd=obs.matrix[,4]
obs.ht=obs.matrix[,1:2]

data = list(obs.sd=obs.sd,obs.ht=obs.ht,seed.count=seed.count)

inits= list(list(beta1=13,beta0=110,alpha0=0,change.p=11,sigma=1),list(beta0=100,beta1=10,change.p=8,alpha0=2,sigma=5),list(beta0=120,beta1=18,change.p=16,alpha0=1,sigma=.5))
#alpha needs to start higher than 0.

n.adapt = 10000
n.update = 10000
n.iter =50000

jM=jags.model("/Users/raiho/Documents/Willows/seedJAGS6.R", data=data, n.chain=length(inits), inits=inits,n.adapt=n.adapt)

update(jM,n.iter=n.update)

zm=coda.samples(jM,variable.names=c("change.p","beta1","beta0","alpha0","sigma","p.value.mse"),n.iter=n.iter,n.thin=1)

zjm=jags.samples(jM,variable.names=c("mu"), n.iter=n.iter,n.thin=1)

zm$beta1

gelman.diag(zm)
summary(zm)
plot(zm)

plot(zjm$mse.y,zjm$mse.ypred)

hat=summary(zjm$mu,quantile,c(.025,.5,.975))$stat
#jpeg("/Users/raiho/Documents/Willows/seed.est2.jpeg")
#plot(sort(s.obs.list))
plot(obs.ht[,2],seed.count)
points(obs.ht[,2],hat[2,],col="red",typ="l")
points(obs.ht[,2],hat[3,],col="blue",typ="l")
points(obs.ht[,2],hat[1,],col="blue",typ="l")
#dev.off()

dim(zjm$s)
s=as.mcmc.list(zjm$s)

alpha=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
p=rdirichlet(1,alpha)
plot(p[1,])
rcat(1,p)


