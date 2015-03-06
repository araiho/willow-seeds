change.pt.pois.regMCMC.no.inflat1 <- function(n.mcmc,X,y.seed,opt.ht.strt,hshape1,hshape2,beta,beta.mns,beta.tune,beta.var,error.var,error.tune){
	
	###
	### Subroutines and Libraries 
	###
	
	library(mvtnorm)
	library(data.table)
	dunifdisc<-function(x, min=0, max=1) ifelse(x>=min & x<=max & round(x)==x, 1/(max-min+1), 0)
	runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)

	
	###
	### Storage
	###

	opt.ht.save=rep(0,n.mcmc)
	beta.save=matrix(0,n.mcmc,2)
	error.save=matrix(0,21,n.mcmc)
	
	mse.y=rep(0,n.mcmc)
	mse.ypred=rep(0,n.mcmc)
	msediffsave=rep(0,n.mcmc)
	ypred=matrix(0,length(y.seed[,1])*3,n.mcmc)
	msesave=rep(0,n.mcmc)

	
	###
	### Initial Conditions
	###
	
	opt.ht = opt.ht.strt
	error=rep(0,21)	
	error.mns=rep(0,21)
	W = matrix(c(.5,.25,.25),1,3)
	X1 = data.table(X)
    x = cbind(1,log(X1[,list(A=sum(V3>opt.ht),B=V1),by=V1][,A]))
   
    beta = beta
	
	accept.opt.ht = 0
	accept.beta = 0
	
	###
	### Begin MCMC
	###
		
	for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k," ");flush.console()
    
    ###    
    ### Sample Opt.ht
    ###
    if (TRUE){
    opt.ht.star = runifdisc(1,opt.ht-100,opt.ht+100)
    
    if(opt.ht.star>0 & opt.ht.star<=max(X[,3])){
    
    x.star = cbind(1,log(X1[,list(A=sum(V3>opt.ht.star),B=V1),by=V1][,A]))
    
 	  mh1.ht = sum(dpois(as.vector(y.seed[x.star[,2]>=0,]),as.vector(exp(x.star[x.star[,2]>=0,]%*%beta+error[x.star[,2]>=0])%*%W),log=TRUE)) + dunifdisc(opt.ht.star,hshape1,hshape2) + dunifdisc(opt.ht,opt.ht.star-100,opt.ht.star+100)
	  mh2.ht = sum(dpois(as.vector(y.seed[x.star[,2]>=0,]),as.vector(exp(x[x.star[,2]>=0,]%*%beta+error[x.star[,2]>=0])%*%W),log=TRUE)) + dunifdisc(opt.ht,hshape1,hshape2) + dunifdisc(opt.ht.star,opt.ht-100,opt.ht+100)
	  mh.ht = exp(mh1.ht - mh2.ht)
    
    if(mh.ht > runif(1)){
      opt.ht=opt.ht.star
      accept.opt.ht=accept.opt.ht+1
      x=x.star
    }
    }
    }
   
   
   
    ###
    ### Sample Betas    
    ###
    
    if (TRUE){
	  beta.star = rnorm(2,beta,beta.tune)

	  mh1.beta = sum(dpois(as.vector(y.seed[x[,2]>=0,]),as.vector(exp(x[x[,2]>=0,]%*%beta.star+error[x[,2]>=0])%*%W),log=TRUE))+dmvnorm(beta.star,beta.mns,beta.var*diag(2),log=TRUE)
	  mh2.beta = sum(dpois(as.vector(y.seed[x[,2]>=0,]),as.vector(exp(x[x[,2]>=0,]%*%beta+error[x[,2]>=0])%*%W),log=TRUE))+dmvnorm(beta,beta.mns,beta.var*diag(2),log=TRUE) #rep(beta.var,2)*diag(2)#solve(t(x[x[,2]>0,])%*%x[x[,2]>0,])
	  mh.beta = exp(mh1.beta-mh2.beta)

   if(mh.beta > runif(1)){
      beta = beta.star
      accept.beta = accept.beta+1
      
    }
    }
    
    ###
    ### Sample Error    
    ###
    
    if (TRUE){
	  error.star = rnorm(21,error,error.tune)
	  
	  if (min(exp(x[x[,2]>=0,]%*%beta.star+error.star[x[,2]>=0])%*%W)>=0){
	  	

	  mh1.error = dpois(y.seed[x[,2]>=0,],exp(x[x[,2]>=0,]%*%beta.star+error.star[x[,2]>=0])%*%W,log=TRUE)+dnorm(error.star[x[,2]>=0],error.mns[x[,2]>=0],error.var,log=TRUE)
	  mh2.error = dpois(y.seed[x[,2]>=0,],exp(x[x[,2]>=0,]%*%beta+error[x[,2]>=0])%*%W,log=TRUE)+dnorm(error[x[,2]>=0],error.mns[x[,2]>=0],error.var,log=TRUE) #rep(beta.var,2)*diag(2)#solve(t(x[x[,2]>0,])%*%x[x[,2]>0,])
	  mh.error = exp(mh1.error-mh2.error)

	  mh.index=which(mh.error > runif(length(mh.error[,1])*length(mh.error[1,])),arr.ind=TRUE)
  	  error[unique(mh.index[,1])] = error.star[unique(mh.index[,1])]
    }
    }
    
    ###
    ### Calculate Posterior Predictive Checks
    ###
    

	mn.nbinom=exp(x%*%beta+error)%*%W
    ypred[,k]=rpois(length(mn.nbinom),mn.nbinom)
    
    mse.y[k]=mean((as.vector(y.seed)-mn.nbinom)^2)
    mse.ypred[k]=mean((ypred-as.vector(mn.nbinom))^2)
   

    ###
    ### Save Values
    ###
        
	opt.ht.save[k] = opt.ht
	beta.save[k,] = beta    
	msediffsave[k]=mse.ypred[k]-mse.y[k]
	error.save[,k]=error
        
  }
    
    ###
    ### Create Plots
    ###
	n.burn=round(.2*n.mcmc)
	p.value=sum(mse.ypred[n.burn:n.mcmc] > mse.y[n.burn:n.mcmc]) / (n.mcmc - n.burn)
	
	
	layout(matrix(1:3,3,1))
	plot(opt.ht.save[n.burn:n.mcmc],type="l",main="Threshold Height")
	plot(beta.save[n.burn:n.mcmc,1],type="l",main=expression(paste(beta[0])))
	plot(beta.save[n.burn:n.mcmc,2],type="l",main=expression(paste(beta[1]))) 
	

    
    list(error.save=error.save,ypred=ypred,msediffsave=msediffsave,p.value=p.value,opt.ht.save=opt.ht.save,beta.save=beta.save,accept.beta=accept.beta,accept.opt.ht=accept.opt.ht)
}